//! Here we construct a polynomial commitment that enables users to commit to a
//! single polynomial `p`, and then later provide an evaluation proof that
//! convinces verifiers that a claimed value `v` is the true evaluation of `p`
//! at a chosen point `x`. Our construction follows the template of the construction
//! proposed by Kate, Zaverucha, and Goldberg ([KZG11](http://cacr.uwaterloo.ca/techreports/2010/cacr2010-10.pdf)).
//! This construction achieves extractability in the algebraic group model (AGM).

use ark_ec::{
    msm::{FixedBaseMSM, VariableBaseMSM},
    AffineCurve, PairingEngine, ProjectiveCurve,
};
use ark_ff::{to_bytes, Field, One, PrimeField, ToBytes, UniformRand, Zero};
use ark_poly::{polynomial::univariate::DensePolynomial, Polynomial, UVPolynomial};
use ark_serialize::*;
use ark_std::{cfg_iter, io, rand::RngCore};
use core::marker::PhantomData;
use core::ops::{Add, AddAssign};
use rand::Rng;

use ark_poly_commit::{
    Error, 
    PCRandomness,
    kzg10::{Powers, VerifierKey, Proof, Commitment, Randomness, KZG10}
};

fn check_degree_is_within_bounds(degree: usize, powers: usize) -> Result<(), Error> {
    if degree < 1 {
        Err(Error::DegreeIsZero)
    } else if degree > powers {
        Err(Error::PolynomialDegreeTooLarge { poly_degree: degree, supported_degree: powers, label: " ".to_string() })
    } else {
        Ok(())
    }
}

pub fn batch_open<'a, E: PairingEngine>(
    powers: &Powers<'_, E>,
    polynomials: &[DensePolynomial<E::Fr>],
    point: E::Fr,
    opening_challenge: E::Fr,
    rands: &Vec<Randomness<E::Fr, DensePolynomial<E::Fr>>>,
) -> Result<Proof<E>, Error> {
    let mut p = DensePolynomial::zero();
    let mut r = Randomness::empty();
    // let mut shifted_w = DensePolynomial::zero();
    // let mut shifted_r = Randomness::empty();
    // let mut shifted_r_witness = DensePolynomial::zero();

    let mut challenge_j = E::Fr::one();

    for (polynomial, rand) in polynomials.into_iter().zip(rands) {
        check_degree_is_within_bounds(polynomial.degree(), powers.size())?;
        // compute challenge^j and challenge^{j+1}.
        //let challenge_j = opening_challenge.pow([2 * j as u64]);
        p += (challenge_j, polynomial);
        r += (challenge_j, rand);

        challenge_j *= &opening_challenge.square();
    }

    //let proof_time = start_timer!(|| "Creating proof for unshifted polynomials");
    let proof = KZG10::open(powers, &p, point, &r)?;
    let w = proof.w.into_projective();
    let random_v = proof.random_v;
    //end_timer!(proof_time);

    Ok(Proof {
        w: w.into_affine(),
        random_v,
    })
}

/// Check that each `proof_i` in `proofs` is a valid proof of evaluation for
/// `commitment_i` at `point_i`.
pub fn batch_check_to_mul_values<E: PairingEngine, R: RngCore>(
    vk: &VerifierKey<E>,
    commitments: &[Commitment<E>],
    points: &[E::Fr],
    values: &[E::Fr],
    proofs: &[Proof<E>],
    rng: &mut R,
) -> Result<bool, Error> {
    // let check_time =
    //     start_timer!(|| format!("Checking {} evaluation proofs", commitments.len()));
    let g = vk.g.into_projective();
    let gamma_g = vk.gamma_g.into_projective();
    let mut total_c = <E::G1Projective>::zero();
    let mut total_w = <E::G1Projective>::zero();

    //let combination_time = start_timer!(|| "Combining commitments and proofs");
    let mut randomizer = E::Fr::one();
    // Instead of multiplying g and gamma_g in each turn, we simply accumulate
    // their coefficients and perform a final multiplication at the end.
    let mut g_multiplier = E::Fr::zero();
    let mut gamma_g_multiplier = E::Fr::zero();
    for (((c, z), v), proof) in commitments.iter().zip(points).zip(values).zip(proofs) {
        let w = proof.w;
        let mut temp = w.mul(*z);
        temp.add_assign_mixed(&c.0);
        let c = temp;
        g_multiplier += &(randomizer * v);
        if let Some(random_v) = proof.random_v {
            gamma_g_multiplier += &(randomizer * &random_v);
        }
        total_c += &c.mul(randomizer.into_repr());
        total_w += &w.mul(randomizer.into_repr());
        // We don't need to sample randomizers from the full field,
        // only from 128-bit strings.
        randomizer = u128::rand(rng).into();
    }
    total_c -= &g.mul(g_multiplier.into_repr());
    total_c -= &gamma_g.mul(gamma_g_multiplier.into_repr());
    //end_timer!(combination_time);

    //let to_affine_time = start_timer!(|| "Converting results to affine for pairing");
    let affine_points = E::G1Projective::batch_normalization_into_affine(&[-total_w, total_c]);
    let (total_w, total_c) = (affine_points[0], affine_points[1]);
    //end_timer!(to_affine_time);

    //let pairing_time = start_timer!(|| "Performing product of pairings");
    let result = E::product_of_pairings(&[
        (total_w.into(), vk.prepared_beta_h.clone()),
        (total_c.into(), vk.prepared_h.clone()),
    ])
    .is_one();
    //end_timer!(pairing_time);

    //end_timer!(check_time, || format!("Result: {}", result));
    Ok(result)
}


  