//! Here we construct a polynomial commitment that enables users to commit to a
//! single polynomial `p`, and then later provide an evaluation proof that
//! convinces verifiers that a claimed value `v` is the true evaluation of `p`
//! at a chosen point `x`. Our construction follows the template of the construction
//! proposed by Kate, Zaverucha, and Goldberg ([KZG11](http://cacr.uwaterloo.ca/techreports/2010/cacr2010-10.pdf)).
//! This construction achieves extractability in the algebraic group model (AGM).

use ark_ec::{
    AffineCurve, PairingEngine, ProjectiveCurve,
};
use ark_ff::{Field, One, Zero};
use ark_poly::{polynomial::univariate::DensePolynomial, Polynomial};

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



fn accumulate_commitments_and_values<'a, E: PairingEngine>(
    _vk: &VerifierKey<E>,
    commitments: &[Commitment<E>],
    values: &[E::Fr],
    opening_challenge: E::Fr,
) -> Result<(E::G1Projective, E::Fr), Error> {
    //let acc_time = start_timer!(|| "Accumulating commitments and values");
    let mut combined_comm = E::G1Projective::zero();
    let mut combined_value = E::Fr::zero();
    let mut challenge_i = E::Fr::one();
    for (commitment, value) in commitments.into_iter().zip(values) {
        combined_comm += &commitment.0.mul(challenge_i);
        combined_value += &(*value * &challenge_i);
        challenge_i *= &opening_challenge.square();
    }

    //end_timer!(acc_time);
    Ok((combined_comm, combined_value))
}


pub fn batch_check<'a, E: PairingEngine>(
    vk: &VerifierKey<E>,
    commitments: &[Commitment<E>],
    point: E::Fr,
    values: &[E::Fr],
    proof: &Proof<E>,
    opening_challenge: E::Fr,
) -> Result<bool, Error> {
    //let check_time = start_timer!(|| "Checking evaluations");
    let (combined_comm, combined_value) =
        accumulate_commitments_and_values(vk, commitments, values, opening_challenge)?;
    let combined_comm = Commitment(combined_comm.into());
    let result = KZG10::<E, DensePolynomial<E::Fr>>::check(vk, &combined_comm, point, combined_value, proof)?;
    //end_timer!(check_time);
    Ok(result)
}


  