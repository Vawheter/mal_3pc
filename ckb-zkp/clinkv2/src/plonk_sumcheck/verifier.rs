#![allow(non_snake_case)]

use ark_ec::PairingEngine;
use ark_ff::One;
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};

use ark_poly_commit::kzg10::KZG10;
pub type Kzg10VerKey<E> = ark_poly_commit::kzg10::VerifierKey<E>;
use crate::batch_kzg10::batch_check;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::plonk_sumcheck::{Proof, VerMsg};

pub fn gen_vermsg<E: PairingEngine> (
    inputs: &Vec<Vec<Vec<E::Fr>>>,
    eta: E::Fr,
    r: E::Fr,
    z: E::Fr,
) -> VerMsg<E> {
    let L = inputs[0].len();
    let M = inputs[0][0].len();

    // Compute shares of p(z)
    let domain_M: GeneralEvaluationDomain<E::Fr> =
        EvaluationDomain::<E::Fr>::new(M).unwrap();

    let domain_L: GeneralEvaluationDomain<E::Fr> =
        EvaluationDomain::<E::Fr>::new(L).unwrap();

    let lag_r_M = domain_M.evaluate_all_lagrange_coefficients(r);
    let lag_z_L = domain_L.evaluate_all_lagrange_coefficients(z);

    let mut eta_power = E::Fr::one();
    let mut eta_powers = vec![];
    for _ in 0..M {
        eta_powers.push(eta_power);
        eta_power *= eta;
    }

    let mut f_r_shares = vec![];
    
    for k in 0..6 {
        if k == 2 || k == 3 {
            let fk_r_shares = (0..L).map(|l|
                (0..M).map(|j|
                    inputs[k][l][j] * lag_r_M[j]
                ).sum()
            ).collect::<Vec<E::Fr>>();
            f_r_shares.push(fk_r_shares);
        } else {
            let fk_r_shares = (0..L).map(|l|
                (0..M).map(|j|
                    eta_powers[j] * inputs[k][l][j] * lag_r_M[j]
                ).sum()
            ).collect::<Vec<E::Fr>>();
            f_r_shares.push(fk_r_shares);
        }
    }

    let F_z_shares = (0..6).map(|k|
        (0..L).map(|l|
            f_r_shares[k][l] * lag_z_L[l]
        ).sum()
    ).collect::<Vec<_>>();

    VerMsg {
        F_z_shares,
    }
}

pub fn verify_bgin19_proof<E: PairingEngine>(
    p_vermsg: VerMsg<E>,
    proof: Proof<E>,
    inputs: &Vec<Vec<Vec<E::Fr>>>,
    kzg10_vk: Kzg10VerKey<E>,
    eta: E::Fr,
    r: E::Fr,
    z: E::Fr,
) -> bool {
    let L = inputs[0].len();
    let M = inputs[0][0].len();

    // Proof Structure:
    //     poly_comms, // q(x), U(x), Q(x)
    //     open_values, // q(r), U(gz), U(z), Q(z)
    //     open_proofs, // q(r), U(gz), batch(U(z), Q(z))
    //     open_challenge 

    // Check poly openings: q(r)
    assert!(KZG10::<E, DensePolynomial<E::Fr>>::check(&kzg10_vk, &proof.poly_comms[0], r, proof.open_values[0], &proof.open_proofs[0]).unwrap());
    // Check poly openings: U(gz)
    let domain_L: GeneralEvaluationDomain<E::Fr> =
        EvaluationDomain::<E::Fr>::new(L).unwrap();
    let gz = z * domain_L.group_gen();
    assert!(KZG10::<E, DensePolynomial<E::Fr>>::check(&kzg10_vk, &proof.poly_comms[1], gz, proof.open_values[1], &proof.open_proofs[1]).unwrap());
    // Check poly openings: Q(z), U(z)
    assert!(batch_check(&kzg10_vk, &proof.poly_comms[1..], z, &proof.open_values[2..], &proof.open_proofs[2], proof.open_challenge).unwrap());

    let my_vermsg: VerMsg<E> = gen_vermsg(inputs, eta, r, z);
    // Reconstruct F(z)
    let F_z_values = (0..6).map(|k|
        my_vermsg.F_z_shares[k] + p_vermsg.F_z_shares[k]
    ).collect::<Vec<_>>();
    let domain_M: GeneralEvaluationDomain<E::Fr> =
        EvaluationDomain::<E::Fr>::new(M).unwrap();
    let v = proof.open_values[0] * domain_M.evaluate_vanishing_polynomial(r) * domain_L.size_inv();
    let p_z_value = F_z_values[0] * (F_z_values[2] + F_z_values[3]) + F_z_values[1] * F_z_values[2] + F_z_values[4] - F_z_values[5] - v;

    // Check U(gz) - U(z) = p(z) + Q(z)T(z)
    let T_z_value = domain_L.evaluate_vanishing_polynomial(z);
    assert_eq!(proof.open_values[1] - proof.open_values[2], p_z_value + proof.open_values[3] * T_z_value);
    
    true
}
