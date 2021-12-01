use std::vec;
use ark_ec::PairingEngine;
use ark_ff::{Field, One, ToBytes, UniformRand, Zero};
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain, Polynomial, UVPolynomial};
use ark_std::{cfg_iter, cfg_iter_mut};
use merlin::Transcript;
use rand::Rng;

// DEV
//use std::time::{Duration, Instant};
use ark_poly_commit::kzg10::KZG10;
pub type Kzg10VerKey<E: PairingEngine> = ark_poly_commit::kzg10::VerifierKey<E>;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::flpcp_opt::{Proof, VerMsg};

pub fn gen_vermsg<E: PairingEngine> (
    inputs: &Vec<Vec<E::Fr>>,
    r: E::Fr,
) -> VerMsg<E> {
    let m = inputs[0].len();

    let mut f_r_shares = vec![];
    let zero = E::Fr::zero();

    let domain: GeneralEvaluationDomain<E::Fr> =
        EvaluationDomain::<E::Fr>::new(m+1).unwrap();
    let lag_r_values = domain.evaluate_all_lagrange_coefficients(r);
    for i in 0..6 {
        let mut f_r_share = zero;
        for j in 0..m {
            f_r_share += &(lag_r_values[j] * &inputs[i][j]);
        }
        f_r_shares.push(f_r_share);
    }
    VerMsg {f_r_shares}
}

pub fn verify_bgin19_proof<E: PairingEngine>(
    kzg10_vk: &Kzg10VerKey<E>,
    p_vermsg: VerMsg<E>,
    proof: Proof<E>,
    inputs: &Vec<Vec<E::Fr>>,
    r: E::Fr,
) -> bool {
    let my_vermsg = gen_vermsg::<E>(inputs, r);
    let m = inputs[0].len();

    assert!(KZG10::<E, DensePolynomial<E::Fr>>::check(&kzg10_vk, &proof.q_comm, r, proof.q_r_value, &proof.q_r_proof).unwrap());

    let domain: GeneralEvaluationDomain<E::Fr> =
        EvaluationDomain::<E::Fr>::new(m+1).unwrap();
    let vanishing_poly = domain.vanishing_polynomial();
    let vanishing_value = vanishing_poly.evaluate(&r);
    

    let mut f_r_values:Vec<E::Fr> = vec![];
    cfg_iter!(my_vermsg.f_r_shares)
        .zip(&p_vermsg.f_r_shares)
        .for_each(|(f_r_share1, f_r_share2)| f_r_values.push(*f_r_share1 + f_r_share2));
    
    let p_r_value = f_r_values[0]*f_r_values[2] + f_r_values[0]*f_r_values[3] + f_r_values[1]*f_r_values[2] + f_r_values[4] - f_r_values[5];
    // assert_eq!(p_r_value, proof.q_r_value*vanishing_value);
    true
}
