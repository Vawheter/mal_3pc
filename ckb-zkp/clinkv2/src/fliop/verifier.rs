use ark_ec::PairingEngine;
use ark_ff::{Field, One, ToBytes, UniformRand, Zero};
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain, Polynomial, UVPolynomial};
use ark_std::{cfg_iter, cfg_iter_mut};
use merlin::Transcript;
use rand::Rng;

// DEV
//use std::time::{Duration, Instant};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::fliop::{Proof, VerMsg};

pub fn gen_vermsg<E: PairingEngine> (
    // ws: Vec<E::Fr>,
    proof: &Proof<E>,
    inputs: &Vec<Vec<E::Fr>>,
    gammas: &Vec<E::Fr>,
    lag_vals_domain_3: &Vec<Vec<E::Fr>>,
    lag_vals_domain_5: &Vec<Vec<E::Fr>>,
    ws: &Vec<E::Fr>,
    rs: &Vec<E::Fr>,
) -> VerMsg<E> {

    let mut c = inputs[5].iter().map(|zi| zi).sum();
    let mut f_r_values = inputs.clone();
    let L = proof.q_shares.len() - 1;
    let mut b = E::Fr::zero();
    for l in 0..L {
        let b_l = c - proof.q_shares[l][0] - proof.q_shares[l][1];
        b += gammas[l] * b_l;
        // q_r_share
        c = (0..3).map(|j| 
            lag_vals_domain_3[l][j] * proof.q_shares[l][j]
        ).sum();

        let len_by_2 = f_r_values[0].len() / 2;
        f_r_values = (0..5).map(|k|
            (0..len_by_2).map(|j|
                lag_vals_domain_3[l][0] * f_r_values[k][j] + lag_vals_domain_3[l][1] * f_r_values[k][j+len_by_2]
            ).collect::<Vec<_>>()
        ).collect::<Vec<Vec<_>>>();
    }
    c = (0..5).map(|j| 
        lag_vals_domain_5[L][j] * proof.q_shares[L][j]
    ).sum();
    let b_L = c - proof.q_shares[L][0] - proof.q_shares[L][1] - proof.q_shares[L][2] - proof.q_shares[L][3];
    b += gammas[L] * b_L;

    let f_r_values = (0..5).map(|k|
        lag_vals_domain_5[L][0] * ws[k] + lag_vals_domain_5[L][1] * f_r_values[k][0] + lag_vals_domain_5[L][2] * f_r_values[k][1]
    ).collect::<Vec<E::Fr>>();

    VerMsg {
        b_share: b,
        f_r_shares: f_r_values,
        q_r_share: c,
    }
}

pub fn verify_bgin19_proof<E: PairingEngine>(
    p_vermsg: VerMsg<E>,
    proof: Proof<E>,
    inputs: &Vec<Vec<E::Fr>>,
    gammas: &Vec<E::Fr>,
    lag_vals_domain_3: &Vec<Vec<E::Fr>>,
    lag_vals_domain_5: &Vec<Vec<E::Fr>>,
    ws: &Vec<E::Fr>,
    rs: &Vec<E::Fr>,
) -> bool {

    let my_vermsg = gen_vermsg::<E>(&proof, &inputs, &gammas, lag_vals_domain_3, lag_vals_domain_5, ws, rs);
    let b = p_vermsg.b_share + my_vermsg.b_share;
    // assert_eq!(b, E::Fr::zero());

    let q_r_value = p_vermsg.q_r_share + my_vermsg.q_r_share;
    let f_r_values = (0..5).map(|i| 
        p_vermsg.f_r_shares[i] + my_vermsg.f_r_shares[i]
    ).collect::<Vec<E::Fr>>();
    let f_r_value = f_r_values[0] * (f_r_values[2] + f_r_values[3]) + f_r_values[1] * f_r_values[2] + f_r_values[4]; 
    // assert_eq!(q_r_value, f_r_value);
   
    true
}
