#![allow(non_snake_case)]

use ark_ec::PairingEngine;
use ark_ff::Zero;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};

// DEV
//use std::time::{Duration, Instant};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::fliop::{Proof, VerMsg};

pub fn gen_vermsg<E: PairingEngine> (
    proof: &Proof<E>,
    inputs: &Vec<Vec<E::Fr>>,
    gammas: &Vec<E::Fr>,
    thetas: &Vec<E::Fr>,
    ws: &Vec<E::Fr>,
    rs: &Vec<E::Fr>,
) -> VerMsg<E> {
    let m = inputs[0].len();
    // Make sure that m is a power of 2
    assert_eq!(m, m.next_power_of_two());

    let L = proof.q_shares.len() - 1;
    let mut b = E::Fr::zero();

    let domain_2: GeneralEvaluationDomain<E::Fr> =
        EvaluationDomain::<E::Fr>::new(2).unwrap();

    // evaluations of q(x) are over domain-szie 4
    let domain_4: GeneralEvaluationDomain<E::Fr> =
        EvaluationDomain::<E::Fr>::new(4).unwrap();

    let domain_8: GeneralEvaluationDomain<E::Fr> =
    EvaluationDomain::<E::Fr>::new(8).unwrap();

    let mut f_r_values = inputs.clone();
    for k in 0..2 {
        for j in 0..m {
            f_r_values[k][j] *= thetas[j];
        }
    }
    for k in 4..6 {
        for j in 0..m {
            f_r_values[k][j] *= thetas[j];
        }
    }

    let mut c: E::Fr = f_r_values[5].iter().map(|zi| zi).sum();

    // Common rounds
    for l in 0..L {
        let b_l = c - proof.q_shares[l][4] - proof.q_shares[l][5];
        b += gammas[l] * b_l;

        let lag_vals_domain_4 = domain_4.evaluate_all_lagrange_coefficients(rs[l]);

        // q_r_share
        c = (0..4).map(|j| 
            lag_vals_domain_4[j] * proof.q_shares[l][j]
        ).sum();

        let len_by_2 = f_r_values[0].len() / 2;
        let lag_vals_domain_2 = domain_2.evaluate_all_lagrange_coefficients(rs[l]);
        f_r_values = (0..5).map(|k|
            (0..len_by_2).map(|j|
                lag_vals_domain_2[0] * f_r_values[k][j] + lag_vals_domain_2[1] * f_r_values[k][j+len_by_2]
            ).collect::<Vec<_>>()
        ).collect::<Vec<Vec<_>>>();
    }

    // The last round
    // IMPORTANT: the last round [c] should still be equal to [q(1)] + [q(2)], which are indexed by 8, 9 respectively
    let b_L = c - proof.q_shares[L][8] - proof.q_shares[L][9];
    b += gammas[L] * b_L;

    let lag_vals_domain_8 = domain_8.evaluate_all_lagrange_coefficients(rs[L]);
    c = (0..8).map(|j| 
        lag_vals_domain_8[j] * proof.q_shares[L][j]
    ).sum();

    let lag_vals_domain_4 = domain_4.evaluate_all_lagrange_coefficients(rs[L]);
    let f_r_values = (0..5).map(|k|
        lag_vals_domain_4[0] * ws[k] + lag_vals_domain_4[1] * f_r_values[k][0] + lag_vals_domain_4[2] * f_r_values[k][1]
    ).collect::<Vec<E::Fr>>();

    VerMsg {
        b_share: b,
        f_r_shares: f_r_values,
        q_r_share: c,
    }
}

pub fn verify_bgin19_proof<E: PairingEngine>(
    p_vermsg: VerMsg<E>,
    proof: &Proof<E>,
    inputs: &Vec<Vec<E::Fr>>,
    gammas: &Vec<E::Fr>,
    thetas: &Vec<E::Fr>,
    ws: &Vec<E::Fr>,
    rs: &Vec<E::Fr>,
) -> bool {
    
    let my_vermsg = gen_vermsg::<E>(&proof, &inputs, &gammas, &thetas, ws, rs);

    let b = p_vermsg.b_share + my_vermsg.b_share;
    assert_eq!(b, E::Fr::zero());

    let q_r_value = p_vermsg.q_r_share + my_vermsg.q_r_share;
    let f_r_values = (0..5).map(|i| 
        p_vermsg.f_r_shares[i] + my_vermsg.f_r_shares[i]
    ).collect::<Vec<E::Fr>>();
    let f_r_value = f_r_values[0] * (f_r_values[2] + f_r_values[3]) + f_r_values[1] * f_r_values[2] + f_r_values[4]; 
    assert_eq!(q_r_value, f_r_value);
   
    true
}
