#![allow(non_snake_case)]

use ark_ec::PairingEngine;
use ark_ff::{One, Zero};
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain, Polynomial, UVPolynomial};

// DEV
//use std::time::{Duration, Instant};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::flpcp::{Proof, VerMsg};

pub fn gen_vermsg<E: PairingEngine> (
    proof: Proof<E>,
    inputs: &Vec<Vec<Vec<E::Fr>>>,
    beta: E::Fr,
    ws: &Vec<Vec<E::Fr>>,
    r: E::Fr,
) -> VerMsg<E> {
    // let L6 = proof.ws.len();
    assert_eq!(inputs.len(), 6);
    let L = inputs[0].len();
    let M = inputs[0][0].len();
    
    // let one = E::Fr::one();
    // let zero = E::Fr::zero();

    let domain_M: GeneralEvaluationDomain<E::Fr> =
        EvaluationDomain::<E::Fr>::new(M + 1).unwrap();

    let lar_r_values = domain_M.evaluate_all_lagrange_coefficients(r);
    let mut f_r_shares = vec![];
    for k in 0..6 {
        let mut fk_r_shares = vec![];
        for l in 0..L {
            let mut fkl_r_shares = ws[k][l] * lar_r_values[0];
            for j in 0..M {
                fkl_r_shares += inputs[k][l][j] * lar_r_values[j + 1];
            } 
            fk_r_shares.push(fkl_r_shares);
        }
        f_r_shares.push(fk_r_shares);
    }

    let p_poly = DensePolynomial::from_coefficients_vec(proof.p_coeffs_shares);
    let p_r_share = p_poly.evaluate(&r);

    // p(x) degree: 2M, evaluates to zero over [1, M]
    // let domain_2M: GeneralEvaluationDomain<E::Fr> =
    //     EvaluationDomain::<E::Fr>::new(2 * M + 1).unwrap();
    // let p_ploy_evals = p_poly.evaluate_over_domain_by_ref(domain_2M);
    // println!("p_ploy_evals_over_domain: {:?}", p_ploy_evals);

    let mut b = E::Fr::zero();
    // IMPORTANT: should evaluate over domain-M since f(x) are calculated over domain-M and so is p(x)
    let g = domain_M.group_gen();
    // IMPORTANT: should begin with g^1 = g, since the first interpolated piont is random mask, p(x) evaluates to zero over domain [1, M]
    let mut g_power = g;
    let mut beta_power = E::Fr::one();
    for _ in 0..M {
        let p_eval_i = p_poly.evaluate(&g_power);
        b += beta_power * p_eval_i;
        beta_power *= beta;
        g_power *= g;
    }
    VerMsg {
        f_r_shares,
        p_r_share,
        b,
    }
}

pub fn verify_bgin19_proof<E: PairingEngine>(
    p_vermsg: &VerMsg<E>,
    proof: Proof<E>,
    inputs: &Vec<Vec<Vec<E::Fr>>>,
    beta: E::Fr,
    theta: E::Fr,
    ws: &Vec<Vec<E::Fr>>,
    r: E::Fr,
) -> bool {
    let L = inputs[0].len();

    let my_vermsg = gen_vermsg(proof, inputs, beta, ws, r);

    // b = 0
    assert_eq!(my_vermsg.b + p_vermsg.b, E::Fr::zero());

    let f_r_values = (0..6).map(|k| 
        (0..L).map(|l| 
            my_vermsg.f_r_shares[k][l] + p_vermsg.f_r_shares[k][l]
        ).collect::<Vec<_>>()
    ).collect::<Vec<_>>();
    
    let mut g_output_r = E::Fr::zero();
    let mut theta_power = E::Fr::one();
    for l in 0..L {
        g_output_r += theta_power * (f_r_values[0][l] * (f_r_values[2][l] + f_r_values[3][l]) + f_r_values[1][l] * f_r_values[2][l] + f_r_values[4][l] - f_r_values[5][l]);
        theta_power *= theta;
    }
    
    let p_r_value = my_vermsg.p_r_share + p_vermsg.p_r_share;
    // p(r) = g(r)
    assert_eq!(p_r_value, g_output_r);

    true
}
