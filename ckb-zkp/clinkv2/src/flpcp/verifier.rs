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

use crate::flpcp::{Proof, VerMsg};

pub fn gen_vermsg<E: PairingEngine> (
    proof: Proof<E>,
    f_polys: &Vec<DensePolynomial<E::Fr>>,
    betas: &Vec<E::Fr>,
    r: E::Fr,
    M: usize,
) -> VerMsg<E> {
    let L6 = proof.ws.len();
    let f_r_shares = f_polys.iter().map(|f| f.evaluate(&r)).collect();
    let p_poly = DensePolynomial::from_coefficients_vec(proof.p_coeffs);
    assert!(r >= E::Fr::from(p_poly.len() as u64));
    println!("r: {:?}", r);
    let p_r_value = p_poly.evaluate(&r);
    let mut b = E::Fr::zero();
    // for i in 1..M+1 {
    //     let p_poly_i = p_poly.evaluate(&E::Fr::from(i as u32));
    //     // println!("p_poly_i: {:?}", p_poly_i);
    //     b += betas[i-1]*p_poly_i;
    // }
    let domain: GeneralEvaluationDomain<E::Fr> =
        EvaluationDomain::<E::Fr>::new(2*(M+1)).unwrap(); 
    let p_ploy_evals = p_poly.evaluate_over_domain_by_ref(domain);
    // println!("p_ploy_evals_over_domain: {:?}", p_ploy_evals);
    for i in 0..p_ploy_evals.evals.len() {
        b += p_ploy_evals.evals[i];
    }
    VerMsg {
        f_r_shares,
        p_r_value,
        b,
    }
}

pub fn verify_bgin19_proof<E: PairingEngine>(
    p_vermsg: VerMsg<E>,
    proof: Proof<E>,
    f_polys: &Vec<DensePolynomial<E::Fr>>,
    betas: &Vec<E::Fr>,
    thetas: &Vec<E::Fr>,
    r: E::Fr,
    M: usize,
) -> bool {
    let L = f_polys.len() / 6; 
    println!("L: {}", L); 
    let my_vermsg = gen_vermsg(proof, &f_polys, betas, r, M);
    // b = 0
    assert_eq!(my_vermsg.b + p_vermsg.b, E::Fr::zero());

    let mut f_r_values:Vec<E::Fr> = vec![];
    cfg_iter!(my_vermsg.f_r_shares)
        .zip(&p_vermsg.f_r_shares)
        .for_each(|(f_r_share1, f_r_share2)| f_r_values.push(*f_r_share1 + f_r_share2));
    let mut g_output_r = E::Fr::zero();
    for i in 0..L {
        // let id = 6 * i;
        let c_output_r = f_r_values[i]*f_r_values[i+2*L] 
                    + f_r_values[i]*f_r_values[i+3*L] 
                    + f_r_values[i]*f_r_values[i+2*L] 
                    + f_r_values[i] - f_r_values[i+5*L];        
        g_output_r +=  c_output_r * thetas[i];
    }
    println!("my_vermsg.p_r_value: {:?}", my_vermsg.p_r_value);
    println!("p_vermsg.p_r_value: {:?}", p_vermsg.p_r_value);
    assert_eq!(my_vermsg.p_r_value + p_vermsg.p_r_value, g_output_r);
    true
}
