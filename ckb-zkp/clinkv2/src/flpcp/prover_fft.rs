use std::ops::Neg;

use ark_ec::PairingEngine;
use ark_ff::{Field, One, ToBytes, UniformRand, Zero};
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain, Polynomial, UVPolynomial};
use ark_serialize::CanonicalSerialize;
use ark_std::{cfg_iter, cfg_iter_mut};
use rand::Rng;
use ark_std::rand::RngCore;
use rand::prelude::*;

// DEV
//use std::time::{Duration, Instant};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::flpcp::Proof;

pub fn create_bgin19_proof<E: PairingEngine, R: RngCore>(
    inputs: Vec<Vec<Vec<E::Fr>>>,
    M: usize,
    L: usize,
    // theta: E::Fr,
    rng: &mut R,
) -> (Proof<E>, Proof<E>) {
    
    // inputs[0], ..., inputs[5]
    assert_eq!(inputs.len(), 6);
    assert_eq!(inputs[0].len(), L);
    assert_eq!(inputs[0][0].len(), M);
    
    let one = E::Fr::one();
    let zero = E::Fr::zero();

    let domain: GeneralEvaluationDomain<E::Fr> =
        EvaluationDomain::<E::Fr>::new(M).unwrap();
    let domain_size = domain.size();
    println!("domain size: {}", domain_size);

    // Compute p(x)
    // Using FFT - FFT: 6L*MlogM = 6mlogM; Poly Mul: 3L*(MlogM) = 3mlogM -> 9mlogM
    // Using Points - Evaluating M points: 6L*M*M = 6mM; Mul and iFFT: 3L + 2Mlog2M
    let mut f_polys: Vec<Vec<DensePolynomial<E::Fr>>> = vec![];

    // println!("inputs: {:?}", inputs);
    for k in 0..6 {
        let mut fk_polys = vec![];
        for l in 0..L {
            // Compute f(x) polynomial
            let flk_coeffs = domain.ifft(&inputs[k][l]);
            let flk_poly = DensePolynomial::from_coefficients_vec(flk_coeffs);
            // println!("f_poly.len(): {}", f_poly.coeffs.len());
            fk_polys.push(flk_poly);
        }
        f_polys.push(fk_polys);
    }
   
    // Compute p(x) polynomial
    // p(x) = theta1*(f1(x)*f3(x)+f1(x)*f4(x)+f2(x)*f3(x)+f5(x)-f6(x)) + ... 
    let theta = E::Fr::rand(rng);
    let mut theta_power = one;
    let mut p_poly = DensePolynomial::<E::Fr>::zero();
    for l in 0..L {
        let mut c_poly = &f_polys[0][l] * &(&f_polys[2][l] + &f_polys[3][l]) + &f_polys[1][l] * &f_polys[2][l];
        c_poly += &f_polys[4][l];
        c_poly -= &f_polys[5][l];
        p_poly += (theta_power, &c_poly);
        theta_power *= theta;
    }

    // let ws_share1: Vec<E::Fr> = (0..6*L).map(|_| E::Fr::rand(rng)).collect();
    let p_coeffs_share1: Vec<E::Fr> = (0..domain_size).map(|_| E::Fr::rand(rng)).collect();

    // let ws_share2 = (0..6*L).map(|i| ws[i] - ws_share1[i]).collect();
    let p_coeffs_share2 = (0..domain_size).map(|i| p_poly.coeffs[i] - p_coeffs_share1[i]).collect();

    let proof1 = Proof {
        // ws: ws_share1,
        p_coeffs_shares: p_coeffs_share1,
    };
    
    let proof2 = Proof {
        // ws: ws_share2,
        p_coeffs_shares: p_coeffs_share2,
    };

    (proof1, proof2)
}
