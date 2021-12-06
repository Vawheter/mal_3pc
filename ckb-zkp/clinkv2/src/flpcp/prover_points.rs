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
    lag_values: &Vec<Vec<E::Fr>>,
    rng: &mut R,
) -> (Proof<E>, Proof<E>) {
    
    // inputs[0], ..., inputs[5]
    assert_eq!(inputs.len(), 6);
    assert_eq!(inputs[0].len(), L);
    assert_eq!(inputs[0][0].len(), M);
    
    let one = E::Fr::one();
    let zero = E::Fr::zero();

    let domain: GeneralEvaluationDomain<E::Fr> =
        EvaluationDomain::<E::Fr>::new(2*M).unwrap();
    let domain_size = domain.size();
    println!("domain size: {}", domain_size);

    // Compute p(x)
    // Need M more points
    let theta = E::Fr::rand(rng);
    let mut theta_power = one;
    let mut p_points = vec![zero; 2 * M];
    for i in 0..M {
        let mut p_point = zero;
        for l in 0..L {
            let mut fl_pioints = vec![];
            for k in 0..6 {
                let flk_point: E::Fr = inputs[k][l].iter()
                    .zip(lag_values[i].iter()) // lagrange polynomial evaluations at M new points
                    .map(|(flk_value, lag_value)| *flk_value * lag_value)
                    .sum();
                fl_pioints.push(flk_point);
            }
            p_point += theta_power * (fl_pioints[0] * (fl_pioints[2] + fl_pioints[3]) + fl_pioints[1] * fl_pioints[2] + fl_pioints[4] - fl_pioints[5]);
        }
        theta_power *= theta;
        p_points[M+i] = p_point;
    }
    // println!("p_points: {:?}", p_points);
    let p_coeffs = domain.ifft(&p_points);
    let p_poly = DensePolynomial::from_coefficients_vec(p_coeffs);
    println!("p_poly: {:?}", p_poly);

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
