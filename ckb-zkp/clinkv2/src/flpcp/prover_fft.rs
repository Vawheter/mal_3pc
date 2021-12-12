#![allow(non_snake_case)]

use std::vec;

use ark_ec::PairingEngine;
use ark_ff::{One, UniformRand, Zero};
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain, UVPolynomial};
use ark_std::rand::RngCore;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::flpcp::Proof;

pub fn create_bgin19_proof<E: PairingEngine, R: RngCore>(
    inputs: Vec<Vec<Vec<E::Fr>>>,
    theta: E::Fr,
    ws: &Vec<Vec<E::Fr>>,
    rng: &mut R,
) -> (Proof<E>, Proof<E>) {
    
    // inputs[0], ..., inputs[5]
    assert_eq!(inputs.len(), 6);
    let L = inputs[0].len();
    let M = inputs[0][0].len();
    
    let zero = E::Fr::zero();
    let one = E::Fr::one();

    let domain: GeneralEvaluationDomain<E::Fr> =
        EvaluationDomain::<E::Fr>::new(M + 1).unwrap();
    let domain_size = domain.size();
    println!("domain size: {}", domain_size);

    // Compute p(x)
    // Using FFT - FFT: 6L*MlogM = 6mlogM; Poly Mul: 2L*(3*2Mlog2M) = 12mlog2M -> 18mlogM
    // Using Points - Evaluating M points: 6L*M*M = 6mM; Mul and iFFT: 3LM + 2Mlog2M -> 3m + 6mM + 2MlogM
    let mut f_polys: Vec<Vec<DensePolynomial<E::Fr>>> = vec![];

    // println!("inputs: {:?}", inputs);
    for k in 0..6 {
        let mut fk_polys = vec![];
        for l in 0..L {
            // Compute f(x) polynomial
            let points = [&[ws[k][l]], &inputs[k][l][..]].concat();
            let flk_coeffs = domain.ifft(&points);
            let flk_poly = DensePolynomial::from_coefficients_vec(flk_coeffs);
            // println!("f_poly.len(): {}", f_poly.coeffs.len());
            fk_polys.push(flk_poly);
        }
        f_polys.push(fk_polys);
    }
   
    // Compute p(x) polynomial, p(x) = theta1*(f1(x)*f3(x)+f1(x)*f4(x)+f2(x)*f3(x)+f5(x)-f6(x)) + ... 

    // let mut theta_power = one;
    // let mut p_poly = DensePolynomial::<E::Fr>::zero();
    // for l in 0..L {
    //     let mut c_poly = &f_polys[0][l] * &(&f_polys[2][l] + &f_polys[3][l]) + &f_polys[1][l] * &f_polys[2][l];
    //     c_poly += &f_polys[4][l];
    //     c_poly -= &f_polys[5][l];
    //     p_poly += (theta_power, &c_poly);
    //     theta_power *= theta;
    // }

    // Optimization
    let domain_2M: GeneralEvaluationDomain<E::Fr> =
        EvaluationDomain::<E::Fr>::new(2 * M - 1).unwrap();
    let domain_2M_size = domain_2M.size();
    let mut p_values = vec![zero; domain_2M_size];
    let mut theta_power = one;
    for l in 0..L {
        let fl_2M_values = (0..6).map(|k| 
            f_polys[k][l].evaluate_over_domain_by_ref(domain_2M).evals
        ).collect::<Vec<_>>();
        for j in 0..domain_2M_size {
            p_values[j] += theta_power * (fl_2M_values[0][j] * (fl_2M_values[2][j] + fl_2M_values[3][j]) + fl_2M_values[1][j] * fl_2M_values[2][j] + fl_2M_values[4][j] - fl_2M_values[5][j])
        }
        theta_power *= theta;
    }
    let p_coeffs = domain_2M.ifft(&p_values);

    // p(x) should be degree 2M , 2*(domain_size-1)+1 after padding
    let p_coeffs_share1: Vec<E::Fr> = (0..domain_2M_size).map(|_| E::Fr::rand(rng)).collect();
    let p_coeffs_share2: Vec<E::Fr> = (0..domain_2M_size).map(|i| p_coeffs[i] - p_coeffs_share1[i]).collect();

    let proof1 = Proof {
        p_coeffs_shares: p_coeffs_share1,
    };
    
    let proof2 = Proof {
        p_coeffs_shares: p_coeffs_share2,
    };

    (proof1, proof2)
}
