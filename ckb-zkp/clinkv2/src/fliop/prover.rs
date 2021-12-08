use std::slice::RSplit;
use std::vec;

use ark_ec::PairingEngine;
use ark_ff::{Field, One, ToBytes, UniformRand, Zero};
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain, Polynomial, UVPolynomial};
use ark_serialize::CanonicalSerialize;
use ark_std::{cfg_iter, cfg_iter_mut};
use rand::Rng;
use rand::prelude::*;
use ark_std::rand::RngCore;

// DEV
//use std::time::{Duration, Instant};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::fliop::Proof;

pub fn create_bgin19_proof<E: PairingEngine, R: RngCore>(
    inputs: Vec<Vec<E::Fr>>,
    rs: &Vec<E::Fr>,
    ws: &Vec<E::Fr>,
    rng: &mut R,
) -> (Proof<E>, Proof<E>) {
    
    // inputs[0]: x_i, x_i^{1,m}
    // inputs[1]: xi_1
    // inputs[2]: yi
    // inputs[3]: yi_1
    // inputs[4]: alpha_i
    // inputs[5]: zi
    let m = inputs[0].len();
    // Make sure that m is a power of 2
    assert_eq!(m, m.next_power_of_two());

    let zero = E::Fr::zero();
    let one = E::Fr::one();
    let two = one + one;
    let three = two + one;

    // Check inputs, no problem
    // for j in 0..m {
    //     let r = 
    //         inputs[0][j]*inputs[2][j] + inputs[0][j]*inputs[3][j] + inputs[1][j]*inputs[2][j] + inputs[4][j] - inputs[5][j];
    //     assert_eq!(r, zero);
    // }
    let max_degree = 4;

    let domain_2: GeneralEvaluationDomain<E::Fr> =
        EvaluationDomain::<E::Fr>::new(2).unwrap();
    let domain_2_size = domain_2.size();
    println!("domain size: {}", domain_2_size);

    // println!("inputs: {:?}", inputs);

    let mut L = m/2;
    let mut cnt = 0;

    let mut f_r_values = inputs.clone();
    let mut c = zero;

    let mut q_shares_1: Vec<Vec<E::Fr>> = vec![];
    let mut q_shares_2: Vec<Vec<E::Fr>> = vec![];
    // let mut c_shares_1: Vec<E::Fr> = vec![];
    // let mut c_shares_2: Vec<E::Fr> = vec![];

    loop {
        let mut q_L_poly = DensePolynomial::<E::Fr>::zero();
        let mut f_polys = vec![];
        for k in 0..5 {
            let mut fk_polys: Vec<DensePolynomial<E::Fr>> = vec![];
            for l in 0..L {
                // length of f_r_values: 2L (m -> 4)
                let points = [f_r_values[k][l], f_r_values[k][l+L]];
                let flk_coeffs = domain_2.ifft(&points);
                let flk_poly = DensePolynomial::from_coefficients_vec(flk_coeffs);
                fk_polys.push(flk_poly);
            }
            f_polys.push(fk_polys);
        }
        for l in 0..L {
            q_L_poly += &(&f_polys[0][l] * &(&f_polys[2][l] + &f_polys[3][l]) + &f_polys[1][l] * &f_polys[2][l]);
            q_L_poly += &f_polys[4][l];
        }
        // let r = E::Fr::rand(rng);
        // let q_r_value = q_L_poly.evaluate(&r);
        let mut q_L_values = vec![];
        q_L_values.push(q_L_poly.evaluate(&one));
        q_L_values.push(q_L_poly.evaluate(&two));
        q_L_values.push(q_L_poly.evaluate(&three));

        let q_L_shares_1 = (0..3).map(|_| E::Fr::rand(rng) ).collect::<Vec<E::Fr>>();
        let q_L_shares_2 = (0..3).map(|i| q_L_values[i] - q_L_shares_1[i] ).collect::<Vec<E::Fr>>();
        
        q_shares_1.push(q_L_shares_1);
        q_shares_2.push(q_L_shares_2);

        let r = rs[cnt];//domain_2.sample_element_outside_domain(rng);
        cnt += 1; 
        f_r_values = (0..5).map(|k|
            (0..L).map(|l|
                f_polys[k][l].evaluate(&r)
            ).collect::<Vec<_>>()
        ).collect::<Vec<Vec<_>>>();

        // c = q_L_poly.evaluate(&r);
        // let c_L_share_1 = E::Fr::rand(rng);
        // let c_L_share_2 = c - E::Fr::rand(rng);
        // c_shares_1.push(c_L_share_1);
        // c_shares_2.push(c_L_share_2);

        L /= 2;
        if L == 1 {
            let domain_4: GeneralEvaluationDomain<E::Fr> =
                EvaluationDomain::<E::Fr>::new(4).unwrap();
            assert_eq!(f_r_values[0].len(), 2);
            let mut q_L_poly = DensePolynomial::<E::Fr>::zero();
            // 5 polys
            let mut f_polys = vec![];
            for k in 0..5 {
                let mut points = vec![ws[k]; 1];
                points.append(&mut f_r_values[k]);
                let fk_coeffs = domain_4.ifft(&points);
                let fk_poly = DensePolynomial::from_coefficients_vec(fk_coeffs);
                f_polys.push(fk_poly);
            }

            q_L_poly += &(&f_polys[0] * &(&f_polys[2] + &f_polys[3]) + &f_polys[1] * &f_polys[2]);
            q_L_poly += &f_polys[4];

            let mut q_L_values = vec![];
            q_L_values.push(q_L_poly.evaluate(&one));
            q_L_values.push(q_L_poly.evaluate(&two));
            q_L_values.push(q_L_poly.evaluate(&three));
            q_L_values.push(q_L_poly.evaluate(&(two + two)));
            q_L_values.push(q_L_poly.evaluate(&(three + two)));

            let q_L_shares_1 = (0..5).map(|_| E::Fr::rand(rng) ).collect::<Vec<E::Fr>>();
            let q_L_shares_2 = (0..5).map(|i| q_L_values[i] - q_L_shares_1[i] ).collect::<Vec<E::Fr>>();
            
            q_shares_1.push(q_L_shares_1);
            q_shares_2.push(q_L_shares_2);

            // let r = domain_2.sample_element_outside_domain(rng);

            // c = q_L_poly.evaluate(&r);
            // let c_L_share_1 = E::Fr::rand(rng);
            // let c_L_share_2 = c - E::Fr::rand(rng);
            // c_shares_1.push(c_L_share_1);
            // c_shares_2.push(c_L_share_2);

            break;
        }
    }

    let proof1 = Proof {
        q_shares: q_shares_1,
        // c_shares: c_shares_1,
    };
    let proof2 = Proof {
        q_shares: q_shares_2,
        // c_shares: c_shares_2,
    };
    
    (proof1, proof2)
}
