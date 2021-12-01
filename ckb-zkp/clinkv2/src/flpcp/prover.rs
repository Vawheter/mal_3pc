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
    inputs: Vec<Vec<E::Fr>>,
    M: usize,
    L: usize,
    thetas: &Vec<E::Fr>,
    betas: &Vec<E::Fr>,
    rng: &mut R,
) -> (Proof<E>, Proof<E>, Vec<DensePolynomial<E::Fr>>) {
    
    // inputs[0], ..., inputs[5]
    let m = inputs[0].len();
    // let M = (m as f64).sqrt() as usize + 1 as;
    // let L = M;
    // println!("m_io: {:?}, m_mid: {:?}, n: {:?}", m_io, m_mid, n);
    let one = E::Fr::one();
    let zero = E::Fr::zero();

    // Check inputs, no problem
    // for j in 0..m {
    //     let r = 
    //         inputs[0][j]*inputs[2][j] + inputs[0][j]*inputs[3][j] + inputs[1][j]*inputs[2][j] + inputs[4][j] - inputs[5][j];
    //     assert_eq!(r, zero);
    // }

    assert!(m <= M * L);
    let max_degree = 2 * (M + 1);

    let domain: GeneralEvaluationDomain<E::Fr> =
        EvaluationDomain::<E::Fr>::new(M+1).unwrap();
    let domain_size = domain.size();
    println!("domain size: {}", domain_size);

    let mut f_polys: Vec<DensePolynomial<E::Fr>> = vec![];
    let mut ws: Vec<E::Fr> = vec![];
    let mut f_points:Vec<Vec<E::Fr>> = vec![];

    // println!("inputs: {:?}", inputs);
    for i in 0..6 {
        for j in 0..L {
            let mut points = vec![zero; M+1];
            // f0(x)
            let w = E::Fr::rand(rng);
            // let w = zero;
            points[0] = w; 
            // Split inputs[i] into L vectors of size M
            for k in 0..M {
                points[k+1] = inputs[i][j+k*L];
            }
            ws.push(w);
            // Compute f(x) polynomial
            // println!("i:{}, j:{}, points:{:?}", i, j, points);
            let f_coeffs = domain.ifft(&points);
            let f_poly = DensePolynomial::from_coefficients_vec(f_coeffs);
            // println!("f_poly.len(): {}", f_poly.coeffs.len());
            f_polys.push(f_poly);
            f_points.push(points);
        }
    }
    // let f0_ploy_evals = f_polys[0].evaluate_over_domain_by_ref(domain);
    // println!("f0_ploy_evals_over_domain: {:?}", f0_ploy_evals);
    // println!("f0_points: {:?}", f_points[0]);
    // let f0_ploy_eval_zero = f_polys[0].evaluate(&zero);
    // println!("f0_ploy_eval_at_zero: {:?}", f0_ploy_eval_zero);
    // let f0_ploy_eval_one = f_polys[0].evaluate(&one);
    // println!("f0_ploy_eval_at_one: {:?}", f0_ploy_eval_one);

    // Check points, no problem
    // for i in 0..L {
    //     for j in 1..M+1 {
    //         let r = 
    //             f_points[i][j]*f_points[i+2*L][j] + f_points[i][j]*f_points[i+3*L][j] + f_points[i+L][j]*f_points[i+2*L][j] + f_points[i+4*L][j] - f_points[i+5*L][j];
    //         assert_eq!(r, zero);
    //     }
    // }

    // println!("f_polys[0]: {:?}", f_polys[0]);
    // println!("f_polys[1]: {:?}", f_polys[1]);

    // let thetas:Vec<E::Fr> = (0..L).map(|_| E::Fr::rand(rng)).collect();
    // let betas:Vec<E::Fr> = (0..6*L).map(|_| E::Fr::rand(rng)).collect();
    
    // Compute p(x) polynomial
    // p(x) = theta1*(f1(x)*f3(x)+f1(x)*f4(x)+f2(x)*f3(x)+f5(x)-f6(x)) + ... 
    // let p_coeffs = vec![zero; domain_size];
    let mut p_poly = DensePolynomial::<E::Fr>::zero();
    for i in 0..L {
        // p(x) evaluates to zero for all x in [M], should work on coset
        // let mut f_coset_points:Vec<Vec<E::Fr>> = vec![];
        // for j in 0..6 {
        //     let coset_points = domain.coset_fft(&f_polys[i+j*L].coeffs);
        //     f_coset_points.push(coset_points);
        // }

        // let mut f1f3_coset_mul = domain.mul_polynomials_in_evaluation_domain(&f_coset_points[0], &f_coset_points[2]);
        // println!("f1f3_coset_mul.len(): {}", f1f3_coset_mul.len());
        // Check poly mul
        // domain.coset_ifft_in_place(&mut f1f3_coset_mul);

        // let f1f3_mul_poly = DensePolynomial::from_coefficients_vec(f1f3_coset_mul.clone());
        // let f1f3_mul_poly = &f_polys[i] * &f_polys[i+2*L];
        // let f1f3_mul_eta = f1f3_mul_poly.evaluate(&eta);
        // assert_eq!(f1_eta*f3_eta, f1f3_mul_eta);

        // let mut f1f4_coset_mul = domain.mul_polynomials_in_evaluation_domain(&f_coset_points[0], &f_coset_points[3]);
        // // domain.coset_ifft_in_place(&mut f1f4_coset_mul);
        // let mut f2f3_coset_mul = domain.mul_polynomials_in_evaluation_domain(&f_coset_points[1], &f_coset_points[2]);
        // // domain.coset_ifft_in_place(&mut f2f3_coset_mul);
        // let mut c_coset_outputs = vec![zero; domain_size];
        // for j in 0..domain_size {
        //     c_coset_outputs[j] = 
        //         f1f3_coset_mul[j] + f1f4_coset_mul[j] + f2f3_coset_mul[j] + f_polys[id+4][j] - f_polys[id+5][j];
        // }

        // Check points, no problem
        // for j in 1..M+1 {
        //     let r = f_points[i][j]*f_points[i+2*L][j] + f_points[i][j]*f_points[i+3*L][j] + f_points[i+L][j]*f_points[i+2*L][j] + f_points[i+4*L][j] - f_points[i+5*L][j];
        //     assert_eq!(r, zero);
        // }
        // println!("i: {}", i);
        // let fi_ploy_evals = f_polys[i].evaluate_over_domain_by_ref(domain);
        // println!("fi_ploy_evals_over_domain: {:?}", fi_ploy_evals);
        // println!("fi_points: {:?}", f_points[i]);
        // assert_eq!(f_points[i][0], f_polys[i].evaluate(&one));
        // // assert_eq!(f_points[i][1], f_polys[i].evaluate(&(one+one)));
        // assert_eq!(f_points[i][2], f_polys[i].evaluate(&(E::Fr::from(3 as u64))));

        // Check polys
        // let eta = one;
        // let mut fs_eta = vec![];
        // for j in 0..6 {
        //     fs_eta.push(f_polys[i+j*L].evaluate(&eta));
        // }
        // let r = 
        //     fs_eta[0]*fs_eta[2] + fs_eta[0]*fs_eta[3] + fs_eta[1]*fs_eta[2] + fs_eta[4] - fs_eta[5];
        // assert_eq!(r, zero);

        let f1f3_mul_poly = &f_polys[i] * &f_polys[i+2*L];
        // let f1f3_mul_eta = f1f3_mul_poly.evaluate(&eta);

        // assert_eq!(fs_eta[0]*fs_eta[2], f1f3_mul_eta);
        // let f1f3_mul_poly_evals = f1f3_mul_poly.evaluate_over_domain_by_ref(domain);
        // println!("f1f3_mul_poly_evals_over_domain: {:?}", f1f3_mul_poly_evals);
        // let mut f1f3_mul_points:Vec<E::Fr> = vec![];
        // for j in 0..M+1 {
        //     f1f3_mul_points.push(f_points[i][j]*f_points[i+2*L][j]);
        // }
        // println!("f1f3_mul_points: {:?}", f1f3_mul_points);

        let f1f4_mul_poly = &f_polys[i] * &f_polys[i+3*L];
        let f2f3_mul_poly = &f_polys[i+L] * &f_polys[i+2*L];

        let mut c_poly = DensePolynomial::<E::Fr>::zero();
        c_poly += &f1f3_mul_poly;
        c_poly += &f1f4_mul_poly; 
        c_poly += &f2f3_mul_poly;
        c_poly += &f_polys[i+4*L];
        c_poly -= &f_polys[i+5*L];
        // let p_poly = &c_poly * thetas[i];
        // let c_poly_evals = c_poly.evaluate_over_domain_by_ref(domain);
        // println!("c_ploy_evals_over_domain: {:?}", c_poly_evals);

        p_poly += (thetas[i], &c_poly);
        // for j in 0..domain_size {
        //     c_coset_outputs[j] = 
        //         f1f3_mul_poly[j] + f1f4_mul_poly[j] + f2f3_mul_poly[j] + f_polys[i+4*L][j] - f_polys[i+5*L][j];
        // }
        // domain.coset_ifft_in_place(&mut c_coset_outputs);
        // cfg_iter!(&mut p_coeffs)
        //     .zip(&c_coset_outputs)
        //     .for_each(|(p_coeff, c_coeff)| *p_coeff += (*c_coeff * thetas[i]))
        // if i == 1 {
        //     println!("c_coset_outputs: {:?}", c_coset_outputs);
        // }
        // for j in 0..domain_size {
        //     p_coeffs[j] += c_coset_outputs[j] * thetas[i];
        // }
        // let c_output_poly = DensePolynomial::from_coefficients_vec(c_coset_outputs);
        // println!("c_output_poly_2: {:?}", c_output_poly.evaluate(&(E::Fr::one()+E::Fr::one())));
    }
    // println!("p_poly_coeffs: {:?}", p_coeffs);
    // println!("p_poly.coeffs.len(): {}", p_poly.coeffs.len());
    // // let p_poly = DensePolynomial::from_coefficients_vec(p_coeffs);
    // let p_ploy_evals = p_poly.evaluate_over_domain_by_ref(domain);
    // println!("p_ploy_evals_over_domain: {:?}", p_ploy_evals);
    // for i in 1..M+1 {
    //     let p_poly_i = p_poly.evaluate(&E::Fr::from(i as u32));
    //     println!("p_poly_i: {:?}", p_poly_i);
    // } 

    let ws_share1: Vec<E::Fr> = (0..6*L).map(|_| E::Fr::rand(rng)).collect();
    let p_coeffs_share1: Vec<E::Fr> = (0..p_poly.coeffs.len()).map(|_| E::Fr::rand(rng)).collect();

    let ws_share2 = (0..6*L).map(|i| ws[i] - ws_share1[i]).collect();
    let p_coeffs_share2 = (0..p_poly.coeffs.len()).map(|i| p_poly.coeffs[i] - p_coeffs_share1[i]).collect();

    let proof1 = Proof {
        ws: ws_share1,
        p_coeffs: p_coeffs_share1,
    };
    
    let proof2 = Proof {
        ws: ws_share2,
        p_coeffs: p_coeffs_share2,
    };

    (proof1, proof2, f_polys)
}
