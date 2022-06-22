#![allow(non_snake_case)]

use std::vec;

use ark_ec::PairingEngine;
use ark_ff::{to_bytes, One, UniformRand, Zero, Field, PrimeField};
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain, Polynomial, UVPolynomial, domain};
use ark_std::rand::RngCore;
use merlin::Transcript;

use crate::batch_kzg10::batch_open;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::bbc19_fliop::{Proof, Kzg10ComKey};
use ark_poly_commit::kzg10::KZG10;

// #[inline]
// fn dot_prod_subcircuit_P1<E: PairingEngine>(
//     eta: E::Fr,
//     inputs: Vec<E::Fr>,
// ) -> E::Fr {   
//     // eta * (delta_x * (delta_y + r_y1) + delta_y * r_x1 + rxy1 - rz1)
//     eta * (inputs[0] + (inputs[2] + inputs[3]) + inputs[2] * inputs[1] + inputs[4] - inputs[5])
// }

pub fn prove_dot_prod<E: PairingEngine, R: RngCore>(
    inputs: &Vec<Vec<E::Fr>>,
    rhos: &Vec<E::Fr>, // randomness for hiding f(X)
    k: usize, // compression parameter
    sid: usize, // session id
    rng: &mut R,
) -> Proof<E> {
    // Inputs, T * (4n+3)
    let T: usize = inputs[0].len(); // number of dot products
    let L: usize = inputs.len(); // number of variables, (should be 4n+3, n for number of multiplications in one dot product)
    let n: usize = (L - 3) / 4;
    assert_eq!(4 * n + 3, L);
    assert_eq!(rhos.len(), L);

    let domain_k: GeneralEvaluationDomain<E::Fr> =
        EvaluationDomain::<E::Fr>::new(k).unwrap();
    // println!("domain_k: {:?}", domain_k);
    // println!("domain_k size: {}", domain_size);

    let one = E::Fr::one();
    let zero = E::Fr::zero();

    let S: usize = (T as f64 / k as f64).ceil() as usize;
    let mut s = S;
    let mut s0 = T;
    // println!("T: {:?}, k = {:?}", T, k);

    let mut transcript = Transcript::new(b"DZKP_DotProd");
    // Add public information
    transcript.append_message(b"k", &to_bytes!(k as u64).unwrap());
    transcript.append_message(b"sid", &to_bytes!(sid as u64).unwrap());
    
    // Inputs to circuits in iterations
    let mut f_evals_r: Vec<Vec<E::Fr>> = inputs.clone();
    let mut f_evals_r0: Vec<Vec<E::Fr>> = vec![];

    let mut buf = [0u8; 31];

    // Randomness for linear combination
    let mut eta = zero;
    // Powers of eta
    let mut eta_powers: Vec<E::Fr> = vec![];
    let mut f_eta_evals_r: Vec<E::Fr> = vec![];
    let mut f_eta_evals_r0 =  vec![];
    let mut f_eta_polys: Vec<DensePolynomial<E::Fr>> = vec![];

    // let eta = E::Fr::rand(rng);
    transcript.challenge_bytes(b"eta", &mut buf);
    eta = <E::Fr>::from_random_bytes(&buf).unwrap();
    let mut eta_power = one;
    
    for _ in 0..s {
        eta_powers.push(eta_power);
        eta_power *= eta;
    }

    f_evals_r = inputs.clone();
    f_eta_evals_r = eta_powers.clone();

    let mut G_eval_avg = zero;

    let mut p_coeffs_shares1: Vec<Vec<E::Fr>> = vec![];
    let mut p_coeffs_shares2: Vec<Vec<E::Fr>> = vec![];
    let mut G_evals = vec![];

    let mut r = zero;
    // For Fiat-Shamir
    let mut q_eval_r = zero;

    // println!("Begin Iteration...");
    let mut cnt = 0;
    loop {
        // println!("Iteration: {:?}", cnt);
        // cnt += 1;
        // println!("s0: {:?}", s0);
        // println!("s: {:?}", s);
        // println!("L: {:?}", L);

        // f(X), polynomials for input variables, s groups in total
        let mut f_polys: Vec<Vec<DensePolynomial<E::Fr>>> = vec![];
        let mut qt_eval_r = zero;

        for l in 0..L {
            let mut f_l_polys: Vec<DensePolynomial<E::Fr>> = vec![];
            for i in 0..s {
                let cur = i * k;
                let end = if i == s - 1 {s0} else {cur + k}; // automatic padding in FFT
                
                let coeffs = domain_k.ifft(&f_evals_r[l][cur..end]);
                let mut poly = DensePolynomial::from_coefficients_vec(coeffs);
                f_l_polys.push(poly);
            }
            f_polys.push(f_l_polys);
        }

            // let f_eta_poly = domain_k.ifft(&eta_powers);
            // let f_eta_cos_evals = domain_k.coset_fft(&f_eta_poly);

            let domain2k: GeneralEvaluationDomain<E::Fr> =
                EvaluationDomain::<E::Fr>::new(2 * k).unwrap();
            let domain2k_size = domain2k.size();
            let mut p_2k_evals = vec![zero; domain2k_size];
            
            let f_2k_evals = (0..L).map(|l| 
                (0..s).map(|i|
                    f_polys[l][i].evaluate_over_domain_by_ref(domain2k).evals
                ).collect::<Vec<_>>()
            ).collect::<Vec<_>>();
            
            for t in 0..domain2k_size {
                let mut p_t_eval = zero;
                for i in 0..s {
                    let mut ip_res = zero;
                    for j in 0..n {
                        let cur = 4 * j;
                        ip_res += f_2k_evals[cur][i][t] * (f_2k_evals[cur + 2][i][t] + f_2k_evals[cur + 3][i][t]) + f_2k_evals[cur + 2][i][t] * f_2k_evals[cur + 1][i][t];
                    }
                    ip_res += f_2k_evals[L - 3][i][t] + f_2k_evals[L - 2][i][t] - f_2k_evals[L - 1][i][t];
                    p_t_eval += f_eta_evals_r[i] * ip_res;
                }
                p_2k_evals.push(p_t_eval);
            }
            let p_coeffs = domain2k.ifft(&p_2k_evals);
            let p_coeffs_share1: Vec<E::Fr> = (0..domain2k_size).map(|_| E::Fr::rand(rng)).collect();
            let p_coeffs_share2: Vec<E::Fr> = (0..domain2k_size).map(|i| p_coeffs[i] - p_coeffs_share1[i]).collect();

            p_coeffs_shares1.push(p_coeffs_share1);
            p_coeffs_shares2.push(p_coeffs_share2);
              
            // Opening at random point r
            // let r = E::Fr::rand(rng);
            transcript.challenge_bytes(b"r", &mut buf);
            // Todo: Check r outside domain
            r = <E::Fr>::from_random_bytes(&buf).unwrap();
            
            let p_eval_r = DensePolynomial::from_coefficients_vec(p_coeffs).evaluate(&r);
            G_evals.push(p_eval_r);

            if s == 1 {
                break;
            }

        
        // Compute new f(r)
        f_evals_r = vec![];
        for l in 0..L {
            f_evals_r.push((0..s).map(|i| 
                f_polys[l][i].evaluate(&r)
            ).collect());
        }

        s0 = s;
        s = (s0 as f64 / k as f64).ceil() as usize;
   
        // Compute f_eta(r)
        f_eta_evals_r0 = f_eta_evals_r.clone();
        for i in 0..s {
            let cur = i * k;
            let end = if i == s - 1 {s0} else {cur + k}; // automatic padding in FFT
            
            let coeffs = domain_k.ifft(&f_eta_evals_r0[cur..end]);
            let mut poly = DensePolynomial::from_coefficients_vec(coeffs);
            f_eta_polys.push(poly);
        }

        f_eta_evals_r = vec![];
        for i in 0..s {
            f_eta_evals_r.push(f_eta_polys[i].evaluate(&r));
        }
        
    }

    // Out of loop
    Proof { 
        p_coeffs_shares1,
        p_coeffs_shares2,
        G_evals,
    }
}
