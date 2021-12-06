// use ark_ec::PairingEngine;
// use ark_ff::{Field, One, ToBytes, UniformRand, Zero};
// use ark_poly::polynomial::univariate::DensePolynomial;
// use ark_poly::{EvaluationDomain, GeneralEvaluationDomain, Polynomial, UVPolynomial};
// use ark_serialize::CanonicalSerialize;
// use ark_std::{cfg_iter, cfg_iter_mut};
// use rand::Rng;
// use ark_std::rand::RngCore;
// use rand::prelude::*;

// // DEV
// //use std::time::{Duration, Instant};

// #[cfg(feature = "parallel")]
// use rayon::prelude::*;

// use crate::flpcp::Proof;

// pub fn create_bgin19_proof<E: PairingEngine, R: RngCore>(
//     inputs: Vec<Vec<Vec<E::Fr>>>,
//     M: usize,
//     L: usize,
//     thetas: &Vec<E::Fr>,
//     betas: &Vec<E::Fr>,
//     lag_values: &Vec<Vec<E::Fr>>,
//     rng: &mut R,
// ) -> (Proof<E>, Proof<E>, Vec<DensePolynomial<E::Fr>>) {
    
//     // inputs[0], ..., inputs[5]
//     assert_eq!(inputs.len(), 6);
//     assert_eq!(inputs[0].len(), L);
//     assert_eq!(inputs[0][0].len(), M);
    
//     let one = E::Fr::one();
//     let zero = E::Fr::zero();

//     // Check inputs, no problem
//     // for j in 0..m {
//     //     let r = 
//     //         inputs[0][j]*inputs[2][j] + inputs[0][j]*inputs[3][j] + inputs[1][j]*inputs[2][j] + inputs[4][j] - inputs[5][j];
//     //     assert_eq!(r, zero);
//     // }

//     let domain: GeneralEvaluationDomain<E::Fr> =
//         EvaluationDomain::<E::Fr>::new(M+1).unwrap();
//     let domain_size = domain.size();
//     println!("domain size: {}", domain_size);

//     // Compute p(x), needs M more points
//     let mut p_poly = DensePolynomial::<E::Fr>::zero();
//     for i in 0..6 {

//     }

//     // let ws_share1: Vec<E::Fr> = (0..6*L).map(|_| E::Fr::rand(rng)).collect();
//     let p_coeffs_share1: Vec<E::Fr> = (0..p_poly.coeffs.len()).map(|_| E::Fr::rand(rng)).collect();

//     // let ws_share2 = (0..6*L).map(|i| ws[i] - ws_share1[i]).collect();
//     let p_coeffs_share2 = (0..p_poly.coeffs.len()).map(|i| p_poly.coeffs[i] - p_coeffs_share1[i]).collect();

//     let proof1 = Proof {
//         // ws: ws_share1,
//         p_coeffs_shares: p_coeffs_share1,
//     };
    
//     let proof2 = Proof {
//         // ws: ws_share2,
//         p_coeffs_shares: p_coeffs_share2,
//     };

//     (proof1, proof2, f_polys)
// }
