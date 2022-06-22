use std::vec;

use ark_ec::PairingEngine;
use ark_ff::FftField;
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain, Polynomial, UVPolynomial};
use ark_std::rand::RngCore;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::flpcp_opt::Proof;
use ark_poly_commit::kzg10::KZG10;

pub type Kzg10ComKey<'a, E> = ark_poly_commit::kzg10::Powers<'a, E>;

pub fn create_bgin19_proof<E: PairingEngine, R: RngCore>(
    inputs: Vec<Vec<E::Fr>>,
    kzg10_ck: &Kzg10ComKey<'_, E>,
    // thetas: &Vec<E::Fr>,
    // betas: &Vec<E::Fr>,
    r: E::Fr,
    rng: &mut R,
) -> Proof<E> {
    
    // inputs[0], ..., inputs[5]
    let m = inputs[0].len();

    // Check inputs, no problem
    // for j in 0..m {
    //     let r = 
    //         inputs[0][j]*inputs[2][j] + inputs[0][j]*inputs[3][j] + inputs[1][j]*inputs[2][j] + inputs[4][j] - inputs[5][j];
    //     assert_eq!(r, zero);
    // }

    let domain: GeneralEvaluationDomain<E::Fr> =
        EvaluationDomain::<E::Fr>::new(m+1).unwrap();
    let domain_size = domain.size();
    println!("domain size: {}", domain_size);

    // println!("inputs: {:?}", inputs);
    // for i in 0..6 {
    //     // Compute f(x) polynomial
    //     let f_coeffs = domain.ifft(&inputs[i]);
    //     let f_poly = DensePolynomial::from_coefficients_vec(f_coeffs);
    //     // println!("f_poly.len(): {}", f_poly.coeffs.len());
    //     f_polys.push(f_poly);
    // }
    
    // Compute q(x)
    let mut q_points = vec![];
    let multiplicative_generator = E::Fr::multiplicative_generator();
    let coset_points = inputs.iter().map(|input_i|
        input_i.iter().map(|input_ij|
            *input_ij * multiplicative_generator
        ).collect::<Vec<_>>()
    ).collect::<Vec<_>>();

    for j in 0..m {
        let qi = coset_points[0][j] * coset_points[2][j] + coset_points[0][j] * coset_points[3][j] + coset_points[1][j] * coset_points[2][j] + coset_points[4][j] - coset_points[5][j];
        q_points.push(qi)
    }
    domain.divide_by_vanishing_poly_on_coset_in_place(&mut q_points);
    let q_coeffs = domain.coset_ifft(&q_points);
    let q_poly = DensePolynomial::from_coefficients_vec(q_coeffs);

    let hiding_bound = Some(1);
    let (q_comm, q_rand) = KZG10::<E, DensePolynomial<E::Fr>>::commit(&kzg10_ck, &q_poly, hiding_bound, Some(rng)).unwrap();
    let q_r_value = q_poly.evaluate(&r);
    let q_r_proof = KZG10::<E, DensePolynomial<E::Fr>>::open(&kzg10_ck, &q_poly, r, &q_rand).unwrap();

    let proof = Proof {  
        q_comm: q_comm, 
        q_r_value: q_r_value, 
        q_r_proof: q_r_proof,
    };

    proof
}
