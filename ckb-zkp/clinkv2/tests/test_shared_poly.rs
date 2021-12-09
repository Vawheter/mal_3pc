
use ark_bls12_381::{Bls12_381, Fr};
use ark_ff::{Field, One, Zero};
use ark_ec::PairingEngine;
use rand::Rng;
use ark_std::test_rng;
use std::time::{Duration, Instant};
use ark_std::UniformRand;
use ark_std::rand::RngCore;
use ark_std::rand::rngs::StdRng;
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain, Polynomial, UVPolynomial};

#[test]
fn test_shared_poly_1() {

    let mut m: usize = 8;
    m = m.next_power_of_two();
    let rng = &mut test_rng();

    let inputs: Vec<Fr> = (0..m).map(|_| Fr::rand(rng) ).collect();

    let inputs_share1: Vec<Fr> = (0..m).map(|_| Fr::rand(rng) ).collect();
    let inputs_share2: Vec<Fr> = (0..m).map(|j| inputs[j] - inputs_share1[j] ).collect();

    let domain_m: GeneralEvaluationDomain<Fr> =
    EvaluationDomain::<Fr>::new(m).unwrap();

    let q_coeffs = domain_m.ifft(&inputs);

    // println!("q_coeffs: {:?}", &q_coeffs);
    let q_coeffs_share1 = domain_m.ifft(&inputs_share1);
    let q_coeffs_share2 = domain_m.ifft(&inputs_share2);
    let q_coeffs_recons = (0..m).map(|j| q_coeffs_share1[j] + q_coeffs_share2[j]).collect::<Vec<Fr>>();
    // println!("q_coeffs_recons: {:?}", &q_coeffs_recons);
    assert_eq!(q_coeffs, q_coeffs_recons);
}


#[test]
fn test_shared_poly_2() {

    let rng = &mut test_rng();

    let inputs: Vec<Fr> = (0..1).map(|_| Fr::rand(rng) ).collect();
    let domain_3: GeneralEvaluationDomain<Fr> =
        EvaluationDomain::<Fr>::new(3).unwrap();
    let q_coeffs = domain_3.ifft(&inputs); // IMPORTANT: get 4 coeffs
    println!("q_coeffs: {:?}", &q_coeffs);
    println!("q_coeffs.len(): {:?}", &q_coeffs.len());
    let q_poly = DensePolynomial::from_coefficients_vec(q_coeffs.clone());

    let g_zero = Fr::one();
    let g_one = domain_3.group_gen();
    let g_two = g_one * g_one;

    let mut q_values = vec![];
    q_values.push(q_poly.evaluate(&g_zero));
    q_values.push(q_poly.evaluate(&g_one));
    q_values.push(q_poly.evaluate(&g_two));

    let q_coeffs_recons_full = domain_3.ifft(&q_values);
    println!("q_coeffs_recons_full.len(): {:?}", &q_coeffs_recons_full.len());

    // println!("q_coeffs_recons_full: {:?}", &q_coeffs_recons_full);

    let q_values_share1: Vec<Fr> = (0..3).map(|_| Fr::rand(rng) ).collect();
    let q_values_share2: Vec<Fr> = (0..3).map(|j| q_values[j] - q_values_share1[j] ).collect();

    let q_coeffs_share1 = domain_3.ifft(&q_values_share1);
    let q_coeffs_share2 = domain_3.ifft(&q_values_share2);
    
    let q_coeffs_recons_share = (0..q_coeffs_share1.len()).map(|j| q_coeffs_share1[j] + q_coeffs_share2[j]).collect::<Vec<Fr>>();
    // println!("q_coeffs_recons_share: {:?}", &q_coeffs_recons_share);
    println!("q_coeffs_recons_share.len(): {:?}", &q_coeffs_recons_share.len());

    assert_eq!(q_coeffs_recons_full, q_coeffs_recons_share);
}

#[test]
fn test_group_gen() {

    let domain_2: GeneralEvaluationDomain<Fr> =
        EvaluationDomain::<Fr>::new(2).unwrap();
    println!("domain_2.group_gen(): {:?}", domain_2.group_gen());
    

    let domain_3: GeneralEvaluationDomain<Fr> =
        EvaluationDomain::<Fr>::new(3).unwrap();
    println!("domain_3.group_gen(): {:?}", domain_3.group_gen());

    let domain_5: GeneralEvaluationDomain<Fr> =
        EvaluationDomain::<Fr>::new(5).unwrap();
    println!("domain_5.group_gen(): {:?}", domain_5.group_gen());

}