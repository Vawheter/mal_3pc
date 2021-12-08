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

use crate::fliop::{
    Proof, 
    create_bgin19_proof,
    gen_vermsg, verify_bgin19_proof,
};

// We'll use these interfaces to construct our circuit.

fn mul_local<F: Field>(xi: &F, xi_1: &F, yi: &F, yi_1: &F, alphai: &F) -> F {
    let mul1 = xi.mul(yi);
    let mul2 = xi.mul(yi_1);
    let mul3 = yi.mul(xi_1);
    let z1 = mul1.add(mul2).add(mul3).add(alphai);
    z1
}

#[test]
fn bgin19_mul_fliop() {
    use ark_serialize::*;

    let mut m: usize = 100000;
    m = m.next_power_of_two();
    let logm = ark_std::log2(m) as usize;
    let rng = &mut test_rng();

    let mut inputs: Vec<Vec<Fr>> = (0..5).map(|_| (0..m).map(|_| Fr::rand(rng) ).collect()).collect();
    let zis = (0..m).map(|j| 
        mul_local(&inputs[0][j], &inputs[1][j], &inputs[2][j], &inputs[3][j], &inputs[4][j])
    ).collect::<Vec<_>>();
    inputs.push(zis);

    let mut inputsi_1: Vec<Vec<Fr>> = (0..6).map(|_| (0..m).map(|_| Fr::rand(rng) ).collect()).collect();
    let mut inputsi_2: Vec<Vec<Fr>> = (0..6).map(|k| (0..m).map(|j| inputs[k][j] - inputsi_1[k][j] ).collect()).collect();
    
    let rs= (0..logm).map(|_| Fr::rand(rng)).collect::<Vec<_>>();

    let domain_3: GeneralEvaluationDomain<Fr> =
        EvaluationDomain::<Fr>::new(3).unwrap();
    let lag_vals_domain_3 = (0..logm).map(|i| 
        domain_3.evaluate_all_lagrange_coefficients(rs[i])
    ).collect::<Vec<_>>();
    
    let domain_5: GeneralEvaluationDomain<Fr> =
    EvaluationDomain::<Fr>::new(5).unwrap();
    let lag_vals_domain_5 = (0..logm).map(|i| 
        domain_5.evaluate_all_lagrange_coefficients(rs[i])
    ).collect::<Vec<_>>();

    let ws= (0..logm).map(|_| Fr::rand(rng)).collect::<Vec<_>>();
    let wsi_1= (0..logm).map(|_| Fr::rand(rng)).collect::<Vec<_>>();
    let wsi_2= (0..logm).map(|i| ws[i] - wsi_1[i]).collect::<Vec<_>>();
    let gammas= (0..logm).map(|_| Fr::rand(rng)).collect::<Vec<_>>();

    let prove_start = Instant::now();
    // let (proofi_1, proofi_2) = create_bgin19_proof_points::<Bls12_381, _>(inputs, M, L, &lag_values, rng);
    let (proofi_1, proofi_2) = create_bgin19_proof::<Bls12_381, _>(inputs, &rs, &ws, rng);
    let prove_time = prove_start.elapsed();
    
    let mut proofi_1_bytes = vec![];
    proofi_1.serialize(&mut proofi_1_bytes).unwrap();
    let mut proofi_2_bytes = vec![];
    proofi_2.serialize(&mut proofi_2_bytes).unwrap();
    println!("[FLIOP] Proof length: {}", proofi_1_bytes.len() + proofi_2_bytes.len());
    println!("Proving time: {:?}", prove_time);

    // // Two verifiers
    let verify_start = Instant::now();
    let pi_1_vermsg = gen_vermsg(&proofi_1, &inputsi_1, &gammas, &lag_vals_domain_3, &lag_vals_domain_5, &wsi_1, &rs);
    let result = verify_bgin19_proof(pi_1_vermsg, proofi_2, &inputsi_2, &gammas, &lag_vals_domain_3, &lag_vals_domain_5, &wsi_2, &rs);
    let verify_time = verify_start.elapsed();
    assert!(result);

    println!("Verifying time: {:?}", verify_time);
    // print!("Proof verified")
}