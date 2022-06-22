use ark_bls12_381::{Bls12_381, Fr};
use ark_ff::{Field, One, Zero};
use ark_ec::PairingEngine;
use rand::Rng;
use ark_std::test_rng;
use std::time::{Duration, Instant};
use crate::flpcp::{
    Proof, 
    create_bgin19_proof_fft,
    create_bgin19_proof_points,
    gen_vermsg, verify_bgin19_proof,
};
use ark_std::UniformRand;
use ark_std::rand::RngCore;
use ark_std::rand::rngs::StdRng;
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain, Polynomial, UVPolynomial};

// We'll use these interfaces to construct our circuit.

fn mul_local<F: Field>(xi: &F, xi_1: &F, yi: &F, yi_1: &F, alphai: &F) -> F {
    let mul1 = xi.mul(yi);
    let mul2 = xi.mul(yi_1);
    let mul3 = yi.mul(xi_1);
    let z1 = mul1.add(mul2).add(mul3).add(alphai);
    z1
}

#[test]
fn bgin19_mul_flpcp_original() {
    use ark_serialize::*;

    // 2^20 = 1048576
    // 2^23 = 8388608
    let _m:usize = 1000000;
    let M: usize = 1000;
    let L: usize = 1000;
    let rng = &mut test_rng();

    let mut inputs: Vec<Vec<Vec<Fr>>> = (0..5).map(|_| (0..L).map(|_| (0..M).map(|_| Fr::rand(rng)).collect() ).collect()).collect();
    let zis = (0..L).map(|l| 
        (0..M).map(|i| 
            mul_local(&inputs[0][l][i], &inputs[1][l][i], &inputs[2][l][i], &inputs[3][l][i], &inputs[4][l][i])
        ).collect::<Vec<_>>()
    ).collect::<Vec<_>>();
    inputs.push(zis);

    let inputsi_1: Vec<Vec<Vec<Fr>>> = (0..6).map(|_| (0..L).map(|_| (0..M).map(|_| Fr::rand(rng)).collect() ).collect()).collect();
    let inputsi_2: Vec<Vec<Vec<Fr>>> = (0..6).map(|k| (0..L).map(|l| (0..M).map(|j| inputs[k][l][j] - inputsi_1[k][l][j] ).collect() ).collect()).collect();
    
    let domain: GeneralEvaluationDomain<Fr> =
        EvaluationDomain::<Fr>::new(M).unwrap();

    let theta: Fr = Fr::rand(rng);
    let beta: Fr = Fr::rand(rng);
    let r: Fr = domain.sample_element_outside_domain(rng);
    let ws: Vec<Vec<Fr>> = (0..6).map(|_| (0..L).map(|_| Fr::rand(rng) ).collect() ).collect();
    let wsi_1: Vec<Vec<Fr>> = (0..6).map(|_| (0..L).map(|_| Fr::rand(rng) ).collect() ).collect();
    let wsi_2: Vec<Vec<Fr>> = (0..6).map(|k| (0..L).map(|l| ws[k][l] - wsi_1[k][l] ).collect() ).collect();

    // let lag_values =  (M..2*M).map(|i| 
    //     domain.evaluate_all_lagrange_coefficients(Fr::from(i as u64))
    // ).collect::<Vec<_>>();

    let prove_start = Instant::now();
    // let (proofi_1, proofi_2) = create_bgin19_proof_points::<Bls12_381, _>(inputs, M, L, &lag_values, rng);
    let (proofi_1, proofi_2) = create_bgin19_proof_fft::<Bls12_381, _>(inputs, theta, &ws, rng);
    let prove_time = prove_start.elapsed();
    println!("Proving time: {:?}", prove_time);

    let mut proofi_1_bytes = vec![];
    proofi_1.serialize(&mut proofi_1_bytes).unwrap();
    let mut proofi_2_bytes = vec![];
    proofi_2.serialize(&mut proofi_2_bytes).unwrap();

    // Two verifiers
    let verify_start = Instant::now();
    let pi_1_vermsg = gen_vermsg(proofi_1, &inputsi_1, beta, &wsi_1, r);
    let result = verify_bgin19_proof(&pi_1_vermsg, proofi_2, &inputsi_2, beta, theta, &wsi_2, r);
    let verify_time = verify_start.elapsed();
    assert!(result);
    println!("Verifying time: {:?}", verify_time);

    let mut pi_1_vermsg_bytes = vec![];
    pi_1_vermsg.serialize(&mut pi_1_vermsg_bytes).unwrap();

    println!("[FLPCP] Proof length: {}", proofi_1_bytes.len() + proofi_2_bytes.len() + pi_1_vermsg_bytes.len());

    print!("Proof verified")
}