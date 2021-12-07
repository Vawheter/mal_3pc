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
     //gen_vermsg, verify_bgin19_proof,
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

// const SAMPLES: usize = 100;

// pub fn prove_mul<E: PairingEngine, R: Rng> (
//     inputs:Vec<Vec<E::Fr>>,
//     M: usize,
//     L: usize,
//     thetas: &Vec<E::Fr>,
//     betas: &Vec<E::Fr>,
//     rng: &mut R,
// ) -> (Proof<E>, Proof<E>, Vec<> {
//     let rng = &mut thread_rng();

//     println!("Creating proof...");
    
// }


// pub fn verify_mul<E: PairingEngine> (
//     proof: Proof<E>,
//     n: usize,
// ) -> bool {
//     // Verifier
// }

#[test]
fn bgin19_mul_flpcp_original() {
    use ark_serialize::*;

    let m:usize = 1000;
    let M: usize = 100;
    let L: usize = 100;
    let rng = &mut test_rng();

    let mut inputs: Vec<Vec<Vec<Fr>>> = (0..5).map(|_| (0..L).map(|_| (0..M).map(|_| Fr::rand(rng)).collect() ).collect()).collect();
    let zis = (0..L).map(|l| 
        (0..M).map(|i| 
            mul_local(&inputs[0][l][i], &inputs[1][l][i], &inputs[2][l][i], &inputs[3][l][i], &inputs[4][l][i])
        ).collect::<Vec<_>>()
    ).collect::<Vec<_>>();
    inputs.push(zis);
    
    let domain: GeneralEvaluationDomain<Fr> =
        EvaluationDomain::<Fr>::new(2*M).unwrap();

    let theta= Fr::rand(rng);
    let lag_values =  (M..2*M).map(|i| 
        domain.evaluate_all_lagrange_coefficients(Fr::from(i as u64))
    ).collect::<Vec<_>>();

    let prove_start = Instant::now();
    // let (proofi_1, proofi_2) = create_bgin19_proof_points::<Bls12_381, _>(inputs, M, L, &lag_values, rng);
    let (proofi_1, proofi_2) = create_bgin19_proof_fft::<Bls12_381, _>(inputs, M, L, rng);
    let prove_time = prove_start.elapsed();
    
    let mut proofi_1_bytes = vec![];
    proofi_1.serialize(&mut proofi_1_bytes).unwrap();
    let mut proofi_2_bytes = vec![];
    proofi_2.serialize(&mut proofi_2_bytes).unwrap();
    println!("[FLPCP] Proof length: {}", proofi_1_bytes.len() + proofi_2_bytes.len());
    println!("Proving time: {:?}", prove_time);

    // // Two verifiers

    // let (_, _, fi_1_polys) = create_bgin19_proof::<Bls12_381, _>(inputsi_1, M, L, &thetas, &betas, rng);
    // let (_, _, fi_2_polys) = create_bgin19_proof::<Bls12_381, _>(inputsi_2, M, L, &thetas, &betas, rng);
    
    // let min = Fr::from((M+1).next_power_of_two() as u64);
    // let mut r: Fr = Fr::rand(rng);
    // while r <= min {
    //     r = Fr::rand(rng);
    // }
    // let verify_start = Instant::now();
    // let pi_1_vermsg = gen_vermsg(proofi_1, &fi_1_polys, &betas, r, M);
    // let result = verify_bgin19_proof(pi_1_vermsg, proofi_2, &fi_2_polys, &betas, &thetas, r, M);
    // let verify_time = verify_start.elapsed();
    // assert!(result);

    // println!("Verifying time: {:?}", verify_time);
    // print!("Proof verified")
}