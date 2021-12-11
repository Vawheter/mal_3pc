use ark_bls12_381::{Bls12_381, Fr};
use ark_ff::{Field, One, Zero};
use ark_ec::PairingEngine;
use rand::Rng;
use ark_std::test_rng;
use std::time::{Duration, Instant};
use crate::plonk_sumcheck::{
    Proof, create_bgin19_proof, 
    gen_vermsg, verify_bgin19_proof,
};
use ark_std::UniformRand;
use ark_std::rand::RngCore;
use ark_std::rand::rngs::StdRng;

use ark_poly_commit::kzg10::KZG10;
use ark_poly::polynomial::univariate::DensePolynomial;

// We'll use these interfaces to construct our circuit.

fn mul_local<F: Field>(xi: &F, xi_1: &F, yi: &F, yi_1: &F, alphai: &F) -> F {
    let mul1 = xi.mul(yi);
    let mul2 = xi.mul(yi_1);
    let mul3 = yi.mul(xi_1);
    let z1 = mul1.add(mul2).add(mul3).add(alphai);
    z1
}

#[test]
fn bgin19_mul_plonk() {
    use ark_serialize::*;

    let m:usize = 100;
    let M: usize = 10;
    let L: usize = 10;
    let rng = &mut test_rng();

    type KZG_Bls12_381 = KZG10<Bls12_381, DensePolynomial<Fr>>;

    let degree =  (2*L).next_power_of_two();
    let kzg10_pp = KZG_Bls12_381::setup(degree, false, rng).unwrap();
    let (kzg10_ck, kzg10_vk) = KZG_Bls12_381::trim(&kzg10_pp, degree).unwrap();

    let mut inputs: Vec<Vec<Vec<Fr>>> = (0..5).map(|_| (0..L).map(|_| (0..M).map(|_| Fr::rand(rng)).collect() ).collect()).collect();
    let zis = (0..L).map(|l| 
        (0..M).map(|i| 
            mul_local(&inputs[0][l][i], &inputs[1][l][i], &inputs[2][l][i], &inputs[3][l][i], &inputs[4][l][i])
        ).collect::<Vec<_>>()
    ).collect::<Vec<_>>();
    inputs.push(zis);

    let inputsi_1: Vec<Vec<Vec<Fr>>> = (0..6).map(|_| (0..L).map(|_| (0..M).map(|_| Fr::rand(rng)).collect() ).collect()).collect();
    let inputsi_2: Vec<Vec<Vec<Fr>>> = (0..6).map(|k| (0..L).map(|l| (0..M).map(|j| inputs[k][l][j] - inputsi_1[k][l][j] ).collect() ).collect()).collect();

    let r = Fr::rand(rng);
    let z = Fr::rand(rng);
    // let eta = Fr::rand(rng);
    let eta = Fr::one();

    let prove_start = Instant::now();
    let proof = create_bgin19_proof::<Bls12_381, _>(inputs, &kzg10_ck, eta, r, z, rng);
    let prove_time = prove_start.elapsed();
    println!("Proving time: {:?}", prove_time);

    let mut proof_bytes = vec![];
    proof.serialize(&mut proof_bytes).unwrap();

    // Two verifiers
    let verify_start = Instant::now();
    let pi_1_vermsg = gen_vermsg(&inputsi_1, eta, r, z);
    let result = verify_bgin19_proof(pi_1_vermsg, proof, &inputsi_2, kzg10_vk, eta, r, z);
    let verify_time = verify_start.elapsed();
    assert!(result);
    println!("Verifying time: {:?}", verify_time);

    println!("[PLONK] Proof length: {}", proof_bytes.len());

    print!("Proof verified")
}