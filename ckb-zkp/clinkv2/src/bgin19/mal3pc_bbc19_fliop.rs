use ark_bls12_381::{Bls12_381, Fr};
use ark_ff::{Field, One, Zero};
use ark_ec::PairingEngine;
use rand::Rng;
use ark_std::test_rng;
use std::cmp::max_by;
use std::time::{Duration, Instant};
use crate::bbc19_fliop::{
    Proof, prove_dot_prod, 
    // verify_dot_prod, verify_dot_prod_helper,
};
use ark_std::UniformRand;
use ark_std::rand::RngCore;
use ark_std::rand::rngs::StdRng;
use ark_std::rand::Rng as OtherRng;

use ark_poly_commit::kzg10::KZG10;
use ark_poly::polynomial::univariate::DensePolynomial;

// We'll use these interfaces to construct our circuit.

#[test]
#[cfg(feature = "parallel")]
fn mal3pc_mul_bbc19_fliop() {
    use ark_serialize::*;

    let T: usize = 10000; // 
    let k: usize = 8; 
    let n = 25;
    let L = 4 * n + 3;

    let rng = &mut test_rng();
    let sid: usize = rng.gen();

    type KZG_Bls12_381 = KZG10<Bls12_381, DensePolynomial<Fr>>;

    let zero = Fr::zero();
    // let degree = k.next_power_of_two();
    // let setup_start = Instant::now();
    // let kzg10_pp = KZG_Bls12_381::setup(degree, false, rng).unwrap();
    // let (kzg10_ck, kzg10_vk) = KZG_Bls12_381::trim(&kzg10_pp, degree).unwrap();
    // let setup_time = setup_start.elapsed();
    // println!("Setup time: {:?}", setup_time);

    let mut inputs: Vec<Vec<Fr>> = (0..L - 1).map(|_| (0..T).map(|_| Fr::rand(rng)).collect() ).collect();
    let mut results = vec![];
    for t in 0..T {
        let mut res = zero;
        for i in 0..n {
            let cur = 4 * i;
            res += inputs[cur][t] * (inputs[cur + 2][t] + inputs[cur + 3][t]) + inputs[cur][t] * inputs[cur + 2][t];
        }
        res += inputs[L - 3][t] + inputs[L - 2][t];
        results.push(res);
    }
    inputs.push(results);

    let inputsi_1: Vec<Vec<Fr>> = (0..L).map(|_| (0..T).map(|_| Fr::rand(rng)).collect() ).collect();
    let inputsi_2: Vec<Vec<Fr>> = (0..L).map(|l| (0..T).map(|t| inputs[l][t] - inputsi_1[l][t] ).collect() ).collect();

    let rhos = vec![zero; L];

    let prove_start = Instant::now();
    let proof = prove_dot_prod::<Bls12_381, _>(&inputs, &rhos, k, sid, rng);
    let prove_time = prove_start.elapsed();
    println!("Proving time: {:?}", prove_time);

    let mut proof_bytes = vec![];
    proof.serialize(&mut proof_bytes).unwrap();
    println!("[BBC+19 FLIOP] Proof length: {}", proof_bytes.len());


    // Two verifiers
    // let verify_start = Instant::now();
    // let pi_1_vermsg = verify_dot_prod_helper(&proof, inputsi_1, &rhos, k, sid);
    // let result = verify_dot_prod(pi_1_vermsg, &proof, k, &inputsi_2, &rhos, &kzg10_vk, sid);
    // let verify_time = verify_start.elapsed();
    // assert!(result);
    // println!("Verifying time: {:?}", verify_time);


    print!("Proof verified")
}



#[test]
#[cfg(feature = "parallel")]
fn sh3pc_dotprod() {
    use ark_serialize::*;

    let rng = &mut test_rng();

    let T: usize = 3000000;
    let xs: Vec<u64> = (0..T).map(|_| rng.gen()).collect();
    let ys: Vec<u64> = (0..T).map(|_| rng.gen()).collect();

    let start = Instant::now();
    let mut res = 0u64;
    for i in 0..T {
        res += xs[i] * ys[i];
    }
    let end = start.elapsed();
    println!("res: {}", res);
    println!("time: {:?}", end);

}