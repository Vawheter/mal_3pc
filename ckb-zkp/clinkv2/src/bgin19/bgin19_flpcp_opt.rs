// use ark_bls12_381::{Bls12_381, Fr};
// use ark_ff::{Field, One, Zero};
// use ark_ec::PairingEngine;
// use rand::prelude::*;
// use rand::Rng;
// use std::time::{Duration, Instant};
// use crate::flpcp_opt::{
//     Proof, create_bgin19_proof, gen_vermsg, verify_bgin19_proof,
// };
// use ark_poly::polynomial::univariate::DensePolynomial;
// use ark_poly_commit::kzg10::KZG10;
// use ark_std::UniformRand;
// use ark_std::test_rng;

// // We'll use these interfaces to construct our circuit.

// fn mul_local<F: Field>(xi: &F, xi_1: &F, yi: &F, yi_1: &F, alphai: &F) -> F {
//     let mul1 = xi.mul(yi);
//     let mul2 = xi.mul(yi_1);
//     let mul3 = yi.mul(xi_1);
//     let z1 = mul1.add(mul2).add(mul3).add(alphai);
//     z1
// }

// // const SAMPLES: usize = 100;

// // pub fn prove_mul<E: PairingEngine, R: Rng> (
// //     inputs:Vec<Vec<E::Fr>>,
// //     M: usize,
// //     L: usize,
// //     thetas: &Vec<E::Fr>,
// //     betas: &Vec<E::Fr>,
// //     rng: &mut R,
// // ) -> (Proof<E>, Proof<E>, Vec<> {
// //     let rng = &mut thread_rng();

// //     println!("Creating proof...");
    
// // }


// // pub fn verify_mul<E: PairingEngine> (
// //     proof: Proof<E>,
// //     n: usize,
// // ) -> bool {
// //     // Verifier
// // }

// #[test]
// fn bgin19_mul_flpcp_opt() {
//     use ark_serialize::*;

//     let m: usize = 1000000;
//     // let M: usize = 1000;
//     // let L = 1000;
//     // let M: usize = 100;
//     // let L = 100;
//     let rng = &mut test_rng();
//     type KZG_Bls12_381 = KZG10<Bls12_381, DensePolynomial<Fr>>;

//     let degree =  (m + 1).next_power_of_two();
//     let kzg10_pp = KZG_Bls12_381::setup(degree, false, rng).unwrap();
//     let (kzg10_ck, kzg10_vk) = KZG_Bls12_381::trim(&kzg10_pp, degree).unwrap();

//     let xis:Vec<Fr> = (0..m).map(|_| Fr::rand(rng)).collect();
//     let xi_1s:Vec<Fr> = (0..m).map(|_| Fr::rand(rng)).collect();
//     let yis:Vec<Fr> = (0..m).map(|_| Fr::rand(rng)).collect();
//     let yi_1s:Vec<Fr> = (0..m).map(|_| Fr::rand(rng)).collect();
//     let alphais:Vec<Fr> = (0..m).map(|_| Fr::rand(rng)).collect();
//     let zis:Vec<Fr> = (0..m).map(|i| mul_local(&xis[i], &xi_1s[i], &yis[i], &yi_1s[i], &alphais[i])).collect();

//     let mut inputsi: Vec<Vec<Fr>> = vec![];
//     inputsi.push(xis.clone());
//     inputsi.push(xi_1s.clone());
//     inputsi.push(yis.clone());
//     inputsi.push(yi_1s.clone());
//     inputsi.push(alphais.clone());
//     inputsi.push(zis.clone());

//     let mut inputsi_1: Vec<Vec<Fr>> = vec![];
//     let zero = Fr::zero();
//     let zeros = vec![zero; m];
//     let alphai_1s:Vec<Fr> = (0..m).map(|_| Fr::rand(rng)).collect();
//     let input_alphai_1s:Vec<Fr> = (0..m).map(|i| -alphai_1s[i]).collect();
//     inputsi_1.push(xis);
//     inputsi_1.push(zeros.clone());
//     inputsi_1.push(yis);
//     inputsi_1.push(zeros.clone());
//     inputsi_1.push(input_alphai_1s);
//     inputsi_1.push(zis.clone());

//     let mut inputsi_2: Vec<Vec<Fr>> = vec![];
//     let alphai_2s:Vec<Fr> = (0..m).map(|i| zero - alphais[i] - alphai_1s[i]).collect();
//     let input_alphai_2s:Vec<Fr> = (0..m).map(|i| -alphai_2s[i]).collect();
//     inputsi_2.push(zeros.clone());
//     inputsi_2.push(xi_1s);
//     inputsi_2.push(zeros.clone());
//     inputsi_2.push(yi_1s);
//     inputsi_2.push(input_alphai_2s);
//     inputsi_2.push(zeros);

//     // let thetas: Vec<Fr> = (0..L).map(|_| Fr::rand(rng)).collect();
//     // let betas: Vec<Fr> = (0..M).map(|_| Fr::rand(rng)).collect();

//     let min = Fr::from((m+1).next_power_of_two() as u64);
//     let mut r: Fr = Fr::rand(rng);
//     while r <= min {
//         r = Fr::rand(rng);
//     }

//     let prove_start = Instant::now();
//     // let (proofi_1, proofi_2, _) = create_bgin19_proof::<Bls12_381, _>(inputsi, &kzg10_ck, M, L, &thetas, &betas, r, rng);
//     let proof = create_bgin19_proof::<Bls12_381, _>(inputsi, &kzg10_ck, r, rng);
//     let prove_time = prove_start.elapsed();

//     println!("Proving time: {:?}", prove_time);

//     // // Two verifiers

//     // let (_, _, fi_1_polys) = create_bgin19_proof::<Bls12_381, _>(inputsi_1, &kzg10_ck, M, L, &thetas, &betas, r, rng);
//     // let (_, _, fi_2_polys) = create_bgin19_proof::<Bls12_381, _>(inputsi_2, &kzg10_ck, M, L, &thetas, &betas, r, rng);
    

//     let verify_start = Instant::now();
//     let pi_1_vermsg = gen_vermsg(&inputsi_1, r);
//     let result = verify_bgin19_proof(&kzg10_vk, &pi_1_vermsg, &proof, &inputsi_2, r);
//     let verify_time = verify_start.elapsed();
//     assert!(result);

//     let mut vermsg_bytes = vec![];
//     pi_1_vermsg.serialize(&mut vermsg_bytes).unwrap();
//     let mut proof_bytes = vec![];
//     proof.serialize(&mut proof_bytes).unwrap();
//     println!("[FLPCP_OPT] communication: {}", vermsg_bytes.len() + proof_bytes.len());

//     println!("Verifying time: {:?}", verify_time);
//     // print!("Proof verified")
// }