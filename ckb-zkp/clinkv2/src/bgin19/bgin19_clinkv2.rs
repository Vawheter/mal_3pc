use ark_bls12_381::{Bls12_381, Fr};
use ark_ff::{Field, One, Zero};
use ark_ec::PairingEngine;
use rand::prelude::*;
use rand::Rng;
use std::time::{Duration, Instant};
use crate::kzg10::kzg10::Error;
use crate::kzg10::{
    create_random_proof, verify_proof, ProveAssignment, VerifyAssignment, KZG10, Proof, UniversalParams,
};
use crate::r1cs::{ConstraintSynthesizer, ConstraintSystem, SynthesisError};

// We'll use these interfaces to construct our circuit.

fn mul_local<F: Field>(xi: &F, xi_1: &F, yi: &F, yi_1: &F, alphai: &F) -> F {
    let mul1 = xi.mul(yi);
    let mul2 = xi.mul(yi_1);
    let mul3 = yi.mul(xi_1);
    let z1 = mul1.add(mul2).add(mul3).add(alphai);
    z1
}

struct MulLocalCircuit<F: Field> {
    xi: Option<F>,
    xi_1: Option<F>,
    yi: Option<F>,
    yi_1: Option<F>,
    alphai: Option<F>,
}

impl<'a, F: Field> ConstraintSynthesizer<F> for MulLocalCircuit<F> {
    fn generate_constraints<CS: ConstraintSystem<F>>(
        self,
        cs: &mut CS,
        index: usize,
    ) -> Result<(), SynthesisError> {

        cs.alloc_input(|| "constant one", || Ok(F::one()), index)?;

        // Alloc variable xi
        let xi_val = self.xi;
        let xi_var = cs.alloc(
            || "preimage xi_val",
            || xi_val.ok_or(SynthesisError::AssignmentMissing),
            index,
        )?;

        // Alloc variable xi_1
        let xi_1_val = self.xi;
        let xi_1_var = cs.alloc(
            || "preimage xi_1_val",
            || xi_1_val.ok_or(SynthesisError::AssignmentMissing),
            index,
        )?;

        // Alloc variable yi
        let yi_val = self.xi;
        let yi_var = cs.alloc(
            || "preimage yi_val",
            || yi_val.ok_or(SynthesisError::AssignmentMissing),
            index,
        )?;

        // Alloc variable yi_1
        let yi_1_val = self.xi;
        let yi_1_var = cs.alloc(
            || "preimage yi_1_val",
            || yi_1_val.ok_or(SynthesisError::AssignmentMissing),
            index,
        )?;

        // Alloc variable alphai
        let alphai_val = self.xi;
        let alphai_var = cs.alloc(
            || "preimage alphai_val",
            || alphai_val.ok_or(SynthesisError::AssignmentMissing),
            index,
        )?;
    
        // mul1 = xi * yi
        let mul1_val = xi_val.map(|mut e| {
            e.mul_assign(&yi_val.unwrap());
            e
        });
        let mul1_var = cs.alloc(
            || "preimage mul1_var",
            || mul1_val.ok_or(SynthesisError::AssignmentMissing),
            index,
        )?;
        if index == 0 {
            cs.enforce(
                || "mul1 = xi * yi",
                |lc| lc + xi_var,
                |lc| lc + yi_var,
                |lc| lc + mul1_var,
            );
        }

        // mul2 = xi * yi_1
        let mul2_val = xi_val.map(|mut e| {
            e.mul_assign(&yi_1_val.unwrap());
            e
        });
        let mul2_var = cs.alloc(
            || "preimage mul2_val",
            || mul2_val.ok_or(SynthesisError::AssignmentMissing),
            index,
        )?;
        if index == 0 {
            cs.enforce(
                || "mul2 = xi * yi_1",
                |lc| lc + xi_var,
                |lc| lc + yi_1_var,
                |lc| lc + mul2_var,
            );
        }

        // mul3 = xi_1 * yi
        let mul3_val = xi_1_val.map(|mut e| {
            e.mul_assign(&yi_val.unwrap());
            e
        });
        let mul3_var = cs.alloc(
            || "preimage mul2_val",
            || mul3_val.ok_or(SynthesisError::AssignmentMissing),
            index,
        )?;
        if index == 0 {
            cs.enforce(
                || "mul3 = xi_1 * yi",
                |lc| lc + xi_1_var,
                |lc| lc + yi_var,
                |lc| lc + mul3_var,
            );
        }

        // zi = mul1 + mul2 + mul3 + alphai
        let zi_val = mul1_val.map(|mut e| {
            e.add_assign(&mul2_val.unwrap());
            e.add_assign(&mul3_val.unwrap());
            e.add_assign(&alphai_val.unwrap());
            e
        });
        let zi_var = cs.alloc(
            || "preimage zi_val",
            || zi_val.ok_or(SynthesisError::AssignmentMissing),
            index,
        )?;
        if index == 0 {
            cs.enforce(
                || "zi = (mul1 + mul2 + mul3 + alphai) * one",
                |lc| lc + mul1_var + mul2_var + mul3_var + alphai_var,
                |lc| lc + (F::one(), CS::one()),
                |lc| lc + zi_var,
            );
        }  

        Ok(())
    }
}

// const SAMPLES: usize = 100;

pub fn preprocess<E: PairingEngine, R: Rng> (
    degree: usize,
    rng: &mut R,
) -> Result<UniversalParams<E>, Error> {
    // let degree: usize = n.next_power_of_two();

    println!("Creating KZG10 parameters...");
    // Create parameters for our circuit
    let kzg10_pp = KZG10::<E>::setup(degree, false, rng).unwrap();
    Ok(kzg10_pp)
}


pub fn prove_mul<E: PairingEngine> (
    kzg10_pp: UniversalParams<E>,
    xis: Vec<E::Fr>,
    xi_1s: Vec<E::Fr>,
    yis: Vec<E::Fr>,
    yi_1s: Vec<E::Fr>,
    alphais: Vec<E::Fr>,
) -> Result<Proof<E>, Error> {
    let rng = &mut thread_rng();

    println!("Generating constraints (prover)...");
    // Prover
    let mut prover_pa = ProveAssignment::<E>::default();
    // Todo: check n
    let n = xis.len();
    for i in 0..n {
        // compute the output
        let zi = mul_local(&xis[i], &xi_1s[i], &yis[i], &yi_1s[i], &alphais[i]);

        // Create an instance of our circuit (with the witness)
        let c = MulLocalCircuit {
            xi: Some(xis[i]),
            xi_1: Some(xi_1s[i]),
            yi: Some(yis[i]),
            yi_1: Some(yi_1s[i]),
            alphai: Some(alphais[i]),
        };
        c.generate_constraints(&mut prover_pa, i).unwrap();
    }

    println!("Creating proof...");
    let degree: usize = n.next_power_of_two();
    let (ck, _) = KZG10::<E>::trim(&kzg10_pp, degree).unwrap();
    // Create a clinkv2 proof with our parameters.
    let proof = create_random_proof(&prover_pa, &ck, rng).unwrap();
    Ok(proof)
}


pub fn verify_mul<E: PairingEngine> (
    kzg10_pp: UniversalParams<E>,
    proof: Proof<E>,
    n: usize,
) -> Result<bool, Error> {
    // Verifier
    println!("Generating constraints (verifier)...");

    let mut verifier_pa = VerifyAssignment::<E>::default();

    // Create an instance of our circuit (with the witness)
    let verify_c = MulLocalCircuit {
        xi: None,
        xi_1: None,
        yi: None,
        yi_1: None,
        alphai: None,
    };
    verify_c
        .generate_constraints(&mut verifier_pa, 0usize)
        .unwrap();

    println!("Verifying proof...");

    // Todo: check n
    let mut io: Vec<Vec<E::Fr>> = vec![];
    let one = vec![E::Fr::one(); n];
    io.push(one);
    // Check the proof
    let degree: usize = n.next_power_of_two();
    let (_, vk) = KZG10::<E>::trim(&kzg10_pp, degree).unwrap();
    let result = verify_proof(&verifier_pa, &vk, &proof, &io).unwrap();
    Ok(result)
}

#[test]
fn bgin19_mul_clinkv2() {
    use ark_serialize::*;

    let n:usize = 1000000;
    let rng = &mut thread_rng();
    let degree: usize = n.next_power_of_two();
    let pp = preprocess::<Bls12_381, _>(degree, rng).unwrap();
    let pp2 = pp.clone();
    // let (kzg10_ck, kzg10_vk) = KZG10::<Bls12_381>::trim(&pp, degree).unwrap();

    let xis:Vec<Fr> = (0..n).map(|_| rng.gen()).collect();
    let xi_1s:Vec<Fr> = (0..n).map(|_| rng.gen()).collect();
    let yis:Vec<Fr> = (0..n).map(|_| rng.gen()).collect();
    let yi_1s:Vec<Fr> = (0..n).map(|_| rng.gen()).collect();
    let alphais:Vec<Fr> = (0..n).map(|_| rng.gen()).collect();

    let prove_start = Instant::now();
    let proof = prove_mul(pp, xis, xi_1s, yis, yi_1s, alphais).unwrap();
    let prove_time = prove_start.elapsed();
    
    let mut proof_bytes = vec![];
    proof.serialize(&mut proof_bytes).unwrap();
    println!("[Clinkv2 Kzg10] Proof length: {}", proof_bytes.len());

    let verify_start = Instant::now();
    let result = verify_mul(pp2, proof, n).unwrap();
    let verify_time = verify_start.elapsed();
    assert!(result);

    println!("Proving time: {:?}", prove_time);
    println!("Verifying time: {:?}", verify_time);
    // print!("Proof verified")
}