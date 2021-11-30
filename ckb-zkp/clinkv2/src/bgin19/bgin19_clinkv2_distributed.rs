use ark_bls12_381::{Bls12_381, Fr};
use ark_ff::{Field, One, Zero};
use ark_ec::PairingEngine;
use rand::prelude::*;
use rand::Rng;
use std::time::{Duration, Instant};
use crate::kzg10_distributed::{
    Proof, create_random_proof, verify_proof, gen_witness_polys, ProveAssignment, VerifyAssignment, UniversalParams, KZG10, kzg10::Error,
};
use crate::r1cs::{ConstraintSynthesizer, ConstraintSystem, SynthesisError};
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


pub fn prove_mul<E: PairingEngine, R: Rng> (
    kzg10_pp: UniversalParams<E>,
    inputs: Vec<Vec<E::Fr>>,
    M: usize,
    L: usize,
    rhos: Vec<E::Fr>,
    eta: E::Fr,
    r: E::Fr,
    degree: usize,
    rng: &mut R
) -> Result<Proof<E>, Error> {
    // let rng = &mut thread_rng();
    let m = M * L;
    println!("Generating constraints (prover)...");
    // Prover
    let mut prover_pa = ProveAssignment::<E>::default();
    for i in 0..m {
        // compute the output
        let zi = mul_local(&inputs[0][i], &inputs[1][i], &inputs[2][i], &inputs[3][i], &inputs[4][i]);

        // Create an instance of our circuit (with the witness)
        let c = MulLocalCircuit {
            xi: Some(inputs[0][i]),
            xi_1: Some(inputs[1][i]),
            yi: Some(inputs[2][i]),
            yi_1: Some(inputs[3][i]),
            alphai: Some(inputs[4][i]),
        };
        c.generate_constraints(&mut prover_pa, i).unwrap();
    }

    println!("Creating proof...");
    let (ck, _) = KZG10::<E>::trim(&kzg10_pp, degree).unwrap();
    // Create a clinkv2 proof with our parameters.
    let proof = create_random_proof(&prover_pa, &ck, M, L, rhos, eta, r, rng).unwrap();
    Ok(proof)
}

pub fn verify_mul<E: PairingEngine> (
    kzg10_pp: UniversalParams<E>,
    proof: Proof<E>,
    witnesses: Vec<Vec<E::Fr>>,
    vanishing_r_value: E::Fr,
    eta: E::Fr,
    r: E::Fr,
    M: usize,
    L: usize,
) -> Result<bool, Error> {
    // Verifier
    println!("Generating constraints (verifier)...");

    let mut verifier_pa = VerifyAssignment::<E>::default();

    let w_polys = gen_witness_polys::<E>(&witnesses, M, L);

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

    // Check the proof
    let degree: usize = (2 * (M + 1)).next_power_of_two();
    let (_, vk) = KZG10::<E>::trim(&kzg10_pp, degree).unwrap();
    let result = verify_proof(&verifier_pa, &vk, &proof, r, &w_polys, vanishing_r_value, eta).unwrap();
    Ok(result)
}

#[test]
fn bgin19_mul_flpcp_opt() {
    use ark_serialize::*;

    let m:usize = 1000000;
    let M: usize = 100;
    let L = 10000;
    // let M: usize = 100;
    // let L = 100;
    let T = 10;
    let rng = &mut thread_rng();

    let degree =  2 * ((M + 1).next_power_of_two());
    let pp = preprocess::<Bls12_381, _>(degree, rng).unwrap();
    let (kzg10_ck, kzg10_vk) = KZG10::<Bls12_381>::trim(&pp, degree).unwrap();

    let xis: Vec<Fr> = (0..m).map(|_| rng.gen()).collect();
    let xi_1s: Vec<Fr> = (0..m).map(|_| rng.gen()).collect();
    let yis: Vec<Fr> = (0..m).map(|_| rng.gen()).collect();
    let yi_1s: Vec<Fr> = (0..m).map(|_| rng.gen()).collect();
    let alphais: Vec<Fr> = (0..m).map(|_| rng.gen()).collect();
    let zis: Vec<Fr> = vec![];
    let zis:Vec<Fr> = (0..m).map(|i| mul_local(&xis[i], &xi_1s[i], &yis[i], &yi_1s[i], &alphais[i])).collect();

    let mut inputsi: Vec<Vec<Fr>> = vec![];
    inputsi.push(xis.clone());
    inputsi.push(xi_1s.clone());
    inputsi.push(yis.clone());
    inputsi.push(yi_1s.clone());
    inputsi.push(alphais.clone());
    inputsi.push(zis.clone());

    let rhos_i: Vec<Fr> = (0..L*T).map(|_| rng.gen()).collect();
    let rhos_i_1: Vec<Fr> = (0..L*T).map(|_| rng.gen()).collect();
    let rhos_i_2: Vec<Fr> = (0..L*T).map(|i| rhos_i[i] - rhos_i_1[i]).collect();

    let eta: Fr = rng.gen();
    let r: Fr = rng.gen();

    let domain: GeneralEvaluationDomain<Fr> =
        EvaluationDomain::<Fr>::new(M).ok_or(SynthesisError::PolynomialDegreeTooLarge).unwrap(); 
    let lag_values = domain.evaluate_all_lagrange_coefficients(r);

    let vanishing_poly = domain.vanishing_polynomial();
    let vanishing_r_value = vanishing_poly.evaluate(&r);


    let prove_start = Instant::now();
    let proof = prove_mul(pp.clone(), inputsi, M, L, rhos_i_1, eta, r, degree, rng).unwrap();
    let prove_time = prove_start.elapsed();
    
    let mut proof_bytes = vec![];
    proof.serialize(&mut proof_bytes).unwrap();
    println!("[Clinkv2 distributed] Proof length: {}", proof_bytes.len());

    // // Two verifiers

    let mut inputs_i_1: Vec<Vec<Fr>> = vec![];
    let zero = Fr::zero();
    let zeros = vec![zero; m];
    let alphai_1s:Vec<Fr> = (0..m).map(|_| rng.gen()).collect();
    let input_alphai_1s:Vec<Fr> = (0..m).map(|i| -alphai_1s[i]).collect();
    inputs_i_1.push(xis);
    inputs_i_1.push(zeros.clone());
    inputs_i_1.push(yis);
    inputs_i_1.push(zeros.clone());
    inputs_i_1.push(input_alphai_1s);
    inputs_i_1.push(zis.clone());

    let mut inputs_i_2: Vec<Vec<Fr>> = vec![];
    let alphai_2s:Vec<Fr> = (0..m).map(|i| zero - alphais[i] - alphai_1s[i]).collect();
    let input_alphai_2s:Vec<Fr> = (0..m).map(|i| -alphai_2s[i]).collect();
    inputs_i_2.push(zeros.clone());
    inputs_i_2.push(xi_1s);
    inputs_i_2.push(zeros.clone());
    inputs_i_2.push(yi_1s);
    inputs_i_2.push(input_alphai_2s);
    inputs_i_2.push(zeros);

    let w_polys_i_1 = gen_witness_polys::<Bls12_381>(&inputs_i_1, M, L); 
    // To P_{i+1}, has no x_{i-1} and y_{i-1}, index 1 and 3
    // let mut w_sec_polys_i_2 = vec![];
    // w_sec_polys_i_2.push(w_sec_polys_i_2[1].clone());
    // w_sec_polys_i_2.push(w_sec_polys_i_2[3].clone());
    // let w_sec_r_values = 
    let w_polys_i_2 = gen_witness_polys::<Bls12_381>(&inputs_i_2, M, L); 
    
    let verify_start = Instant::now();
    let result = verify_mul(pp, proof, inputs_i_1, vanishing_r_value, eta, r, M, L).unwrap();
    let verify_time = verify_start.elapsed();
    assert!(result);

    println!("Verifying time: {:?}", verify_time);
    // print!("Proof verified")
}