use ark_bls12_381::{Bls12_381, Fr};
use ark_ff::{Field, One, Zero};
use rand::Rng;
use rand::prelude::*;
use std::time::{Duration, Instant};
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain, Polynomial, UVPolynomial};

use crate::{
    kzg10::{Kzg10Proof, ProveAssignment, ProveKey, KZG10},
};

#[test]
fn test_fft() {

    let _m:usize = 100000000;
    let M = 10000;
    let L = 10000;
    let rng = &mut thread_rng();

    let domain: GeneralEvaluationDomain<Fr> =
        EvaluationDomain::<Fr>::new(M+1).unwrap();

    let mut fft_time = Duration::new(0, 0);
    for i in 0..L {
        let points:Vec<Fr> = (0..M).map(|_| rng.gen()).collect();

        let fft_start = Instant::now();
        let _ = domain.ifft(&points);
        fft_time += fft_start.elapsed();
    }

    println!("Verifying time: {:?}", fft_time);
    // print!("Proof verified")
}


#[test]
fn test_kzg10() {

    let m:usize = 1000000;
    let L = 1;
    let rng = &mut thread_rng();

    let degree =  (m + 1).next_power_of_two();
    let pp = KZG10::<Bls12_381>::setup(degree, false, rng).unwrap();
    let (kzg10_ck, _) = KZG10::<Bls12_381>::trim(&pp, degree).unwrap();

    let domain: GeneralEvaluationDomain<Fr> =
        EvaluationDomain::<Fr>::new(m+1).unwrap();

    let mut kzg10_commit_time = Duration::new(0, 0);
    let hiding_bound = Some(1);
    for i in 0..L {
        let r_poly = DensePolynomial::<Fr>::rand(m, rng);
        
        let kzg10_commit_start = Instant::now();
        let (_, _) = KZG10::<Bls12_381>::commit(&kzg10_ck, &r_poly, hiding_bound, Some(rng)).unwrap();
        kzg10_commit_time += kzg10_commit_start.elapsed();
    }

    println!("Verifying time: {:?}", kzg10_commit_time);
    // print!("Proof verified")
}