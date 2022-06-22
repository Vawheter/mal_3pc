use ark_bls12_381::{Bls12_381, Fr};
use ark_ff::{Field, One, Zero};
use rand::Rng;
use rand::prelude::*;
use std::time::{Duration, Instant};
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain, Polynomial, UVPolynomial};
use ark_poly_commit::kzg10::KZG10;
use ark_std::test_rng;
use ark_std::UniformRand;

#[test]
fn test_fft_time() {

    let m: usize = 1<<20;
    let cnt = 10;
    let rng = &mut test_rng();

    let domain: GeneralEvaluationDomain<Fr> =
        EvaluationDomain::<Fr>::new(m).unwrap();

    let mut fft_time = Duration::new(0, 0);
    for _ in 0..cnt {
        let points:Vec<Fr> = (0..m).map(|_| Fr::rand(rng)).collect();

        let fft_start = Instant::now();
        let _ = domain.ifft(&points);
        fft_time += fft_start.elapsed();
    }
    println!("size: {:?}", m);
    println!("fft_time: {:?}", fft_time/cnt);
}


// #[test]
// fn test_commit_time() {

//     let m: usize = 1<<20;
//     let cnt = 20;
//     let rng = &mut test_rng();

//     let degree =  m.next_power_of_two();
//     let pp = KZG10::<Bls12_381, DensePolynomial<Fr>>::setup(degree, false, rng).unwrap();
//     let (kzg10_ck, _) = KZG10::<Bls12_381, DensePolynomial<Fr>>::trim(&pp, degree).unwrap();

//     let hiding_bound = Some(1);

//     let mut kzg10_com_time = Duration::new(0, 0);
//     for _ in 0..cnt {
//         let r_poly = DensePolynomial::<Fr>::rand(m, rng);
        
//         let kzg10_com_start = Instant::now();
//         let (_, _) = KZG10::<Bls12_381, DensePolynomial<Fr>>::commit(&kzg10_ck, &r_poly, hiding_bound, Some(rng)).unwrap();
//         kzg10_com_time += kzg10_com_start.elapsed();
//     }
//     println!("size: {:?}", m);
//     println!("Committing  time: {:?}", kzg10_com_time/cnt);
// }