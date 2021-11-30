#![allow(non_snake_case)]

use ark_ec::PairingEngine;
use ark_ff::{One, Zero};
use ark_poly::polynomial::univariate::{DensePolynomial, SparsePolynomial};
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain, Polynomial, UVPolynomial};
use rand::Rng;

// DEV
//use std::time::{Duration, Instant};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::{
    kzg10_distributed::{Proof, ProveAssignment, ProveKey, KZG10},
    r1cs::{Index, SynthesisError},
};

pub fn gen_witness_polys<E: PairingEngine> (
    witnesses: &Vec<Vec<E::Fr>>,
    M: usize,
    L: usize,
) -> Vec<Vec<DensePolynomial<E::Fr>>> {
    let T = witnesses.len() - 1;

    let domain: GeneralEvaluationDomain<E::Fr> =
        EvaluationDomain::<E::Fr>::new(M).ok_or(SynthesisError::PolynomialDegreeTooLarge).unwrap(); 

    let mut w_polys:Vec<Vec<DensePolynomial<E::Fr>>> = vec![];
    let one_coeffs = domain.ifft(&witnesses[0][0..M]);
    let one_poly = DensePolynomial::from_coefficients_vec(one_coeffs);
    let mut w_0_polys = vec![];
    for _ in 0..L {
        w_0_polys.push(one_poly.clone());
    }
    w_polys.push(w_0_polys);

    for t in 0..T {
        //let start = Instant::now();
        let mut w_t_polys = vec![];
        for l in 0..L {
            let w_lt_coeffs = domain.ifft(&witnesses[t][l*M..(l+1)*M]);
            let w_lt_poly = DensePolynomial::from_coefficients_vec(w_lt_coeffs);
            w_t_polys.push(w_lt_poly);
        }
        w_polys.push(w_t_polys);
    }

    w_polys
}



pub fn create_random_proof<E: PairingEngine, R: Rng>(
    circuit: &ProveAssignment<E>,
    kzg10_ck: &ProveKey<'_, E>,
    M: usize,
    L: usize,
    rhos: Vec<E::Fr>,//length: TL
    eta: E::Fr,
    r: E::Fr,
    rng: &mut R,
) -> Result<Proof<E>, SynthesisError> {

    let t_io = circuit.input_assignment.len();
    assert_eq!(t_io, 1);

    // Number of all variables
    let T = circuit.aux_assignment.len() + t_io;
    // Number of constraints
    let K = circuit.at.len();
    // Number of copies
    let m = circuit.input_assignment[0].len();
    assert!(m <= M * L);

    // Compute and witness polynomials, degree: M - 1
    let domain: GeneralEvaluationDomain<E::Fr> =
        EvaluationDomain::<E::Fr>::new(M).ok_or(SynthesisError::PolynomialDegreeTooLarge)?;

    let domain_size = domain.size();
    // println!("domain_size: {:?}", domain_size);

    let witnesses = [&circuit.input_assignment[..], &circuit.aux_assignment[..]].concat();
    // Every row contains w_{1t}(x) to w_{lt}(x), in total T rows, L columns
    let mut w_polys = gen_witness_polys::<E>(&witnesses, M, L);

    let zero = E::Fr::zero();
    let one = E::Fr::one();
    let hiding_bound = Some(2);

    //let mut rj_commit_time = Duration::new(0, 0);
    //let mut rj_ifft_time = Duration::new(0, 0);

    // let vanishing_poly = domain.vanishing_polynomial();

    let mut cnt = 0;
    for t in 0..T {
        //let start = Instant::now();
        for l in 0..L {
            let rho_vanishing_coeffs = vec![(0, -rhos[cnt]), (domain_size, rhos[cnt])];
            let rho_vanishing_poly = SparsePolynomial::from_coefficients_vec(rho_vanishing_coeffs);
            w_polys[t][l] += &rho_vanishing_poly.into();
            cnt += 1;
        }
    }

    // for t in 0..T_io {
    //     //let start = Instant::now();
    //     let mut w_t_polys = vec![];
    //     for l in 0..L {
    //         let w_io_lt_coeffs = domain.ifft(&circuit.input_assignment[t][l*M..(l+1)*M]);
    //         let mut w_io_lt_poly = DensePolynomial::from_coefficients_vec(w_io_lt_coeffs);
    //         let rho = E::Fr::rand(rng);
    //         let rho_vanishing_coeffs = vec![(0, -rho), (domain_size, rho)];
    //         let rho_vanishing_poly = SparsePolynomial::from_coefficients_vec(rho_vanishing_coeffs);
    //         w_io_lt_poly += &rho_vanishing_poly;
    //         w_t_polys.push(w_io_lt_poly);
    //     }
    //     w_polys.push(w_t_polys);
    // }

    // for t in 0..T_mid {
    //     //let start = Instant::now();
    //     let mut w_t_polys = vec![];
    //     for l in 0..L {
    //         let w_mid_lt_coeffs = domain.ifft(&circuit.aux_assignment[t][l*M..(l+1)*M]);
    //         let w_mid_lt_poly = DensePolynomial::from_coefficients_vec(w_mid_lt_coeffs);
    //         let rho = E::Fr::rand(rng);
    //         let rho_vanishing_coeffs = vec![(0, -rho), (domain_size, rho)];
    //         let rho_vanishing_poly = SparsePolynomial::from_coefficients_vec(rho_vanishing_coeffs);
    //         w_mid_lt_poly += &rho_vanishing_poly;
    //         w_t_polys.push(w_mid_lt_poly);
    //     }
    //     w_polys.push(w_t_polys);
    // }
    
    //println!("rj_ifft_time: {:?}", rj_ifft_time);
    //println!("rj_commit_time: {:?}", rj_commit_time);

    // Compute and commit quotient polynomials
    // let m_abc = circuit.at.len();

    let mut eta_lk = one;

    //let mut q_commit_time = Duration::new(0, 0);
    //let mut abci_fft_time = Duration::new(0, 0);
    //let start = Instant::now();

    // println!("r_q_polys: {:?} * {:?}", &r_q_polys.len(), &r_q_polys[0].len());
    // println!("r_q_polys: {:?}", &r_q_polys);

    let mut batch_coset_ab = vec![zero; domain_size];
    let mut batch_c_values = vec![zero; domain_size];

    for l in 0.. L {
        for k in 0..K {
            let mut a_lk_coeffs = vec![zero; domain_size];
            for (coeff, index) in (&circuit.at[k]).into_iter() {
                let t = match index {
                    Index::Input(j) => *j,
                    Index::Aux(j) => *j,
                };
                for i in 0..domain_size {
                    a_lk_coeffs[i] += &(w_polys[t][l].coeffs[i] * coeff);
                }
            }
            let mut a_lk_poly = DensePolynomial::from_coefficients_vec(a_lk_coeffs);
    
            let mut b_lk_coeffs = vec![zero; domain_size];
            for (coeff, index) in (&circuit.bt[k]).into_iter() {
                let t = match index {
                    Index::Input(j) => *j,
                    Index::Aux(j) => *j,
                };
                for i in 0..domain_size {
                    b_lk_coeffs[i] += &(w_polys[t][l].coeffs[i] * coeff);
                }
            }
            let mut b_lk_poly = DensePolynomial::from_coefficients_vec(b_lk_coeffs);
    
    
            domain.coset_fft_in_place(&mut a_lk_poly.coeffs);
            domain.coset_fft_in_place(&mut b_lk_poly.coeffs);
    
            // on coset: n values of a*b on coset
            let coset_ab_values = domain.mul_polynomials_in_evaluation_domain(&a_lk_poly, &b_lk_poly);
    
            drop(a_lk_poly);
            drop(b_lk_poly);
    
            // on coset: n values of \sum{eta^i * ab} on coset
            for j in 0..domain_size {
                batch_coset_ab[j] += coset_ab_values[j] * eta_lk;
            }
            // cfg_iter!(coset_ab_values)
            //     .zip(&mut batch_coset_ab)
            //     .for_each(|(coset_ablk_j, batch_coset_ab_j)| *batch_coset_ab_j += &(eta_lk * (*coset_ablk_j)));
    
            let mut c_lk_values = vec![zero; domain_size];
            for (coeff, index) in (&circuit.ct[k]).into_iter() {
                let t = match &index {
                    Index::Input(j) => *j,
                    Index::Aux(j) => *j,
                };
                for j in 0..domain_size {
                    c_lk_values[j] += circuit.aux_assignment[t][j] * coeff;
                }
            }
            
            // on original domain: n values of \sum{eta^i * c} on original domain
            for j in 0..domain_size {
                batch_c_values[j] += c_lk_values[j] * eta_lk;
            }
            // cfg_iter!(c_lk_values)
            //     .zip(&mut batch_c_values)
            //     .for_each(|(c_lk_j, batch_c_j)| *batch_c_j += &(eta_lk * (*c_lk_j)));
    
            eta_lk = eta_lk * &eta;
        }
    }

    domain.ifft_in_place(&mut batch_c_values);
    // on coset: n values of \sum{eta^i * c} on coset
    domain.coset_fft_in_place(&mut batch_c_values);

    // on coset: n values of \sum{eta^i * (ab - c)} on coset
    for j in 0..domain_size {
        batch_coset_ab[j] -= batch_c_values[j];
    }
    // cfg_iter_mut!(sum_coset_ab)
    //     .zip(sum_c)
    //     .for_each(|(sum_coset_ab_j, sum_coset_c_j)| *sum_coset_ab_j -= &sum_coset_c_j);

    domain.divide_by_vanishing_poly_on_coset_in_place(&mut batch_coset_ab);
    domain.coset_ifft_in_place(&mut batch_coset_ab);

    //abci_fft_time += start.elapsed();
    //println!("abci_fft_time: {:?}", abci_fft_time);

    let q_poly = DensePolynomial::from_coefficients_vec(batch_coset_ab);

    // Commit to quotient polynomial
    //let start2 = Instant::now();

    let (q_comm, q_rand) = KZG10::<E>::commit(&kzg10_ck, &q_poly, hiding_bound, Some(rng)).unwrap();
    let q_r_value = q_poly.evaluate(&r);
    let q_r_proof = KZG10::<E>::open(&kzg10_ck, &q_poly, r, &q_rand).unwrap();

    let proof = Proof {
        q_comm,
        q_r_value,
        q_r_proof,
    };

    //q_commit_time += start2.elapsed();
    //println!("q_commit_time: {:?}", q_commit_time);

    // let mut q_comm_bytes = vec![];
    // q_comm.write(&mut q_comm_bytes)?;
    // transcript.append_message(b"quotient polynomial commitments", &q_comm_bytes);

    // Prove
    // Generate a challenge
    // let mut c = [0u8; 31];
    // transcript.challenge_bytes(b"random point", &mut c);
    // let zeta = E::Fr::from_random_bytes(&c).unwrap();

    // r_q_polys.push(q_poly);
    // r_mid_q_rands.push(q_rand);

    //let mut open_r_mid_q_time = Duration::new(0, 0);
    //let start = Instant::now();

    // for j in 0..(m_mid + 1) {
    //     let value = r_q_polys[j + m_io].evaluate(&zeta);
    //     r_mid_q_values.push(value);
    // }

    // let opening_challenge = E::Fr::rand(rng);
    // let r_mid_q_proof = KZG10::<E>::batch_open(
    //     &kzg10_ck,
    //     &r_q_polys[m_io..],
    //     zeta,
    //     opening_challenge,
    //     &r_mid_q_rands,
    // )?;

    //open_r_mid_q_time += start.elapsed();
    //println!("open_r_mid_q_time: {:?}", open_r_mid_q_time);

    // let proof = Proof {
    //     r_mid_comms,
    //     q_comm,
    //     r_mid_q_values,
    //     r_mid_q_proof,
    //     opening_challenge,
    // };

    Ok(proof)
}
