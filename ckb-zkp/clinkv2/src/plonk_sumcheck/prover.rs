use std::ops::Neg;
use std::vec;

use ark_ec::PairingEngine;
use ark_ff::{FftField, Field, One, ToBytes, UniformRand, Zero};
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain, Radix2EvaluationDomain, MixedRadixEvaluationDomain, Polynomial, UVPolynomial, domain};
use ark_serialize::CanonicalSerialize;
use ark_std::{cfg_iter, cfg_iter_mut};
use rand::Rng;
use ark_std::rand::RngCore;
use rand::prelude::*;

use crate::batch_kzg10::batch_open;
// DEV
//use std::time::{Duration, Instant};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::plonk_sumcheck::{Proof, Kzg10ComKey};
use ark_poly_commit::kzg10::KZG10;

pub fn create_bgin19_proof<E: PairingEngine, R: RngCore>(
    inputs: Vec<Vec<Vec<E::Fr>>>,
    kzg10_ck: &Kzg10ComKey<'_, E>,
    L: usize,
    M: usize,
    rng: &mut R,
) -> Proof<E> {
    
    // inputs: L × M × 6
    // inputs[0], ..., inputs[L-1], each of M elements
    assert_eq!(inputs.len(), 6);
    assert_eq!(inputs[0].len(), L);
    assert_eq!(inputs[0][0].len(), M);

    let one = E::Fr::one();
    let zero = E::Fr::zero();

    let domain_M: GeneralEvaluationDomain<E::Fr> =
        EvaluationDomain::<E::Fr>::new(M).unwrap();
    let domain_M_size = domain_M.size();
    println!("domain_M size: {}", domain_M_size);

    // Compute q(x)

    // Use FFT, 6mlogM
    let mut f_polys: Vec<Vec<DensePolynomial<E::Fr>>> = vec![];
    for k in 0..6 {
        let mut f_k_polys: Vec<DensePolynomial<E::Fr>> = vec![];
        for l in 0..L {
            let f_lk_coeffs = domain_M.ifft(&inputs[k][l]);
            let f_lk_poly = DensePolynomial::from_coefficients_vec(f_lk_coeffs);
            f_k_polys.push(f_lk_poly);
        }
        f_polys.push(f_k_polys);
    }

    let mut f14_cos_pts = vec![zero; M];
    let mut f_l5_coeffs = vec![zero; M];
    let mut f_l6_coeffs = vec![zero; M];
    let eta = E::Fr::rand(rng);
    let mut eta_power = one;
    let mut eta_powers: Vec<E::Fr> = vec![];
    for l in 0..L {
        let f_l1_cos_pts = domain_M.coset_fft(&f_polys[0][l]);
        let f_l2_cos_pts = domain_M.coset_fft(&f_polys[1][l]);
        let f_l3_cos_pts = domain_M.coset_fft(&f_polys[2][l]);
        let f_l4_cos_pts = domain_M.coset_fft(&f_polys[3][l]);
        for m in 0..M {
            f14_cos_pts[m] += eta_power * (f_l1_cos_pts[m] * (f_l3_cos_pts[m] + f_l4_cos_pts[m]) + f_l2_cos_pts[m] * f_l3_cos_pts[m]);
            f_l5_coeffs[m] += eta_power * f_polys[4][l].coeffs[m];
            f_l6_coeffs[m] += eta_power * f_polys[5][l].coeffs[m];
        }
        eta_powers.push(eta_power);
        eta_power *= eta;
    }
    let f_l5_cos_pts = domain_M.coset_fft(&f_l5_coeffs);
    let f_l6_cos_pts = domain_M.coset_fft(&f_l6_coeffs);

    for m in 0..M {
        f14_cos_pts[m] = f14_cos_pts[m] + f_l5_cos_pts[m] - f_l6_cos_pts[m];
    }
    domain_M.divide_by_vanishing_poly_on_coset_in_place(&mut f14_cos_pts);
    // Perform iFFT on coset
    let q_coeffs = domain_M.coset_ifft(&f14_cos_pts);
    let q_poly = DensePolynomial::from_coefficients_vec(q_coeffs);
    // println!("q_poly degree: {:?}", q_poly.degree());


    // Use points, 6m mul, 1 MlogM
    // let h = E::Fr::multiplicative_generator(); // coset_offset
    // let h_squre = h * h;
    // let eta = E::Fr::rand(rng);
    // let mut eta_power = one;
    // let mut eta_powers: Vec<E::Fr> = vec![];
    // let mut q_points: Vec<E::Fr> = vec![zero; M];
    // for l in 0..L {
    //     for m in 0..M {
    //         q_points[m] += eta_power * (h_squre * (inputs[0][l][m] * (inputs[2][l][m] + inputs[3][l][m]) + inputs[1][l][m] * inputs[2][l][m]) + h * (inputs[4][l][m] - inputs[5][l][m]));
    //     }
    //     eta_powers.push(eta_power);
    //     eta_power *= eta;
    // }
    // domain.divide_by_vanishing_poly_on_coset_in_place(&mut q_points);
    // // Perform iFFT on coset
    // let q_coeffs = domain.coset_ifft(&q_points);
    // let q_poly = DensePolynomial::from_coefficients_vec(q_coeffs);
    // println!("q_poly degree: {:?}", q_poly.degree());

    // Commit to q(x)
    let hiding_bound = Some(1);
    let (q_comm, q_rand) = KZG10::<E, DensePolynomial<E::Fr>>::commit(&kzg10_ck, &q_poly, hiding_bound, Some(rng)).unwrap();

    let mut poly_comms = vec![];
    let mut poly_rands = vec![];

    // open q(r)
    let r = domain_M.sample_element_outside_domain(rng); 
    let q_r_value = q_poly.evaluate(&r);
    let q_r_proof = KZG10::<E, DensePolynomial<E::Fr>>::open(&kzg10_ck, &q_poly, r, &q_rand).unwrap();
    poly_comms.push(q_comm);
    poly_rands.push(q_rand);

    let mut open_values = vec![];
    let mut open_proofs = vec![];
    open_values.push(q_r_value);
    open_proofs.push(q_r_proof);

    // Compute p(x) = F1(x) * F4(x) + F1(x) * F3(x) + F2(x) * F3(x) + F4(x) - F5(x)
    // Compute flk(r), 6m mul
    let lag_r_values = domain_M.evaluate_all_lagrange_coefficients(r);
    let mut F_points: Vec<Vec<E::Fr>> = vec![]; // 6 * L
    for k in 0..6 {
        if k == 2 || k == 3 {
            F_points.push((0..L).map(|l| 
                (0..M).map(|m|
                    inputs[k][l][m] * lag_r_values[m] 
                ).sum()
            ).collect());
        } else {
            F_points.push((0..L).map(|l| 
                (0..M).map(|m|
                    eta_powers[m] * inputs[k][l][m] * lag_r_values[m] 
                ).sum()
            ).collect());
        }
    }      
    // println!("F_points[0].len(): {:?}", F_points[0].len());

    let domain_L: GeneralEvaluationDomain<E::Fr> =
        EvaluationDomain::<E::Fr>::new(L).unwrap();
    // Compute Fk(x), 6LlogL
    let mut F_polys = vec![];
    for k in 0..6 {
        let F_k_coeffs = domain_L.ifft(&F_points[k]);
        let F_k_poly = DensePolynomial::from_coefficients_vec(F_k_coeffs);
        F_polys.push(F_k_poly);
    }
    // println!("F_polys[0].degree(): {:?}", F_polys[0].degree());

    // Compute p(x), 2 poly mul, 6LlogL
    let mut p_poly = &F_polys[0] * &(&F_polys[2] + &F_polys[3]) + &F_polys[1] * &F_polys[2];
    p_poly += &F_polys[4];
    p_poly -= &F_polys[5];
    let vanishing_poly = domain_L.vanishing_polynomial();
    let vanishing_r_value = vanishing_poly.evaluate(&r);
    let s = q_r_value * vanishing_r_value;
    let v = s * domain_L.size_inv();
    p_poly -= &DensePolynomial::from_coefficients_vec(vec![v]);

    // println!("q_poly: {:?}", q_poly);
    // println!("p_poly: {:?}", p_poly);

    // Compute U(x)
    let mut U_points: Vec<E::Fr> = vec![];
    let U_1 = E::Fr::rand(rng);
    let mut U_i_value = U_1;
    U_points.push(U_1);
    for l in 0..L {
        let p_l_value = F_points[0][l] * (F_points[2][l] + F_points[2][l]) + F_points[1][l] * F_points[2][l] + F_points[4][l] - F_points[5][l];
        U_i_value += p_l_value;
        U_points.push(U_i_value);
    }

    let U_poly = DensePolynomial::from_coefficients_vec(domain_L.ifft(&U_points));

    // let mut U_g_points = U_points[1..].to_vec();
    // U_g_points.push(U_1);
    // let U_g_poly= DensePolynomial::from_coefficients_vec(domain_L.ifft(&U_g_points));
    let g = domain_L.group_gen();
    let mut U_g_coeffs = U_poly.coeffs().clone().to_vec();
    GeneralEvaluationDomain::<<E as PairingEngine>::Fr>::distribute_powers(&mut U_g_coeffs, g);
    let U_g_poly = DensePolynomial::from_coefficients_vec(U_g_coeffs);

    // Compute h(x)
    let mut Qt_poly = &U_g_poly - &U_poly.clone();
    Qt_poly -= &p_poly.clone();
    let mut Qt_coset_points = domain_L.coset_fft(Qt_poly.coeffs());
    domain_L.divide_by_vanishing_poly_on_coset_in_place(&mut Qt_coset_points);
    domain_L.coset_ifft_in_place(&mut Qt_coset_points);
    let Q_poly = DensePolynomial::from_coefficients_vec(Qt_coset_points);    

    // Commit to U(x), h(x), p(x)
    let (p_comm, p_rand) = KZG10::<E, DensePolynomial<E::Fr>>::commit(&kzg10_ck, &p_poly, hiding_bound, Some(rng)).unwrap();
    let (U_comm, U_rand) = KZG10::<E, DensePolynomial<E::Fr>>::commit(&kzg10_ck, &U_poly, hiding_bound, Some(rng)).unwrap();
    let (Q_comm, Q_rand) = KZG10::<E, DensePolynomial<E::Fr>>::commit(&kzg10_ck, &Q_poly, hiding_bound, Some(rng)).unwrap();

    let z = domain_L.sample_element_outside_domain(rng);
    let gz = g * z;
    open_values.push(p_poly.evaluate(&z));
    open_values.push(U_poly.evaluate(&z));
    open_values.push(Q_poly.evaluate(&z));
    open_values.push(U_poly.evaluate(&gz));

    open_proofs.push(KZG10::<E, DensePolynomial<E::Fr>>::open(&kzg10_ck, &U_poly, gz, &U_rand).unwrap());

    poly_comms.push(p_comm);
    poly_comms.push(U_comm);
    poly_comms.push(Q_comm);
    poly_rands.push(p_rand);
    poly_rands.push(U_rand);
    poly_rands.push(Q_rand);

    let mut polys = vec![];
    println!("p_poly degree: {:?}", p_poly.degree());
    println!("U_poly degree: {:?}", U_poly.degree());
    println!("h_poly degree: {:?}", Q_poly.degree());
    polys.push(p_poly);
    polys.push(U_poly);
    polys.push(Q_poly);

    let open_challenge = E::Fr::rand(rng);
    open_proofs.push(batch_open::<E>(&kzg10_ck, &polys[1..], z, open_challenge, &poly_rands[1..].to_vec()).unwrap());

    Proof { poly_comms, open_values, open_proofs, open_challenge }
}
