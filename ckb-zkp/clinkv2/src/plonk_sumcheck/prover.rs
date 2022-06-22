#![allow(non_snake_case)]

use std::vec;

use ark_ec::PairingEngine;
use ark_ff::{One, UniformRand, Zero};
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain, Polynomial, UVPolynomial};
use ark_std::rand::RngCore;

use crate::batch_kzg10::batch_open;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::plonk_sumcheck::{Proof, Kzg10ComKey};
use ark_poly_commit::kzg10::KZG10;

pub fn create_bgin19_proof<E: PairingEngine, R: RngCore>(
    inputs: Vec<Vec<Vec<E::Fr>>>,
    kzg10_ck: &Kzg10ComKey<'_, E>,
    eta: E::Fr,
    r: E::Fr,
    z: E::Fr,
    rng: &mut R,
) -> Proof<E> {
    
    // inputs: L × M × 6
    // inputs[0], ..., inputs[L-1], each of M elements
    assert_eq!(inputs.len(), 6);
    let L = inputs[0].len();
    let M = inputs[0][0].len();

    let one = E::Fr::one();
    let zero = E::Fr::zero();

    let domain_M: GeneralEvaluationDomain<E::Fr> =
        EvaluationDomain::<E::Fr>::new(M).unwrap();
    let domain_M_size = domain_M.size();
    // println!("domain_M size: {}", domain_M_size);

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

    // IMPORTANT: number of coset points should be domain_M_size
    let mut f14_cos_pts = vec![zero; domain_M_size];
    let mut f_l5_coeffs = vec![zero; domain_M_size];
    let mut f_l6_coeffs = vec![zero; domain_M_size];
    // let eta = E::Fr::rand(rng);
    let mut eta_power = one;
    let mut eta_powers: Vec<E::Fr> = vec![];
    for l in 0..L {
        eta_powers.push(eta_power);
        let f_l1_cos_pts = domain_M.coset_fft(&f_polys[0][l]);
        let f_l2_cos_pts = domain_M.coset_fft(&f_polys[1][l]);
        let f_l3_cos_pts = domain_M.coset_fft(&f_polys[2][l]);
        let f_l4_cos_pts = domain_M.coset_fft(&f_polys[3][l]);
        for m in 0..domain_M_size {
            f14_cos_pts[m] += eta_power * (f_l1_cos_pts[m] * (f_l3_cos_pts[m] + f_l4_cos_pts[m]) + f_l2_cos_pts[m] * f_l3_cos_pts[m]);
            f_l5_coeffs[m] += eta_power * f_polys[4][l].coeffs[m];
            f_l6_coeffs[m] += eta_power * f_polys[5][l].coeffs[m];
        }
        eta_power *= eta;
    }
    let f_l5_cos_pts = domain_M.coset_fft(&f_l5_coeffs);
    let f_l6_cos_pts = domain_M.coset_fft(&f_l6_coeffs);

    for m in 0..domain_M_size {
        f14_cos_pts[m] += f_l5_cos_pts[m] - f_l6_cos_pts[m];
    }
    domain_M.divide_by_vanishing_poly_on_coset_in_place(&mut f14_cos_pts);
    // Perform iFFT on coset
    let q_coeffs = domain_M.coset_ifft(&f14_cos_pts);
    let q_poly = DensePolynomial::from_coefficients_vec(q_coeffs);
    // println!("q_poly degree: {:?}", q_poly.degree());

    // Commit to q(x)
    let hiding_bound = Some(1);
    let (q_comm, q_rand) = KZG10::<E, DensePolynomial<E::Fr>>::commit(&kzg10_ck, &q_poly, hiding_bound, Some(rng)).unwrap();

    let mut poly_comms = vec![];

    // open q(r)
    let q_r_value = q_poly.evaluate(&r);
    let q_r_proof = KZG10::<E, DensePolynomial<E::Fr>>::open(&kzg10_ck, &q_poly, r, &q_rand).unwrap();
    poly_comms.push(q_comm);

    let mut open_values = vec![];
    let mut open_proofs = vec![];
    open_values.push(q_r_value);
    open_proofs.push(q_r_proof);

    // Compute p(x) = F1(x) * F4(x) + F1(x) * F3(x) + F2(x) * F3(x) + F4(x) - F5(x)
    // Compute flk(r), 6m mul
    // let lag_r_values = domain_M.evaluate_all_lagrange_coefficients(r);
    let mut F_points: Vec<Vec<E::Fr>> = vec![]; // 6 * L
    for k in 0..6 {
        if k == 2 || k == 3 {
            F_points.push((0..L).map(|l| 
                f_polys[k][l].evaluate(&r)
            ).collect());
        } else {
            F_points.push((0..L).map(|l| 
                eta_powers[l] * f_polys[k][l].evaluate(&r)
            ).collect());
        }
    }      

    let vanishing_r_value = domain_M.evaluate_vanishing_polynomial(r);
    let s = q_r_value * vanishing_r_value;
    
    // // Check q(x) correctness: OK
    // let mut f_r_value = zero;
    // for l in 0..L {
    //     f_r_value += F_points[0][l] * (F_points[2][l] + F_points[3][l]) + F_points[1][l] * F_points[2][l] + F_points[4][l] - F_points[5][l]
    // }
    // assert_eq!(f_r_value, s);

    let domain_L: GeneralEvaluationDomain<E::Fr> =
        EvaluationDomain::<E::Fr>::new(L).unwrap();
    let domain_L_size = domain_L.size();

    // Compute Fk(x), 6LlogL
    let mut F_polys = vec![];
    for k in 0..6 {
        let F_k_coeffs = domain_L.ifft(&F_points[k]);
        let F_k_poly = DensePolynomial::from_coefficients_vec(F_k_coeffs);
        F_polys.push(F_k_poly);
    }
    // println!("F_polys[0].degree(): {:?}", F_polys[0].degree());

    // Compute P(x), 2 poly mul, 6LlogL
    let domain_2L: GeneralEvaluationDomain<E::Fr> =
        EvaluationDomain::<E::Fr>::new(2 * L - 1).unwrap();
    let domain_2L_size = domain_2L.size();
    let mut p_2L_values = vec![zero; domain_2L_size];
    let F_2L_values = (0..6).map(|k|
        F_polys[k].evaluate_over_domain_by_ref(domain_2L).evals
    ).collect::<Vec<_>>();

    for l in 0..domain_2L_size {
        p_2L_values[l] += F_2L_values[0][l] * (F_2L_values[2][l] + F_2L_values[3][l]) + F_2L_values[1][l] * F_2L_values[2][l] + F_2L_values[4][l] - F_2L_values[5][l]
    }
    let p_coeffs = domain_2L.ifft(&p_2L_values);
    let mut p_poly = DensePolynomial::from_coefficients_vec(p_coeffs);

    // let mut p_poly = &F_polys[0] * &(&F_polys[2] + &F_polys[3]) + &F_polys[1] * &F_polys[2];
    // p_poly += &F_polys[4];
    // p_poly -= &F_polys[5];

    // Check p(x) correctness: OK
    // assert_eq!(p_poly.evaluate(&z), F_polys[0].evaluate(&z) * (F_polys[2].evaluate(&z) + F_polys[3].evaluate(&z)) + F_polys[1].evaluate(&z) * F_polys[2].evaluate(&z) + F_polys[4].evaluate(&z) - F_polys[5].evaluate(&z));
    
    // IMPORTANT: p(x) are now of degree 2*(L-1), but only need to evaluate L points

    let v = s * domain_L.size_inv();
    p_poly.coeffs[0] -= v;

    let mut p_values = vec![];
    for l in 0..L {
        // IMPORTANT: the two polys perform fft on a larger domain of a larger domain which handles their product poly
        p_values.push(F_points[0][l] * (F_points[2][l] + F_points[3][l]) + F_points[1][l] * F_points[2][l] + F_points[4][l] - F_points[5][l] - v);
    }
    for _ in L..domain_L_size {
        p_values.push(-v);
    }

    // Check \sum p(x) = 0: OK
    // assert_eq!(p_values.iter().map(|val| *val).sum::<E::Fr>(), zero);

    // Compute U(x)
    let mut U_points: Vec<E::Fr> = vec![];
    let U_g0 = E::Fr::rand(rng);
    let mut U_gl_value = U_g0;
    U_points.push(U_g0);
    for l in 1..domain_L_size {
        // U(g^l) = U(g^0) + \sum_{i=0}^{l-1}p(g^i)
        // IMPORTANT: should add the constant term term s/domain_L_size 
        U_gl_value += p_values[l - 1];
        U_points.push(U_gl_value);
    }
    // U(g^L) should be equal to U(g^0)

    let U_poly = DensePolynomial::from_coefficients_vec(domain_L.ifft(&U_points));
    
    let g = domain_L.group_gen();
    let gz = g * z;

    // Check U(x) correctness: OK
    // for l in 0..domain_L_size {
    //     let x = domain_L.element(l);
    //     assert_eq!(U_poly.evaluate(&(g*x)), U_poly.evaluate(&x) + p_values[l]);
    // }

    // let mut U_g_points = U_points[1..].to_vec();
    // U_g_points.push(U_g0);
    // let U_g_poly= DensePolynomial::from_coefficients_vec(domain_L.ifft(&U_g_points));
    let mut U_g_coeffs = U_poly.coeffs().to_vec();
    GeneralEvaluationDomain::<<E as PairingEngine>::Fr>::distribute_powers(&mut U_g_coeffs, g);
    let U_g_poly = DensePolynomial::from_coefficients_vec(U_g_coeffs);

    // // Check U(gx) correctness: OK
    // assert_eq!(U_g_poly.evaluate(&z), U_poly.evaluate(&gz));

    // Compute Q(x)
    let mut Qt_poly = &U_g_poly - &U_poly;
    Qt_poly -= &p_poly;
    // IMPORTANT: coset_fft() would resize p(x) to its domain size
    // let mut Qt_coset_points = domain_L.coset_fft(Qt_poly.coeffs());
    // domain_L.divide_by_vanishing_poly_on_coset_in_place(&mut Qt_coset_points);
    // domain_L.coset_ifft_in_place(&mut Qt_coset_points);
    let (Q_poly, _remainder) = Qt_poly.divide_by_vanishing_poly(domain_L).unwrap();
    // Remainder should be zero
    // assert_eq!(Qt_poly, Q_poly.mul_by_vanishing_poly(domain_L) + remainder.clone());
    // assert_eq!(_remainder, DensePolynomial::zero());
    // let Q_poly = DensePolynomial::from_coefficients_vec(Qt_coset_points);    

    // Check Q(x) correctness: OK
    // assert_eq!(U_g_poly.evaluate(&z) - U_poly.evaluate(&z) - p_poly.evaluate(&z), Q_poly.evaluate(&z) * domain_L.evaluate_vanishing_polynomial(z));

    // Commit to U(x), Q(x)
    let (U_comm, U_rand) = KZG10::<E, DensePolynomial<E::Fr>>::commit(&kzg10_ck, &U_poly, hiding_bound, Some(rng)).unwrap();
    let (Q_comm, Q_rand) = KZG10::<E, DensePolynomial<E::Fr>>::commit(&kzg10_ck, &Q_poly, hiding_bound, Some(rng)).unwrap();

    open_values.push(U_poly.evaluate(&gz));
    open_values.push(U_poly.evaluate(&z));
    open_values.push(Q_poly.evaluate(&z));

    open_proofs.push(KZG10::<E, DensePolynomial<E::Fr>>::open(&kzg10_ck, &U_poly, gz, &U_rand).unwrap());

    poly_comms.push(U_comm);
    poly_comms.push(Q_comm);

    // U_rand, p_rand, Q_rand
    let mut poly_rands = vec![];
    poly_rands.push(U_rand);
    poly_rands.push(Q_rand);

    let mut polys = vec![];
    polys.push(U_poly);
    polys.push(Q_poly);

    let open_challenge = E::Fr::rand(rng);
    open_proofs.push(batch_open::<E>(&kzg10_ck, &polys[..], z, open_challenge, &poly_rands[..].to_vec()).unwrap());

    Proof { 
        poly_comms, // q(x), U(x), Q(x)
        open_values, // q(r), U(gz), U(z), Q(z)
        open_proofs, // q(r), U(gz), batch(U(z), Q(z))
        open_challenge 
    }
}
