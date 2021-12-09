#![allow(non_snake_case)]

use std::vec;
use ark_ec::PairingEngine;
use ark_ff::{One, UniformRand, Zero};
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain, Polynomial, UVPolynomial};
use ark_std::rand::RngCore;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::fliop::Proof;

pub fn create_bgin19_proof<E: PairingEngine, R: RngCore>(
    inputs: Vec<Vec<E::Fr>>,
    rs: &Vec<E::Fr>,
    ws: &Vec<E::Fr>,
    rng: &mut R,
) -> (Proof<E>, Proof<E>) {
    
    let m = inputs[0].len();
    // Make sure that m is a power of 2
    assert_eq!(m, m.next_power_of_two());

    // f(x) degree: 1, need 2 points to reconstruct
    let domain_2: GeneralEvaluationDomain<E::Fr> =
        EvaluationDomain::<E::Fr>::new(2).unwrap();

    // q(x) degree: 2, need 3 points to reconstruct, should work on domain of size 4
    let domain_4: GeneralEvaluationDomain<E::Fr> =
        EvaluationDomain::<E::Fr>::new(4).unwrap();

    let g_zero_domain4 = E::Fr::one();
    let g_one_domain4 = domain_4.group_gen();
    let g_two_domain4 = g_one_domain4 * g_one_domain4;

    let mut L = m/2;
    let mut cnt = 0;

    let mut f_r_values = inputs.clone();

    let mut q_shares_1: Vec<Vec<E::Fr>> = vec![];
    let mut q_shares_2: Vec<Vec<E::Fr>> = vec![];

    loop {

        println!("L={}:", L);

        // the last round
        if L == 1 {
            // f(x) degree: 2, should work over domain 4
            // q(x) degree: 4, need 5 points to reconstruct, thus should work over domain 8
            let domain_8: GeneralEvaluationDomain<E::Fr> =
                EvaluationDomain::<E::Fr>::new(8).unwrap();
            assert_eq!(f_r_values[0].len(), 2);
            let mut q_L_poly = DensePolynomial::<E::Fr>::zero();
            // 5 polys of degree 2
            let mut f_polys = vec![];
            for k in 0..5 {
                let mut points = vec![ws[k]; 1];
                points.append(&mut f_r_values[k]); // length 3
                // println!("points.len(): {}", points.len());
                let fk_coeffs = domain_4.ifft(&points);
                // println!("fk_coeffs.len(): {}", fk_coeffs.len());
                let fk_poly = DensePolynomial::from_coefficients_vec(fk_coeffs);
                f_polys.push(fk_poly);
            }
            // shoud be degree-2, 3 after padding
            println!("f_polys[0].coeffs.len(): {:?}", f_polys[0].coeffs.len());
            assert_eq!(f_polys[0].coeffs.len(), 4);

            // shoud be degree-4, 6 after padding
            q_L_poly += &(&f_polys[0] * &(&f_polys[2] + &f_polys[3]) + &f_polys[1] * &f_polys[2]);
            q_L_poly += &f_polys[4];
            println!("q_polys[{:?}]: {:?}", cnt, &q_L_poly.coeffs());
            println!("q_l_coeffs.len():{:?}", q_L_poly.len());
            assert_eq!(q_L_poly.len(), 7);

            // Need 7 points to recover q(x)
            // let g_domain8 = domain_8.group_gen();
            // let mut q_L_values = vec![];
            // let mut g_i = E::Fr::one();
            // for _ in 0..7 {
            //     q_L_values.push(q_L_poly.evaluate(&g_i));
            //     g_i *= g_domain8;
            // }
            // println!("rs[{:?}]: {:?}", cnt, &rs[cnt]);
            // println!("q_r_value[{:?}]: {:?}", cnt, q_L_poly.evaluate(&rs[cnt]));
            let q_L_values = q_L_poly.evaluate_over_domain_by_ref(domain_8);

            let q_L_shares_1 = (0..8).map(|_| E::Fr::rand(rng) ).collect::<Vec<E::Fr>>();
            let q_L_shares_2 = (0..8).map(|i| q_L_values[i] - q_L_shares_1[i] ).collect::<Vec<E::Fr>>();
            
            q_shares_1.push(q_L_shares_1);
            q_shares_2.push(q_L_shares_2);

            break;
        }

        // Common rounds
        let mut q_L_poly = DensePolynomial::<E::Fr>::zero();
        let mut f_polys = vec![];
        for k in 0..5 {
            let mut fk_polys: Vec<DensePolynomial<E::Fr>> = vec![];
            for l in 0..L {
                // length of f_r_values: 2L (m -> 4)
                let points = [f_r_values[k][l], f_r_values[k][l+L]];
                let flk_coeffs = domain_2.ifft(&points);
                let flk_poly = DensePolynomial::from_coefficients_vec(flk_coeffs);
                // L polys of degree 1
                fk_polys.push(flk_poly);
            }
            // 5L polys of degree 1
            f_polys.push(fk_polys);
        }
        // shoud be degree-1
        println!("f_polys[0][0].coeffs.len(): {:?}", f_polys[0][0].coeffs.len());
        assert_eq!(f_polys[0][0].coeffs.len(), 2);

        // degree of q(x): 2
        for l in 0..L {
            q_L_poly += &(&f_polys[0][l] * &(&f_polys[2][l] + &f_polys[3][l]) + &f_polys[1][l] * &f_polys[2][l]);
            q_L_poly += &f_polys[4][l];
        }

        // For debug
        println!("q_polys[{:?}]: {:?}", cnt, &q_L_poly.coeffs());

        // shoud be degree-2
        println!("q_l_coeffs.len():{:?}", q_L_poly.len());
        assert_eq!(q_L_poly.coeffs.len(), 3);

        // Need 3 points to recover q(x), evaluate over domain of size 4
        // let mut q_L_values = vec![];
        // q_L_values.push(q_L_poly.evaluate(&g_zero_domain4));
        // q_L_values.push(q_L_poly.evaluate(&g_one_domain4));
        // q_L_values.push(q_L_poly.evaluate(&g_two_domain4));

        // IMPORTANT: verifier uses ifft to recover q(x), interpolated values are padded with zeros
        // however q(x) doesn't actually evaluate to zeros at the padded points, since it is calculated by poly-mul
        let q_L_values = q_L_poly.evaluate_over_domain_by_ref(domain_4);

        let q_L_shares_1 = (0..4).map(|_| E::Fr::rand(rng) ).collect::<Vec<E::Fr>>();
        let q_L_shares_2 = (0..4).map(|i| q_L_values[i] - q_L_shares_1[i] ).collect::<Vec<E::Fr>>();
        
        q_shares_1.push(q_L_shares_1);
        q_shares_2.push(q_L_shares_2);

        let r = rs[cnt]; // domain_3.sample_element_outside_domain(rng);
        
        // For debug
        // println!("rs[{:?}]: {:?}", cnt, &rs[cnt]);
        // println!("q_r_value[{:?}]: {:?}", cnt, q_L_poly.evaluate(&r));

        cnt += 1; 
        f_r_values = (0..5).map(|k|
            (0..L).map(|l|
                f_polys[k][l].evaluate(&r)
            ).collect::<Vec<_>>()
        ).collect::<Vec<Vec<_>>>();

        L /= 2;
    }

    let proof1 = Proof {
        q_shares: q_shares_1,
        // c_shares: c_shares_1,
    };
    let proof2 = Proof {
        q_shares: q_shares_2,
        // c_shares: c_shares_2,
    };
    
    (proof1, proof2)
}
