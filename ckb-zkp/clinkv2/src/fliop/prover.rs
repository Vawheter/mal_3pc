#![allow(non_snake_case)]

use std::vec;
use ark_ec::PairingEngine;
use ark_ff::{UniformRand, Zero, One};
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
use ark_std::rand::RngCore;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::fliop::Proof;

pub fn create_bgin19_proof<E: PairingEngine, R: RngCore>(
    inputs: Vec<Vec<E::Fr>>,
    rs: &Vec<E::Fr>,
    ws: &Vec<E::Fr>,
    theta: E::Fr,
    rng: &mut R,
) -> (Proof<E>, Proof<E>) {
    
    let m = inputs[0].len();
    // Make sure that m is a power of 2
    assert_eq!(m, m.next_power_of_two());

    let zero = E::Fr::zero();
    let one = E::Fr::one();

    // Common rounds
    // f(x) degree: 1, should work over domain 2
    let domain2: GeneralEvaluationDomain<E::Fr> =
        EvaluationDomain::<E::Fr>::new(2).unwrap();
    // let g_domain2 = domain2.group_gen();

    // q(x) degree: 2, should work over domain 4, need 3 points to reconstruct
    let domain4: GeneralEvaluationDomain<E::Fr> =
        EvaluationDomain::<E::Fr>::new(4).unwrap();

    // Wish to compute q(1), q(g), q(g^2) where g = domain4.group_gen() and q(1) and q(g) where g = domain2.group_gen()
    let mut lag_g4_g2_over_domain2 = vec![];
    let g_domain4 = domain4.group_gen();
    let mut g4_power = one;
    for _ in 0..4 {
        lag_g4_g2_over_domain2.push(domain2.evaluate_all_lagrange_coefficients(g4_power));
        g4_power *= g_domain4;
    }
    // lag_g4_g2_over_domain2.push(domain2.evaluate_all_lagrange_coefficients(g_domain2));

    // The last round
    // f(x) degree: 2, shoud work over domain 4
    // q(x) degree: 4, shoud work over domain 8, need 5 points to reconstruct
    let domain8: GeneralEvaluationDomain<E::Fr> =
        EvaluationDomain::<E::Fr>::new(8).unwrap();

    // Wish to compute q(1), q(g), ..., q(g^4) where g = domain8.group_gen() and q(1), q(g) where g = domain4.group_gen()
    let mut lag_g8_g4_over_domain4 = vec![];
    let g_domain8 = domain8.group_gen();
    let mut g8_power = one;
    for _ in 0..8 {
        lag_g8_g4_over_domain4.push(domain4.evaluate_all_lagrange_coefficients(g8_power));
        g8_power *= g_domain8;
    }
    // lag_g8_g4_over_domain4.push(domain4.evaluate_all_lagrange_coefficients(g_domain4));

    let mut L = m/2;
    let mut cnt = 0;

    let mut f_r_values = inputs.clone();
    // xi, xi-1, yi, yi-1, alphai, zi
    // Linear shift: theta * xi, theta * xi-1, theta * alphai, theta * zi
    // shift index: 0, 1, 4, no need for zi
    let mut theta_power = one;
    for j in 0..m {
        f_r_values[0][j] *= theta_power;
        f_r_values[1][j] *= theta_power;
        f_r_values[4][j] *= theta_power;
        f_r_values[5][j] *= theta_power;
        theta_power *= theta;
    }

    let mut q_shares_1: Vec<Vec<E::Fr>> = vec![];
    let mut q_shares_2: Vec<Vec<E::Fr>> = vec![];

    loop {
        // the last round
        if L == 1 {
            // assert_eq!(f_r_values[0].len(), 2);

            // q(1), q(g), ..., q(g^4) where g = domain8.group_gen() and q(1), q(g) where g = domain4.group_gen()
            let num_points = 8;
            let mut q_g8_g4_values = vec![zero; num_points];
            let mut f_g8_values = vec![];
            for k in 0..5 {
                let fk_g8_value = (0..num_points).map(|i|
                    ws[k] * lag_g8_g4_over_domain4[i][0] + f_r_values[k][0] * lag_g8_g4_over_domain4[i][1] + f_r_values[k][1] * lag_g8_g4_over_domain4[i][2]
                ).collect::<Vec<_>>();
                f_g8_values.push(fk_g8_value);
            }

            for i in 0..num_points {
                q_g8_g4_values[i] = f_g8_values[0][i] * (f_g8_values[2][i] + f_g8_values[3][i]) + f_g8_values[1][i] * f_g8_values[2][i] + f_g8_values[4][i];
            }

            let g1: E::Fr = f_r_values[0][0] * (f_r_values[2][0] + f_r_values[3][0]) + f_r_values[1][0] * f_r_values[2][0] + f_r_values[4][0];
            let g2: E::Fr = f_r_values[0][1] * (f_r_values[2][1] + f_r_values[3][1]) + f_r_values[1][1] * f_r_values[2][1] + f_r_values[4][1];

            q_g8_g4_values.push(g1);
            q_g8_g4_values.push(g2);

            let q_L_shares_1 = (0..num_points + 2).map(|_| E::Fr::rand(rng) ).collect::<Vec<E::Fr>>();
            let q_L_shares_2 = (0..num_points + 2).map(|i| q_g8_g4_values[i] - q_L_shares_1[i] ).collect::<Vec<E::Fr>>();
            
            q_shares_1.push(q_L_shares_1);
            q_shares_2.push(q_L_shares_2);

            break;
        }

        // Common rounds

        // q(1), q(g), q(g^2) where g = domain4.group_gen() and q(1) and q(g) where g = domain2.group_gen()
        let num_points = 4;
        let mut q_g4_g2_values = vec![zero; num_points];
        let mut f_g4_values = vec![];
        for k in 0..5 {
            let mut fk_g4_values = vec![];
            // length of f_r_values: 2L (m -> 4)
            for l in 0..L {
                let fkl_g4_values = (0..num_points).map(|i|
                    f_r_values[k][l] * lag_g4_g2_over_domain2[i][0] + f_r_values[k][l+L] * lag_g4_g2_over_domain2[i][1]
                ).collect::<Vec<_>>();
                fk_g4_values.push(fkl_g4_values);
            }
            f_g4_values.push(fk_g4_values);
        }

        for i in 0..num_points {
            for l in 0..L {
                q_g4_g2_values[i] += f_g4_values[0][l][i] * (f_g4_values[2][l][i] + f_g4_values[3][l][i]) + f_g4_values[1][l][i] * f_g4_values[2][l][i] + f_g4_values[4][l][i];
            }
        }

        let g1: E::Fr = (0..L).map(|l|
            f_r_values[0][l] * (f_r_values[2][l] + f_r_values[3][l]) + f_r_values[1][l] * f_r_values[2][l] + f_r_values[4][l]
        ).sum();
        let g2: E::Fr = (L..2 * L).map(|l|
            f_r_values[0][l] * (f_r_values[2][l] + f_r_values[3][l]) + f_r_values[1][l] * f_r_values[2][l] + f_r_values[4][l]
        ).sum();

        q_g4_g2_values.push(g1);
        q_g4_g2_values.push(g2);

        let q_L_shares_1 = (0..num_points + 2).map(|_| E::Fr::rand(rng) ).collect::<Vec<E::Fr>>();
        let q_L_shares_2 = (0..num_points + 2).map(|i| q_g4_g2_values[i] - q_L_shares_1[i] ).collect::<Vec<E::Fr>>();
        
        q_shares_1.push(q_L_shares_1);
        q_shares_2.push(q_L_shares_2);

        let r = rs[cnt]; // domain3.sample_element_outside_domain(rng);
        cnt += 1;

        let lar_r_over_domain2 = domain2.evaluate_all_lagrange_coefficients(r);
        f_r_values = (0..6).map(|k|
            (0..L).map(|l|
                f_r_values[k][l] * lar_r_over_domain2[0] + f_r_values[k][l+L] * lar_r_over_domain2[1]
            ).collect::<Vec<_>>()
        ).collect::<Vec<Vec<_>>>();

        L /= 2;
    }

    let proof1 = Proof {
        q_shares: q_shares_1,
    };
    let proof2 = Proof {
        q_shares: q_shares_2,
    };
    
    (proof1, proof2)
}
