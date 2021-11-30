#![allow(non_snake_case)]

use ark_ec::PairingEngine;
use ark_ff::{One, Zero};
use ark_poly::Polynomial;
use ark_poly::polynomial::univariate::DensePolynomial;

use crate::{
    kzg10_distributed::{Proof, VerifyAssignment, VerifyKey, KZG10},
    r1cs::{Index, SynthesisError},
};

pub fn gen_ver_msg<E: PairingEngine>(
    w_sec_polys: &Vec<DensePolynomial<E::Fr>>,
    r: E::Fr,
) -> Vec<E::Fr> {
    let mut w_sec_r_values = vec![];
    for i in 0..w_sec_polys.len() {
        w_sec_r_values.push(w_sec_polys[i].evaluate(&r));
    }
    w_sec_r_values
}

pub fn verify_proof<E: PairingEngine> (
    circuit: &VerifyAssignment<E>,
    kzg10_vk: &VerifyKey<E>,
    proof: &Proof<E>,
    r: E::Fr,
    // lag_r_values: Vec<E::Fr>,
    w_polys: &Vec<Vec<DensePolynomial<E::Fr>>>,
    vanishing_r_value: E::Fr,
    eta: E::Fr,
    // io: &Vec<Vec<E::Fr>>,
) -> Result<bool, SynthesisError> {

    let t_io = circuit.input_assignment.len();
    assert_eq!(t_io, 0);

    // Number of all variables
    let T = circuit.aux_assignment.len();
    // Number of constraints
    let K = circuit.at.len();
    let L = w_polys[0].len();

    let zero = E::Fr::zero();
    let one = E::Fr::one();

    assert!(KZG10::<E>::check(
        &kzg10_vk,
        &proof.q_comm,
        r,
        proof.q_r_value,
        &proof.q_r_proof,
    ).unwrap());

    // let domain: GeneralEvaluationDomain<E::Fr> =
    //     EvaluationDomain::<E::Fr>::new(M).ok_or(SynthesisError::PolynomialDegreeTooLarge)?;

    //let domain_size = domain.size();

    // T*L witness polys, thus yielding T*L values at point r
    let mut w_r_values = vec![];
    for l in 0..L {
        for t in 0..T {
            let w_lt_r_value = w_polys[t][l].evaluate(&r);
            w_r_values.push(w_lt_r_value);
        }
    }

    // let lag_values = domain.evaluate_all_lagrange_coefficients(zeta);

    // let vanishing_poly = domain.vanishing_polynomial();
    // let vanishing_value = vanishing_poly.evaluate(&zeta);

    let mut batch_ab_c = zero;
    let mut eta_lk = one;

    for k in 0..K {
        let mut a_lk_r = zero;
        for (coeff, index) in (&circuit.at[k]).into_iter() {
            let t = match index {
                Index::Input(j) => *j,
                Index::Aux(j) => *j,
            };
            a_lk_r += &(w_r_values[t] * coeff);
        }

        let mut b_lk_r = zero;
        for (coeff, index) in (&circuit.bt[k]).into_iter() {
            let t = match index {
                Index::Input(j) => *j,
                Index::Aux(j) => *j,
            };
            b_lk_r += &(w_r_values[t] * coeff);
        }

        let mut c_lk_r = zero;
        for (coeff, index) in (&circuit.ct[k]).into_iter() {
            let t = match index {
                Index::Input(j) => *j,
                Index::Aux(j) => *j,
            };
            c_lk_r += &(w_r_values[t] * coeff);
        }

        batch_ab_c += &(eta_lk * &(a_lk_r * &b_lk_r - &c_lk_r));
        eta_lk *= &eta;
    }
    assert_eq!(batch_ab_c, proof.q_r_value * &vanishing_r_value);

    Ok(true)
}
