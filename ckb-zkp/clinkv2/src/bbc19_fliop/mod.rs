#![allow(non_snake_case)]

use ark_ec::PairingEngine;
use ark_serialize::*;

pub mod prover;
pub mod verifier;

pub use prover::prove_dot_prod;
// pub use verifier::{verify_dot_prod, verify_dot_prod_helper};

pub type Kzg10Comm<E> = ark_poly_commit::kzg10::Commitment<E>;
pub type Kzg10Proof<E> = ark_poly_commit::kzg10::Proof<E>;
pub type Kzg10ComKey<'a, E> = ark_poly_commit::kzg10::Powers<'a, E>;

use crate::Vec;

/// The proof in PLONK FLIOP.
#[derive(Clone, Debug, Eq, PartialEq, CanonicalSerialize, CanonicalDeserialize)]
pub struct Proof<E: PairingEngine> {
    p_coeffs_shares1: Vec<Vec<E::Fr>>,
    p_coeffs_shares2: Vec<Vec<E::Fr>>,
    G_evals: Vec<E::Fr>,
}

/// The verification message in FLPCP.
#[derive(Clone, Debug, Eq, PartialEq, CanonicalSerialize, CanonicalDeserialize)]
pub struct VerMsg<E: PairingEngine> {
    F_z_shares: Vec<E::Fr>,
}