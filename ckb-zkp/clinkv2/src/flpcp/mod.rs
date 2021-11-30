use ark_ec::PairingEngine;
use ark_ff::Field;
use ark_serialize::*;

pub mod prover;
pub mod verifier;

pub use prover::create_bgin19_proof;
pub use verifier::{gen_vermsg, verify_bgin19_proof};

use crate::{String, Vec};

/// The proof in FLPCP.
#[derive(Clone, Debug, Eq, PartialEq, CanonicalSerialize, CanonicalDeserialize)]
pub struct Proof<E: PairingEngine> {
    ws: Vec<E::Fr>,
    p_coeffs: Vec<E::Fr>,
}

/// The verification message in FLPCP.
#[derive(Clone, Debug, Eq, PartialEq, CanonicalSerialize, CanonicalDeserialize)]
pub struct VerMsg<E: PairingEngine> {
    f_r_shares: Vec<E::Fr>,
    p_r_value: E::Fr,
    b: E::Fr,
}