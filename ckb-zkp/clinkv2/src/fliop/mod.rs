use ark_ec::PairingEngine;
use ark_serialize::*;

pub mod prover;
pub mod verifier;

pub use prover::create_bgin19_proof;
pub use verifier::{gen_vermsg, verify_bgin19_proof};

use crate::Vec;

/// The proof in FLPCP.
#[derive(Clone, Debug, Eq, PartialEq, CanonicalSerialize, CanonicalDeserialize)]
pub struct Proof<E: PairingEngine> {
    q_shares: Vec<Vec<E::Fr>>,
    // c_shares: Vec<E::Fr>,
}

/// The verification message in FLPCP.
#[derive(Clone, Debug, Eq, PartialEq, CanonicalSerialize, CanonicalDeserialize)]
pub struct VerMsg<E: PairingEngine> {
    q_r_share: E::Fr,
    f_r_shares: Vec<E::Fr>,
    b_share: E::Fr,
    // For debug
    q_r_shares: Vec<E::Fr>,
    q_coeffs: Vec<Vec<E::Fr>>,
}