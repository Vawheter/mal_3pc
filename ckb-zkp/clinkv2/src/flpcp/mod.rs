use ark_ec::PairingEngine;
use ark_ff::Field;
use ark_serialize::*;

pub mod prover_points;
pub mod prover_fft;
pub mod verifier;

// pub use verifier::{gen_vermsg, verify_bgin19_proof};

use crate::{String, Vec};

pub use prover_fft::create_bgin19_proof as create_bgin19_proof_fft;
pub use prover_points::create_bgin19_proof as create_bgin19_proof_points;

/// The proof in FLPCP.
#[derive(Clone, Debug, Eq, PartialEq, CanonicalSerialize, CanonicalDeserialize)]
pub struct Proof<E: PairingEngine> {
    // ws: Vec<E::Fr>,
    p_coeffs_shares: Vec<E::Fr>,
}

/// The verification message in FLPCP.
#[derive(Clone, Debug, Eq, PartialEq, CanonicalSerialize, CanonicalDeserialize)]
pub struct VerMsg<E: PairingEngine> {
    f_r_shares: Vec<E::Fr>,
    p_r_share: E::Fr,
    b: E::Fr,
}