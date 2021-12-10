use ark_ec::PairingEngine;
use ark_serialize::*;

pub mod prover_points;
pub mod prover_fft;
pub mod verifier;

use crate::Vec;

pub use prover_fft::create_bgin19_proof as create_bgin19_proof_fft;
pub use prover_points::create_bgin19_proof as create_bgin19_proof_points;
pub use verifier::{gen_vermsg, verify_bgin19_proof};

/// The proof in FLPCP.
#[derive(Clone, Debug, Eq, PartialEq, CanonicalSerialize, CanonicalDeserialize)]
pub struct Proof<E: PairingEngine> {
    p_coeffs_shares: Vec<E::Fr>,
}

/// The verification message in FLPCP.
#[derive(Clone, Debug, Eq, PartialEq, CanonicalSerialize, CanonicalDeserialize)]
pub struct VerMsg<E: PairingEngine> {
    f_r_shares: Vec<Vec<E::Fr>>,
    p_r_share: E::Fr,
    b: E::Fr,
}