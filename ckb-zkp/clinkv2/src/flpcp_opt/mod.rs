use ark_ec::PairingEngine;
use ark_serialize::*;

pub mod prover;
pub mod verifier;

pub use prover::{create_bgin19_proof, Kzg10ComKey};
pub use verifier::{gen_vermsg, verify_bgin19_proof, Kzg10VerKey};

use crate::Vec;

pub type Kzg10Comm<E> = ark_poly_commit::kzg10::Commitment<E>;
pub type Kzg10Proof<E> = ark_poly_commit::kzg10::Proof<E>;

/// The proof in FLPCP.
#[derive(Clone, Debug, Eq, PartialEq, CanonicalSerialize, CanonicalDeserialize)]
pub struct Proof<E: PairingEngine> {
    q_comm: Kzg10Comm<E>,
    q_r_value: E::Fr,
    q_r_proof: Kzg10Proof<E>,
}

// #[derive(Clone, Debug, Eq, PartialEq, CanonicalSerialize, CanonicalDeserialize)]
// pub struct Proof2<E: PairingEngine> {
//     ws: Vec<E::Fr>,
// }

/// The verification message in FLPCP.
#[derive(Clone, Debug, Eq, PartialEq, CanonicalSerialize, CanonicalDeserialize)]
pub struct VerMsg<E: PairingEngine> {
    f_r_shares: Vec<E::Fr>,
}