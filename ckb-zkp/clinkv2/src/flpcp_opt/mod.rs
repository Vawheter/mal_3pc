use ark_ec::PairingEngine;
use ark_ff::Field;
use ark_serialize::*;
use crate::{
    kzg10::{Kzg10Proof, Kzg10Comm, ProveAssignment, ProveKey, KZG10},
};


pub mod prover;
pub mod verifier;

pub use prover::create_bgin19_proof;
pub use verifier::{gen_vermsg, verify_bgin19_proof};

use crate::{String, Vec};

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
    // p_r_value: E::Fr,
    // b: E::Fr,
}