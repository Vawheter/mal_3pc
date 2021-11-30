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

// in tatol 3G + 9F

/// The proof in FLPCP_OPT2.
#[derive(Clone, Debug, Eq, PartialEq, CanonicalSerialize, CanonicalDeserialize)]
pub struct Proof<E: PairingEngine> {
    q_comms: Vec<Kzg10Comm<E>>, // 2G
    q_r_values: Vec<E::Fr>, // 2F
    q_r_proof: Kzg10Proof<E>, // 1G
    opening_challenge: E::Fr, // 1F
}

/// The verification message in FLPCP_OPT2.
#[derive(Clone, Debug, Eq, PartialEq, CanonicalSerialize, CanonicalDeserialize)]
pub struct VerMsg<E: PairingEngine> {
    f_r_shares: Vec<E::Fr>, // 6F
}