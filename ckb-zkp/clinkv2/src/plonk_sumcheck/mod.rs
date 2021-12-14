#![allow(non_snake_case)]

use ark_ec::PairingEngine;
use ark_serialize::*;

pub mod prover;
pub mod verifier;

pub use prover::create_bgin19_proof;
pub use verifier::{gen_vermsg, verify_bgin19_proof};

pub type Kzg10Comm<E> = ark_poly_commit::kzg10::Commitment<E>;
pub type Kzg10Proof<E> = ark_poly_commit::kzg10::Proof<E>;
pub type Kzg10ComKey<'a, E> = ark_poly_commit::kzg10::Powers<'a, E>;

use crate::Vec;

/// The proof in PLONK Sumcheck.
#[derive(Clone, Debug, Eq, PartialEq, CanonicalSerialize, CanonicalDeserialize)]
pub struct Proof<E: PairingEngine> {
    poly_comms: Vec<Kzg10Comm<E>>,
    open_values: Vec<E::Fr>,
    open_proofs: Vec<Kzg10Proof<E>>,
    open_challenge: E::Fr,
}

/// The verification message in FLPCP.
#[derive(Clone, Debug, Eq, PartialEq, CanonicalSerialize, CanonicalDeserialize)]
pub struct VerMsg<E: PairingEngine> {
    F_z_shares: Vec<E::Fr>,
}