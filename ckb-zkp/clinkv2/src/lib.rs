//! An implementation of the `CLINKv2`.
#![cfg_attr(not(feature = "std"), no_std)]
#![warn(unused, future_incompatible, nonstandard_style, rust_2018_idioms)]
#![allow(clippy::op_ref, clippy::suspicious_op_assign_impl)]
#![cfg_attr(not(use_asm), forbid(unsafe_code))]
#![cfg_attr(use_asm, feature(llvm_asm))]
#![cfg_attr(use_asm, deny(unsafe_code))]

#[cfg(not(feature = "std"))]
#[macro_use]
extern crate alloc;

#[cfg(not(feature = "std"))]
use alloc::{borrow::Cow, string::String, vec::Vec};

#[cfg(feature = "std")]
use std::vec::Vec;

pub mod batch_kzg10;
pub mod bgin19;
pub mod flpcp;
// pub mod flpcp_opt;
pub mod plonk_sumcheck;
// pub mod kzg10_distributed; // not work
// pub mod flpcp_opt2;
pub mod fliop;

pub mod plonk_fliop;

pub mod bbc19_fliop;

