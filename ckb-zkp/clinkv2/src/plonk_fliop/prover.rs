#![allow(non_snake_case)]

use std::vec;

use ark_ec::PairingEngine;
use ark_ff::{to_bytes, One, UniformRand, Zero, Field, PrimeField};
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain, Polynomial, UVPolynomial, domain};
use ark_std::rand::RngCore;
use merlin::Transcript;

use crate::batch_kzg10::batch_open;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::plonk_fliop::{Proof, Kzg10ComKey};
use ark_poly_commit::kzg10::KZG10;

// #[inline]
// fn dot_prod_subcircuit_P1<E: PairingEngine>(
//     eta: E::Fr,
//     inputs: Vec<E::Fr>,
// ) -> E::Fr {   
//     // eta * (delta_x * (delta_y + r_y1) + delta_y * r_x1 + rxy1 - rz1)
//     eta * (inputs[0] + (inputs[2] + inputs[3]) + inputs[2] * inputs[1] + inputs[4] - inputs[5])
// }

pub fn prove_dot_prod<E: PairingEngine, R: RngCore>(
    inputs: &Vec<Vec<E::Fr>>,
    rhos: &Vec<E::Fr>, // randomness for hiding f(X)
    k: usize, // compression parameter
    kzg10_ck: &Kzg10ComKey<'_, E>,
    sid: usize, // session id
    rng: &mut R,
) -> Proof<E> {
    // Inputs, T * (4n+3)
    let T: usize = inputs[0].len(); // number of dot products
    let L: usize = inputs.len(); // number of variables, (should be 4n+3, n for number of multiplications in one dot product)
    let n: usize = (L - 3) / 4;
    assert_eq!(4 * n + 3, L);
    assert_eq!(rhos.len(), L);

    let domain_k: GeneralEvaluationDomain<E::Fr> =
        EvaluationDomain::<E::Fr>::new(k).unwrap();
    // println!("domain_k: {:?}", domain_k);
    let domain_size = domain_k.size();
    let omega = domain_k.group_gen();
    // println!("domain_k size: {}", domain_size);

    let one = E::Fr::one();
    let zero = E::Fr::zero();

    // Hiding bound for committing
    let hiding_bound = Some(1);

    let S: usize = (T as f64 / k as f64).ceil() as usize;
    let mut s = S;
    let mut s0 = T;
    // println!("T: {:?}, k = {:?}", T, k);

    let mut transcript = Transcript::new(b"DZKP_DotProd");
    // Add public information
    transcript.append_message(b"k", &to_bytes!(k as u64).unwrap());
    transcript.append_message(b"sid", &to_bytes!(sid as u64).unwrap());
    
    // Inputs to circuits in iterations
    let mut f_evals_r: Vec<Vec<E::Fr>> = inputs.clone();
    let mut f_evals_r0 = vec![];

    // Randomness for linear combination
    let mut eta = zero;
    // Powers of eta
    let mut eta_powers: Vec<E::Fr> = vec![];
    let mut f_eta_evals_r: Vec<E::Fr> = vec![];
    let mut f_eta_evals_r0 =  vec![];
    let mut f_eta_polys: Vec<DensePolynomial<E::Fr>> = vec![];

    let mut G_eval_avg = zero;

    // Proof construction
    let mut poly_comms = vec![];
    let mut open_evals = vec![];
    let mut open_proofs = vec![];
    let mut G_evals = vec![];

    let mut r = zero;
    // For Fiat-Shamir
    let mut buf = [0u8; 31];
    let mut q_eval_r = zero;

    // println!("Begin Iteration...");
    let mut cnt = 0;
    loop {
        // println!("Iteration: {:?}", cnt);
        // cnt += 1;
        // println!("s0: {:?}", s0);
        // println!("s: {:?}", s);
        // println!("L: {:?}", L);

        // f(X), polynomials for input variables, s groups in total
        let mut f_polys: Vec<Vec<DensePolynomial<E::Fr>>> = vec![];
        let mut qt_eval_r = zero;

        let mut q_cos_evals = vec![zero; domain_size];
        let mut f_rxy_coeffs = vec![zero; domain_size];
        let mut f_rz_coeffs = vec![zero; domain_size];
        let mut f_deltaz_coeffs = vec![zero; domain_size];

        for l in 0..L {
            let mut f_l_polys: Vec<DensePolynomial<E::Fr>> = vec![];
            for i in 0..s {
                let cur = i * k;
                let end = if i == s - 1 {s0} else {cur + k}; // automatic padding in FFT
                
                let coeffs = domain_k.ifft(&f_evals_r[l][cur..end]);
                let mut poly = DensePolynomial::from_coefficients_vec(coeffs);
                f_l_polys.push(poly);
            }
            f_polys.push(f_l_polys);
        }
        // The last iteration, randomization needed
        // if s == 1 {
        //     for l in 0..L {
        //         let mut coeffs = vec![zero; domain_size];
        //         coeffs[0] = -rhos[l];
        //         coeffs[domain_size - 1] = rhos[l];
        //         f_polys[l][0] += &DensePolynomial::from_coefficients_vec(coeffs);
        //     }
        // }

        // First iteration, no need for computing and committing to U(x)
        if s == S {
            // Compute the quotient polynomial q(x)
            // IMPORTANT: number of coset points should be domain_size

            // let eta = E::Fr::rand(rng);
            transcript.challenge_bytes(b"eta", &mut buf);
            eta = <E::Fr>::from_random_bytes(&buf).unwrap();
            let mut eta_power = one;

            for _ in 0..s {
                eta_powers.push(eta_power);
                eta_power *= eta;
            }

            let f_eta_poly = domain_k.ifft(&eta_powers);
            let f_eta_cos_evals = domain_k.coset_fft(&f_eta_poly);

            for i in 0..s {
                for j in 0..n {
                    let cur = 4 * j;
                    let f_deltax_cos_evals = domain_k.coset_fft(&f_polys[cur][i]);
                    let f_rx_cos_evals = domain_k.coset_fft(&f_polys[cur + 1][i]);
                    let f_deltay_cos_evals = domain_k.coset_fft(&f_polys[cur + 2][i]);
                    let f_ry_cos_evals = domain_k.coset_fft(&f_polys[cur + 3][i]);
                    
                    for t in 0..k {
                        q_cos_evals[t] += f_eta_cos_evals[t] * (f_deltax_cos_evals[t] * (f_deltay_cos_evals[t] + f_ry_cos_evals[t]) + f_deltay_cos_evals[t] * f_rx_cos_evals[t]);
                    }
                    
                }
                for t in 0..k {
                    f_rxy_coeffs[t] += eta_powers[i] * f_polys[L - 3][i].coeffs[t];
                    f_rz_coeffs[t] += eta_powers[i] * f_polys[L - 2][i].coeffs[t];
                    f_deltaz_coeffs[t] += eta_powers[i] * f_polys[L - 1][i].coeffs[t];
                }
            }

            let f_rxy_cos_evals = domain_k.coset_fft(&f_rxy_coeffs);
            let f_rz_cos_evals = domain_k.coset_fft(&f_rz_coeffs);
            let f_deltaz_cos_evals = domain_k.coset_fft(&f_deltaz_coeffs);
            for t in 0..k {
                q_cos_evals[t] += f_rxy_cos_evals[t] + f_rz_cos_evals[t] - f_deltaz_cos_evals[t];
            }
            domain_k.divide_by_vanishing_poly_on_coset_in_place(&mut q_cos_evals);
            // Perform iFFT on coset
            let q_poly = DensePolynomial::from_coefficients_vec(domain_k.coset_ifft(&q_cos_evals));

            // Commit to q(x)
            let (q_comm, q_rand) = KZG10::<E, DensePolynomial<E::Fr>>::commit(&kzg10_ck, &q_poly, hiding_bound, Some(rng)).unwrap();

            transcript.append_message(b"q_comm", &to_bytes!(q_comm).unwrap());
            
            // Opening at random point r
            // let r = E::Fr::rand(rng);
            transcript.challenge_bytes(b"r", &mut buf);
            // Todo: Check r outside domain
            r = <E::Fr>::from_random_bytes(&buf).unwrap();

            q_eval_r = q_poly.evaluate(&r);
            // let q_r_proof = KZG10::<E, DensePolynomial<E::Fr>>::open(&kzg10_ck, &q_poly, r, &q_rand).unwrap();
            
            transcript.append_message(b"q_eval_r", &to_bytes!(q_eval_r).unwrap());
            // transcript.append_message(b"q_r_proof", &to_bytes!(q_r_proof).unwrap());

            poly_comms.push(q_comm);
            open_evals.push(q_eval_r);
            // open_proofs.push(q_r_proof);
            // G_evals.push()

            qt_eval_r = q_eval_r * domain_k.evaluate_vanishing_polynomial(r);
            G_evals.push(qt_eval_r);

            if s == 1 {
                break;
            }

            f_evals_r0 = inputs.clone();
            f_eta_evals_r0 = eta_powers.clone();
        }
        // the second round till the last round  
        else { 
            // let mut P_poly = DensePolynomial::<E::Fr>::zero();
            // for i in 0..s {
            //     for j in 0..n {
            //         let cur = 4 * j;
            //         P_poly += f_polys[cur][i] * (f_polys[cur + 2][i] + f_polys[cur + 3][i]) + f_polys[cur + 2][i] * f_polys[cur + 1][i];
            //     }
            //     P_poly += f_polys[L - 3][i] + f_polys[L - 2][i] - f_polys[L - 1][i];
            //     P_poly *= f_eta_polys[i];
            // }
            
            // P(X) = f_eta(X)\cdot G_i(f_1(X),...., f_{s_i}(X)) - q_{i-1}(r)t(r)/s_i
            // Compute U(X)
            let mut U_evals: Vec<E::Fr> = vec![zero; domain_size];
            let mut P_evals: Vec<E::Fr> = vec![zero; domain_size];
            // let mut U_j_eval = zero;
            // U(omega^0) = random
            let mut U_j_eval = E::Fr::rand(rng);
            U_evals[0] = U_j_eval;
            // // println!("{:?}", f_evals_r);
            // // println!("f_evals_r: {} * {}", f_evals_r.len(), f_evals_r[0].len());
            for j in 1..k {
                let mut P_j_eval = zero;
                for i in 0..s {
                    // let mut G_output = zero;
                    for j in 0..n {
                        let cur = 4 * j;
                        P_j_eval += f_evals_r[cur][i] * (f_evals_r[cur + 2][i] + f_evals_r[cur + 3][i]) + f_evals_r[cur + 2][i] * f_evals_r[cur + 1][i];
                    }
                    P_j_eval += f_evals_r[L - 3][i] + f_evals_r[L - 2][i] - f_evals_r[L - 1][i];
                    P_j_eval *= f_eta_evals_r[i];
                    P_j_eval -= G_eval_avg;
                }
            //         let mut G_output = zero;
            //         // let cur1 = if i == s - 1 {s0 - i * k + j} else {i * k + j};
            //         let mut cur1 = i * k + j;
            //         cur1 = if cur1 > s0 - 1 {s0 - 1} else {cur1};
            //         for t in 0..n {
            //             let cur2 = t * 4;
            //             // println!("[{}][{}]", cur2, cur1);
            //             G_output += f_evals_r[cur2][cur1] + (f_evals_r[cur2 + 2][cur1] + f_evals_r[cur2 + 3][cur1]) + f_evals_r[cur2 + 2][cur1] * f_evals_r[cur2 + 1][cur1];
            //         }
            //         G_output += f_evals_r[L - 3][cur1] + f_evals_r[L - 2][cur1] - f_evals_r[L - 1][cur1];
            //         P_j_eval += f_eta_evals_r[i] * G_output - G_eval_avg;
            //     }
                P_evals[j] = P_j_eval;
                U_j_eval += P_j_eval;
                U_evals[j] = U_j_eval;
            }
            // for j in 1..k {
            //     let cur = j * s;
            //     let P_j_eval = (0..s).map(|i|
            //         f_eta_evals_r[cur + i] * (f_evals_r[cur + i][0] + (f_evals_r[cur + i][2] + f_evals_r[cur + i][3]) + f_evals_r[cur + i][2] * f_evals_r[cur + i][1] + f_evals_r[cur + i][4] - f_evals_r[cur + i][5])
            //     ).sum::<E::Fr>();
            //     // U(omega^{i+1}) = U(omega^i) + P(omega^i)
            //     P_evals[j] = P_j_eval;
            //     U_j_eval += P_j_eval;
            //     U_evals[j] = U_j_eval;
            // }
            // let P_H_poly = DensePolynomial::from_coefficients_vec(domain_k.ifft(&P_evals));
            let mut U_poly = DensePolynomial::from_coefficients_vec(domain_k.ifft(&U_evals));
            
            if s == 1 {
                // let mu1 = E::Fr::rand(rng);
                // let mu2 = E::Fr::rand(rng);
                // let mut coeffs = vec![zero; domain_size + 1];
                // coeffs[0] = -mu2;
                // coeffs[1] = -mu1;
                // coeffs[domain_size - 1] = mu2;
                // coeffs[domain_size] = mu1;
                let mu = E::Fr::rand(rng);
                let mut coeffs = vec![zero; domain_size];
                coeffs[0] = -mu;
                coeffs[domain_size - 1] = mu;
                U_poly += &DensePolynomial::from_coefficients_vec(coeffs);
            }
            // Compute the quotient polynomial q(x)
            // let mut q_cos_evals = vec![zero; domain_size];
            let f_eta_poly = domain_k.ifft(&f_eta_evals_r);
            let f_eta_cos_evals = domain_k.coset_fft(&f_eta_poly);

            let mut U_omega_coeffs = U_poly.coeffs().to_vec();
            GeneralEvaluationDomain::<<E as PairingEngine>::Fr>::distribute_powers(&mut U_omega_coeffs, omega);
            let U_omega_poly = DensePolynomial::from_coefficients_vec(U_omega_coeffs);

            // // let P_cos_evals = domain_k.coset_fft(P_H_poly.coeffs());
            // let U_cos_evals = domain_k.coset_fft(U_poly.coeffs());
            // let U_omega_cos_evals = domain_k.coset_fft(U_omega_poly.coeffs());
            // for i in 0..k {
            //     q_cos_evals[i] = U_omega_cos_evals[i] - U_cos_evals[i] - P_cos_evals[i];
            // }

            for i in 0..s {
                for j in 0..n {
                    let cur = 4 * j;
                    let f_deltax_cos_evals = domain_k.coset_fft(&f_polys[cur][i]);
                    let f_rx_cos_evals = domain_k.coset_fft(&f_polys[cur + 1][i]);
                    let f_deltay_cos_evals = domain_k.coset_fft(&f_polys[cur + 2][i]);
                    let f_ry_cos_evals = domain_k.coset_fft(&f_polys[cur + 3][i]);
                    
                    for t in 0..k {
                        q_cos_evals[t] += f_eta_cos_evals[t] * (f_deltax_cos_evals[t] * (f_deltay_cos_evals[t] + f_ry_cos_evals[t]) + f_deltay_cos_evals[t] * f_rx_cos_evals[t]);
                    }
                    
                }
                for t in 0..k {
                    f_rxy_coeffs[t] += f_eta_evals_r[i] * f_polys[L - 3][i].coeffs[t];
                    f_rz_coeffs[t] += f_eta_evals_r[i] * f_polys[L - 2][i].coeffs[t];
                    f_deltaz_coeffs[t] += f_eta_evals_r[i] * f_polys[L - 1][i].coeffs[t];
                }
            }

            let f_rxy_cos_evals = domain_k.coset_fft(&f_rxy_coeffs);
            let f_rz_cos_evals = domain_k.coset_fft(&f_rz_coeffs);

            let U_cos_evals = domain_k.coset_fft(&U_poly);
            let U_omega_cos_evals = domain_k.coset_fft(&U_omega_poly);
            let f_deltaz_cos_evals = domain_k.coset_fft(&f_deltaz_coeffs);
            for t in 0..k {
                q_cos_evals[t] += U_omega_cos_evals[t] - U_cos_evals[t] - f_rxy_cos_evals[t] + f_rz_cos_evals[t] - f_deltaz_cos_evals[t];
            }
            
            domain_k.divide_by_vanishing_poly_on_coset_in_place(&mut q_cos_evals);
            // Perform iFFT on coset
            let q_poly = DensePolynomial::from_coefficients_vec(domain_k.coset_ifft(&q_cos_evals));

            // println!("q(X) degree: {:?}", q_poly.degree());
            // println!("U(X) degree: {:?}", U_poly.degree());
            // Commit to U(X), q(X)
            let (U_comm, U_rand) = KZG10::<E, DensePolynomial<E::Fr>>::commit(&kzg10_ck, &U_poly, hiding_bound, Some(rng)).unwrap();
            let (q_comm, q_rand) = KZG10::<E, DensePolynomial<E::Fr>>::commit(&kzg10_ck, &q_poly, hiding_bound, Some(rng)).unwrap();

            transcript.append_message(b"U_comm", &to_bytes!(U_comm).unwrap());
            transcript.append_message(b"q_comm", &to_bytes!(q_comm).unwrap());
            
            // Opening at random point r
            transcript.challenge_bytes(b"r", &mut buf);
            r = <E::Fr>::from_random_bytes(&buf).unwrap();
            let omega_r = omega * r;
            // println!("r: {}", r);
    
            let U_omega_eval_r = U_poly.evaluate(&omega_r);
            let U_eval_r = U_poly.evaluate(&r);
            q_eval_r = q_poly.evaluate(&r);

            transcript.append_message(b"U_omega_eval_r", &to_bytes!(U_omega_eval_r).unwrap());
            transcript.append_message(b"U_eval_r", &to_bytes!(U_eval_r).unwrap());
            transcript.append_message(b"q_eval_r", &to_bytes!(q_eval_r).unwrap());

            let U_omega_r_proof = KZG10::<E, DensePolynomial<E::Fr>>::open(&kzg10_ck, &U_poly, omega_r, &U_rand).unwrap();

            let mut polys = vec![];
            polys.push(U_poly);
            polys.push(q_poly);

            let mut poly_rands = vec![];
            poly_rands.push(U_rand);
            poly_rands.push(q_rand);

            transcript.challenge_bytes(b"open_challenge", &mut buf);
            let open_challenge = <E::Fr>::from_random_bytes(&buf).unwrap();
            
            let U_q_r_proof = batch_open::<E>(&kzg10_ck, &polys[..], r, open_challenge, &poly_rands[..].to_vec()).unwrap();
            
            transcript.append_message(b"U_omega_r_proof", &to_bytes!(U_omega_r_proof).unwrap());
            transcript.append_message(b"U_q_r_proof", &to_bytes!(U_q_r_proof).unwrap());

            poly_comms.push(U_comm);
            poly_comms.push(q_comm);
            
            open_evals.push(U_omega_eval_r);
            open_evals.push(U_eval_r);
            open_evals.push(q_eval_r);

            open_proofs.push(U_omega_r_proof);
            open_proofs.push(U_q_r_proof);

            qt_eval_r = q_eval_r * domain_k.evaluate_vanishing_polynomial(r);
            let G_eval = U_omega_eval_r - U_eval_r + G_eval_avg - qt_eval_r;
            G_evals.push(G_eval);

            if s == 1 {
                break;
            }

            f_evals_r0 = f_evals_r.clone();
            // println!("f_evals_r0: {:?}", f_evals_r0);
            f_eta_evals_r0 = f_eta_evals_r.clone();
        }

        // To the next iteration

        // Compute new f(r)
        f_evals_r = vec![];
        for l in 0..L {
            f_evals_r.push((0..s).map(|i| 
                f_polys[l][i].evaluate(&r)
            ).collect());
        }

        s0 = s;
        s = (s0 as f64 / k as f64).ceil() as usize;
        // println!("s: {:?}", s);

        // Compute q(r)t(r)/s
        let s_inv = <E::Fr>::from_repr((s as u64).into()).unwrap().inverse().unwrap();
        G_eval_avg = qt_eval_r * s_inv;

        // let lag_eval_rs = domain_k.evaluate_all_lagrange_coefficients(r);
        // for l in 0..L {
        //     let mut f_l_eval_r = vec![];
        //     for i in 0..s0 {
        //         let cur = i * k;
        //         let len = if i == s - 1 {s0 - cur} else {k}; // automatic padding in FFT
        //         let f_l_i_eval_r = (0..len).map(|j|
        //             f_evals_r0[l][cur + j] * lag_eval_rs[j]
        //         ).sum::<E::Fr>();
        //         f_l_eval_r.push(f_l_i_eval_r);
        //     }
        //     f_evals_r.push(f_l_eval_r);
        // }
        // println!("f_evals_r: {} * {}", f_evals_r.len(), f_evals_r[0].len());

        // for l in 0..L {
        //     let mut f_l_polys: Vec<DensePolynomial<E::Fr>> = vec![];
        //     for i in 0..s {
        //         let cur = i * k;
        //         let end = if i == s - 1 {s0} else {cur + k}; // automatic padding in FFT
                
        //         let coeffs = domain_k.ifft(&f_evals_r[l][cur..end]);
        //         let mut poly = DensePolynomial::from_coefficients_vec(coeffs);
        //         f_l_polys.push(poly);
        //     }
        //     f_polys.push(f_l_polys);
        // }

        // Compute f_eta(r)
        for i in 0..s {
            let cur = i * k;
            let end = if i == s - 1 {s0} else {cur + k}; // automatic padding in FFT
            
            let coeffs = domain_k.ifft(&f_eta_evals_r0[cur..end]);
            let mut poly = DensePolynomial::from_coefficients_vec(coeffs);
            f_eta_polys.push(poly);
        }

        f_eta_evals_r = vec![];
        for i in 0..s {
            f_eta_evals_r.push(f_eta_polys[i].evaluate(&r));
        }

        // println!("s: {:?}", s);
        // for i in 0..s {
        //     let cur = i * k;
        //     let len =  if i == s - 1 {s0 - cur} else {k};
        //     // println!("len: {:?}", len);
        //     f_eta_evals_r.push((0..len).map(|j|
        //         f_eta_evals_r0[cur + j] * lag_eval_rs[j]
        //         ).sum::<E::Fr>()
        //     );
        // }
        
    }

    // Out of loop
    Proof { 
        poly_comms, // U(X), q(X)
        open_evals, // U(omega_r), U(r), q(r)
        open_proofs, // U(omega_r), batch(U(r), q(r))
        G_evals,
    }
}
