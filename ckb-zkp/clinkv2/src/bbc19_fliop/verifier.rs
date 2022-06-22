// #![allow(non_snake_case)]

// use ark_ec::PairingEngine;
// use ark_ff::{to_bytes, One, Zero, Field, PrimeField};
// use ark_poly::polynomial::univariate::DensePolynomial;
// use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};

// use ark_poly_commit::kzg10::KZG10;
// pub type Kzg10VerKey<E> = ark_poly_commit::kzg10::VerifierKey<E>;
// use crate::batch_kzg10::batch_check;

// use merlin::Transcript;

// #[cfg(feature = "parallel")]
// use rayon::prelude::*;

// use crate::plonk_fliop::{Proof, VerMsg};

// pub fn verify_dot_prod_helper<E: PairingEngine>(
//     proof: &Proof<E>,
//     inputs: Vec<Vec<E::Fr>>,
//     rhos: &Vec<E::Fr>, // randomness for hiding f(X)
//     k: usize, // compression parameter
//     sid: usize, // session id
// ) -> Vec<Vec<E::Fr>> {
//    // Inputs shares
//    let T: usize = inputs[0].len(); // number of dot products
//    let L: usize = inputs.len(); // number of variables, (should be 4n+3, n for number of multiplications in one dot product)
//    let n: usize = (L - 3) / 4;
//    assert_eq!(4 * n + 3, L);
//    assert_eq!(rhos.len(), L);

//    let domain_k: GeneralEvaluationDomain<E::Fr> =
//        EvaluationDomain::<E::Fr>::new(k).unwrap();
//    // let domain_size = domain_k.size();
//    let omega = domain_k.group_gen();
//    // println!("domain_k size: {}", domain_size);

//    let one = E::Fr::one();
//    let zero = E::Fr::zero();

//    let mut s: usize = (T as f64 / k as f64).ceil() as usize;

//    // For Fiar-Shamir
//    let mut buf = [0u8; 31];

//    let mut transcript = Transcript::new(b"DZKP_DotProd");
//    // Add public information
//    transcript.append_message(b"k", &to_bytes!(k as u64).unwrap());
//    transcript.append_message(b"sid", &to_bytes!(sid as u64).unwrap());

//    transcript.challenge_bytes(b"eta", &mut buf);
//    let eta = <E::Fr>::from_random_bytes(&buf).unwrap();

//    // Proof Structure:
//    //     poly_comms, // q(x), { U(x), q(x) }
//    //     open_values, // q(r), { U(omega_r), U(r), q(z) }
//    //     open_proofs, // q(r), { pi_U(omega_r), pi_U_q(r) }
//    //     open_challenge 
   
//    let mut r = zero;

// //    let s_inv = <E::Fr>::from_repr((s as u64).into()).unwrap().inverse().unwrap();
// //    let mut G_eval_avg = proof.open_evals[0] * domain_k.evaluate_vanishing_polynomial(r) * s_inv;

//    let mut s0  = T;
//    // Number of recursive iterations
//    let t = (proof.open_proofs.len() - 1) / 2 + 1;

//    let mut eta_powers = vec![];
//    let mut eta_power = one;
   
//    let mut f_eta_evals_r: Vec<E::Fr>  = vec![];
//    let mut f_evals_r: Vec<Vec<E::Fr>> = vec![];

//    // let mut G_eval_avg = zero;
// //    let mut qt_r_eval = zero;

//    let mut f_eta_evals_r0: Vec<E::Fr>  = vec![];
//    let mut f_evals_r0: Vec<Vec<E::Fr>> = vec![];

//    // Check opening proofs
//    for j in 0..t {

//        if j == 0 {
//            transcript.append_message(b"q_comm", &to_bytes!(&proof.poly_comms[0]).unwrap());

//            transcript.challenge_bytes(b"r", &mut buf);
//            r = <E::Fr>::from_random_bytes(&buf).unwrap();
//            // Check opening proof in the first iteration
//         //    assert!(KZG10::<E, DensePolynomial<E::Fr>>::check(&kzg10_vk, &proof.poly_comms[0], r, proof.open_evals[0], &proof.open_proofs[0]).unwrap());
       
//            for _ in 0..s {
//                eta_power *= eta;
//                eta_powers.push(eta_power);
//            }
   
//         //    qt_r_eval = proof.open_evals[0] * domain_k.evaluate_vanishing_polynomial(r);
//         //    assert_eq!(proof.G_evals[0], qt_r_eval);

//            transcript.append_message(b"q_eval_r", &to_bytes!(&proof.open_evals[0]).unwrap());
//            transcript.append_message(b"q_r_proof", &to_bytes!(&proof.open_proofs[0]).unwrap());


//            f_evals_r0 = inputs.clone();
//            f_eta_evals_r0 = eta_powers.clone();
//        }
//        else {
//            let cur1 = 2 * j - 1;
//            let cur2 = 3 * j - 2;

//            transcript.append_message(b"U_comm", &to_bytes!(&proof.poly_comms[cur1]).unwrap());
//            transcript.append_message(b"q_comm", &to_bytes!(&proof.poly_comms[cur1 + 1]).unwrap());

//            transcript.challenge_bytes(b"r", &mut buf);
//            r = <E::Fr>::from_random_bytes(&buf).unwrap();
//            let omega_r = omega * r;
//         //    println!("r: {}", r);

//            transcript.append_message(b"U_omega_eval_r", &to_bytes!(proof.open_evals[cur2]).unwrap());
//            transcript.append_message(b"U_eval_r", &to_bytes!(proof.open_evals[cur2 + 1]).unwrap());
//            transcript.append_message(b"q_eval_r", &to_bytes!(proof.open_evals[cur2 + 2]).unwrap());

//            transcript.challenge_bytes(b"open_challenge", &mut buf);
//            let open_challenge = <E::Fr>::from_random_bytes(&buf).unwrap();

//         //    // Check U(omega_r)
//         //    assert!(KZG10::<E, DensePolynomial<E::Fr>>::check(&kzg10_vk, &proof.poly_comms[cur1], omega_r, proof.open_evals[cur2], &proof.open_proofs[cur1]).unwrap());
//         //    // Check U(r), q(r)
//         //    assert!(batch_check(&kzg10_vk, &proof.poly_comms[cur1..cur1 + 2], r, &proof.open_evals[cur2 + 1..cur2 + 3], &proof.open_proofs[cur1 + 1], open_challenge).unwrap());    
       

//            transcript.append_message(b"U_omega_r_proof", &to_bytes!(&proof.open_proofs[cur1]).unwrap());
//            transcript.append_message(b"U_q_r_proof", &to_bytes!(&proof.open_proofs[cur1 + 1]).unwrap());
           
//         //    qt_r_eval = proof.open_evals[cur2 + 2] * domain_k.evaluate_vanishing_polynomial(r);
//         //    let left = proof.open_evals[cur2] - proof.open_evals[cur2 + 1] - (proof.G_evals[j] - G_eval_avg);
//         //    assert_eq!(left, qt_r_eval);

//            f_evals_r0 = f_evals_r.clone();
//            f_eta_evals_r0 = f_eta_evals_r.clone();
//        }
       
//        let lag_r_evals = domain_k.evaluate_all_lagrange_coefficients(r);
//        f_evals_r = vec![];
       
//        for l in 0..L {
//            let mut f_l_eval_r = vec![];
//            for i in 0..s {
//                let cur = i * k;
//                let len: usize = if i == s - 1 {s0 - cur} else {k}; // automatic padding in FFT
//                let f_l_i_eval_r = (0..len).map(|j|
//                    f_evals_r0[l][cur + j] * lag_r_evals[j]
//                ).sum::<E::Fr>();
//                f_l_eval_r.push(f_l_i_eval_r);
//            }
//            f_evals_r.push(f_l_eval_r);
//        }

//        // f_eta_evals_r = vec![];
//        // for i in 0..s {
//        //     let len =  if i == s {s0 - s * (k - 1)} else {k};
//        //     f_eta_evals_r.push((0..len).map(|q|
//        //         f_eta_evals_r0[q * s + i] * lag_r_evals[q]
//        //         ).sum::<E::Fr>()
//        //     );
//        // }
//        s0 = s;
//        s = (s0 as f64 / k as f64).ceil() as usize;

//        f_eta_evals_r = vec![];
//        // println!("f_eta_evals_r0.len():{}", f_eta_evals_r0.len());
//        // println!("s: {:?}", s);
//        for i in 0..s {
//            let cur = i * k;
//            let len =  if i == s - 1 {s0 - cur} else {k};
//            // println!("len: {:?}", len);
//            f_eta_evals_r.push((0..len).map(|j|
//                f_eta_evals_r0[cur + j] * lag_r_evals[j]
//                ).sum::<E::Fr>()
//            );
//        }

//        f_evals_r0 = f_evals_r.clone();
//        f_eta_evals_r0 = f_eta_evals_r.clone();
//    }

//     f_evals_r
// }

// pub fn verify_dot_prod<E: PairingEngine>(
//     f_evals_r_helper: Vec<Vec<E::Fr>>,
//     proof: &Proof<E>,
//     k: usize, // compression parameter
//     inputs: &Vec<Vec<E::Fr>>,
//     rhos: &Vec<E::Fr>, // randomness for hiding f(X)
//     kzg10_vk: &Kzg10VerKey<E>,
//     sid: usize, // session id
// ) -> bool {
//     // Inputs shares
//     let T: usize = inputs[0].len(); // number of dot products
//     let L: usize = inputs.len(); // number of variables, (should be 4n+3, n for number of multiplications in one dot product)
//     let n: usize = (L - 3) / 4;
//     assert_eq!(4 * n + 3, L);
//     assert_eq!(rhos.len(), L);

//     let domain_k: GeneralEvaluationDomain<E::Fr> =
//         EvaluationDomain::<E::Fr>::new(k).unwrap();
//     // let domain_size = domain_k.size();
//     let omega = domain_k.group_gen();
//     // println!("domain_k size: {}", domain_size);

//     let one = E::Fr::one();
//     let zero = E::Fr::zero();

//     let mut s: usize = (T as f64 / k as f64).ceil() as usize;

//     // For Fiar-Shamir
//     let mut buf = [0u8; 31];

//     let mut transcript = Transcript::new(b"DZKP_DotProd");
//     // Add public information
//     transcript.append_message(b"k", &to_bytes!(k as u64).unwrap());
//     transcript.append_message(b"sid", &to_bytes!(sid as u64).unwrap());

//     transcript.challenge_bytes(b"eta", &mut buf);
//     let eta = <E::Fr>::from_random_bytes(&buf).unwrap();

//     // Proof Structure:
//     //     poly_comms, // q(x), { U(x), q(x) }
//     //     open_values, // q(r), { U(omega_r), U(r), q(z) }
//     //     open_proofs, // q(r), { pi_U(omega_r), pi_U_q(r) }
//     //     open_challenge 
    
//     let mut r = zero;

//     let mut s0  = T;
//     // Number of recursive iterations
//     let t = (proof.open_proofs.len() - 1) / 2 + 1;

//     let mut eta_powers = vec![];
//     let mut eta_power = one;
    
//     let mut f_eta_evals_r: Vec<E::Fr>  = vec![];
//     let mut f_evals_r: Vec<Vec<E::Fr>> = vec![];

//     let mut G_eval_avg = zero;
//     let mut qt_r_eval = zero;

//     let mut f_eta_evals_r0: Vec<E::Fr>  = vec![];
//     let mut f_evals_r0: Vec<Vec<E::Fr>> = vec![];

//     // Check opening proofs
//     for j in 0..t {

//         if j == 0 {
//             transcript.append_message(b"q_comm", &to_bytes!(&proof.poly_comms[0]).unwrap());

//             transcript.challenge_bytes(b"r", &mut buf);
//             r = <E::Fr>::from_random_bytes(&buf).unwrap();
//             // Check opening proof in the first iteration
//             assert!(KZG10::<E, DensePolynomial<E::Fr>>::check(&kzg10_vk, &proof.poly_comms[0], r, proof.open_evals[0], &proof.open_proofs[0]).unwrap());
        
//             for _ in 0..s {
//                 eta_power *= eta;
//                 eta_powers.push(eta_power);
//             }
    
//             qt_r_eval = proof.open_evals[0] * domain_k.evaluate_vanishing_polynomial(r);
//             assert_eq!(proof.G_evals[0], qt_r_eval);

//             transcript.append_message(b"q_eval_r", &to_bytes!(&proof.open_evals[0]).unwrap());
//             transcript.append_message(b"q_r_proof", &to_bytes!(&proof.open_proofs[0]).unwrap());

//             f_evals_r0 = inputs.clone();
//             f_eta_evals_r0 = eta_powers.clone();
//         }
//         else {
//             let cur1 = 2 * j - 1;
//             let cur2 = 3 * j - 2;

//             transcript.append_message(b"U_comm", &to_bytes!(&proof.poly_comms[cur1]).unwrap());
//             transcript.append_message(b"q_comm", &to_bytes!(&proof.poly_comms[cur1 + 1]).unwrap());

//             transcript.challenge_bytes(b"r", &mut buf);
//             r = <E::Fr>::from_random_bytes(&buf).unwrap();
//             let omega_r = omega * r;
//             // println!("r: {}", r);

//             transcript.append_message(b"U_omega_eval_r", &to_bytes!(proof.open_evals[cur2]).unwrap());
//             transcript.append_message(b"U_eval_r", &to_bytes!(proof.open_evals[cur2 + 1]).unwrap());
//             transcript.append_message(b"q_eval_r", &to_bytes!(proof.open_evals[cur2 + 2]).unwrap());

//             transcript.challenge_bytes(b"open_challenge", &mut buf);
//             let open_challenge = <E::Fr>::from_random_bytes(&buf).unwrap();

//             // Check U(omega_r)
//             assert!(KZG10::<E, DensePolynomial<E::Fr>>::check(&kzg10_vk, &proof.poly_comms[cur1], omega_r, proof.open_evals[cur2], &proof.open_proofs[cur1]).unwrap());
//             // Check U(r), q(r)
//             assert!(batch_check(&kzg10_vk, &proof.poly_comms[cur1..cur1 + 2], r, &proof.open_evals[cur2 + 1..cur2 + 3], &proof.open_proofs[cur1 + 1], open_challenge).unwrap());    
        

//             transcript.append_message(b"U_omega_r_proof", &to_bytes!(&proof.open_proofs[cur1]).unwrap());
//             transcript.append_message(b"U_q_r_proof", &to_bytes!(&proof.open_proofs[cur1 + 1]).unwrap());
            
//             qt_r_eval = proof.open_evals[cur2 + 2] * domain_k.evaluate_vanishing_polynomial(r);
//             let left = proof.open_evals[cur2] - proof.open_evals[cur2 + 1] - (proof.G_evals[j] - G_eval_avg);
//             assert_eq!(left, qt_r_eval);

//             f_evals_r0 = f_evals_r.clone();
//             f_eta_evals_r0 = f_eta_evals_r.clone();
//         }
        
//         let lag_r_evals = domain_k.evaluate_all_lagrange_coefficients(r);
//         f_evals_r = vec![];
        
//         for l in 0..L {
//             let mut f_l_eval_r = vec![];
//             for i in 0..s {
//                 let cur = i * k;
//                 let len: usize = if i == s - 1 {s0 - cur} else {k}; // automatic padding in FFT
//                 let f_l_i_eval_r = (0..len).map(|j|
//                     f_evals_r0[l][cur + j] * lag_r_evals[j]
//                 ).sum::<E::Fr>();
//                 f_l_eval_r.push(f_l_i_eval_r);
//             }
//             f_evals_r.push(f_l_eval_r);
//         }

//         // f_eta_evals_r = vec![];
//         // for i in 0..s {
//         //     let len =  if i == s {s0 - s * (k - 1)} else {k};
//         //     f_eta_evals_r.push((0..len).map(|q|
//         //         f_eta_evals_r0[q * s + i] * lag_r_evals[q]
//         //         ).sum::<E::Fr>()
//         //     );
//         // }
//         s0 = s;
//         s = (s0 as f64 / k as f64).ceil() as usize;

//         let s_inv = <E::Fr>::from_repr((s as u64).into()).unwrap().inverse().unwrap();
//         G_eval_avg = qt_r_eval * s_inv;


//         f_eta_evals_r = vec![];
//         // println!("f_eta_evals_r0.len():{}", f_eta_evals_r0.len());
//         // println!("s: {:?}", s);
//         for i in 0..s {
//             let cur = i * k;
//             let len =  if i == s - 1 {s0 - cur} else {k};
//             // println!("len: {:?}", len);
//             f_eta_evals_r.push((0..len).map(|j|
//                 f_eta_evals_r0[cur + j] * lag_r_evals[j]
//                 ).sum::<E::Fr>()
//             );
//         }

//         f_evals_r0 = f_evals_r.clone();
//         f_eta_evals_r0 = f_eta_evals_r.clone();


//     }

//     let mut G_eval_should_be = zero;
//     for j in 0..n {
//         let cur = 4 * j;
//         G_eval_should_be += (f_evals_r[cur][0] + f_evals_r_helper[cur][0]) * (f_evals_r[cur + 2][0] + f_evals_r_helper[cur + 2][0] + f_evals_r[cur + 3][0] + f_evals_r_helper[cur + 3][0]) + (f_evals_r[cur + 2][0] + f_evals_r_helper[cur + 2][0]) * (f_evals_r[cur + 1][0] + f_evals_r_helper[cur + 1][0]);
//     }
//     G_eval_should_be += f_evals_r[L - 3][0] + f_evals_r[L - 2][0] - f_evals_r[L - 1][0];
//     // G_eval_should_be *= f_eta_evals_r[0];
//     assert_eq!(G_eval_should_be, proof.G_evals[t - 1]);

//     true
// }

