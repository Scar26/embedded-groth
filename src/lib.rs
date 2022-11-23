#![no_std]

use pairing::Engine;

#[derive(Clone, Debug)]
pub struct Proof<E: Engine> {
    pub a: E::G1Affine,
    pub b: E::G2Affine,
    pub c: E::G1Affine,
}
#[derive(Clone)]
pub struct VerificationKey<E: Engine, const P: usize> {
    pub alpha_g1: E::G1Affine,

    pub beta_g1: E::G1Affine,
    pub beta_g2: E::G2Affine,

    pub gamma_g2: E::G2Affine,

    pub delta_g1: E::G1Affine,
    pub delta_g2: E::G2Affine,

    // LP_i = [(beta * A_i(tau) + alpha * B_i(tau) + C_i(tau))/gamma]*G_1
    // for all public inputs.
    pub ic: [E::G1Affine; P]
}

/// P: Number of public inputs
///
/// A: Number of aux inputs
/// 
/// Ensure M = P + A
#[derive(Clone)]
pub struct Parameters<E: Engine, const P: usize, const A: usize, const M: usize> {
    pub vk: VerificationKey<E, P>,

    // H query
    // h_i = (tau^i*Z_x(tau)/delta)*G1
    pub h: [E::G1Affine; A],

    // L query
    // l_i = (L_i(tau)/delta)*G1
    pub l: [E::G1Affine; A],

    pub a:  [E::G1Affine; M],
    
    pub b_g1: [E::G1Affine; M],
    pub b_g2: [E::G2Affine; M],
}