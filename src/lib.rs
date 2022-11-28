#![cfg_attr(not(feature = "std"), no_std)]

use pairing::Engine;

mod prover;

pub use prover::*;

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

#[cfg(feature = "std")]
mod assignments {
    use bellman::{ConstraintSystem, LinearCombination, SynthesisError, Variable, Index};
    use pairing::group::ff::PrimeField;

    pub struct ExtractAssignments<S: PrimeField> {
        input_assignment:  Vec<S>,
        aux_assignment: Vec<S>,
    }

    impl<S: PrimeField> ExtractAssignments<S> {
        pub fn get_input_count(&self) -> usize {
            self.input_assignment.len()
        }

        pub fn get_aux_count(&self) -> usize {
            self.aux_assignment.len()
        }

        pub fn get_assignments(self) -> (Vec<S>, Vec<S>) {
            (self.input_assignment, self.aux_assignment)
        }
    }

    impl<S: PrimeField> ConstraintSystem<S> for ExtractAssignments<S> {
        type Root = ExtractAssignments<S>;

        fn alloc<F, A, AR>(&mut self, _: A, f: F) -> Result<Variable, SynthesisError>
        where
            F: FnOnce() -> Result<S, SynthesisError>,
            A: FnOnce() -> AR,
            AR: Into<String>,
        {
            self.aux_assignment.push(f()?);
            Ok(Variable::new_unchecked(Index::Aux(self.get_aux_count() - 1)))
        }

        fn alloc_input<F, A, AR>(&mut self, _: A, f: F) -> Result<Variable, SynthesisError>
        where
            F: FnOnce() -> Result<S, SynthesisError>,
            A: FnOnce() -> AR,
            AR: Into<String>,
        {
            self.input_assignment.push(f()?);
            Ok(Variable::new_unchecked(Index::Aux(self.get_input_count() - 1)))
        }

        fn enforce<A, AR, LA, LB, LC>(&mut self, _: A, _: LA, _: LB, _: LC)
            where
                A: FnOnce() -> AR,
                AR: Into<String>,
                LA: FnOnce(LinearCombination<S>) -> LinearCombination<S>,
                LB: FnOnce(LinearCombination<S>) -> LinearCombination<S>,
                LC: FnOnce(LinearCombination<S>) -> LinearCombination<S>,
        {
            // Do nothing;
        }

        fn push_namespace<NR, N>(&mut self, _: N)
        where
            NR: Into<String>,
            N: FnOnce() -> NR,
        {
            // Do nothing;
        }
    
        fn pop_namespace(&mut self) {
            // Do nothing;
        }
    
        fn get_root(&mut self) -> &mut Self::Root {
            self
        }        
    }
}