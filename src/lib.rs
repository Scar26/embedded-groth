#![cfg_attr(not(feature = "std"), no_std)]

use pairing::Engine;

#[derive(Clone, Debug)]
pub struct Proof<E: Engine> {
    pub a: E::G1Affine,
    pub b: E::G2Affine,
    pub c: E::G1Affine,
}
#[derive(Clone)]
pub struct VerificationKey<'a, E: Engine> {
    pub alpha_g1: E::G1Affine,

    pub beta_g1: E::G1Affine,
    pub beta_g2: E::G2Affine,

    pub gamma_g2: E::G2Affine,

    pub delta_g1: E::G1Affine,
    pub delta_g2: E::G2Affine,

    // LP_i = [(beta * A_i(tau) + alpha * B_i(tau) + C_i(tau))/gamma]*G_1
    // for all public inputs.
    pub ic: &'a [E::G1Affine]
}

#[derive(Clone)]
pub struct Parameters<'a, 'b, E: Engine> {
    pub vk: VerificationKey<'b, E>,

    // H query
    // h_i = (tau^i*Z_x(tau)/delta)*G1
    pub h: &'a [E::G1Affine],

    // L query
    // l_i = (L_i(tau)/delta)*G1
    pub l: &'a [E::G1Affine],

    pub a:  &'a [E::G1Affine],
    
    pub b_g1: &'a [E::G1Affine],
    pub b_g2: &'a [E::G2Affine],
}

pub struct Test<const P: usize, const A: usize, const M: usize> {
    pub a:  [usize; A],
    
    pub b_g1: [usize; P],
    pub b_g2: [usize; M],
}

#[cfg(feature = "std")]
pub mod assignments {
    use bellman::{ConstraintSystem, LinearCombination, SynthesisError, Variable, Index};
    use pairing::group::ff::PrimeField;

    #[derive(Default)]
    pub struct ExtractAssignments<S: PrimeField> {
        input_assignment:  Vec<S>,
        aux_assignment: Vec<S>,
        a_assignments: Vec<usize>,
        b_assignments: Vec<usize>
    }

    impl<S: PrimeField> ExtractAssignments<S> {
        pub fn get_input_count(&self) -> usize {
            self.input_assignment.len()
        }

        pub fn get_aux_count(&self) -> usize {
            self.aux_assignment.len()
        }

        pub fn to_bytes(&self) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
            let input_assignments = self.input_assignment.iter()
                .map(|s| s.to_repr().as_ref().to_owned())
                .collect();
            
            let aux_assignments = self.aux_assignment.iter()
                .map(|s| s.to_repr().as_ref().to_owned())
                .collect();
            
            (input_assignments, aux_assignments)
        }

        pub fn from_bytes(_: Vec<Vec<u8>>, _: Vec<Vec<u8>>) -> Self {
            // [`PrimeField::from_repr`] does not allow conversion from &[u8] slices
            // Deserialization is Engine dependent
            unimplemented!()
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

        fn enforce<A, AR, LA, LB, LC>(&mut self, _: A, a: LA, b: LB, _: LC)
            where
                A: FnOnce() -> AR,
                AR: Into<String>,
                LA: FnOnce(LinearCombination<S>) -> LinearCombination<S>,
                LB: FnOnce(LinearCombination<S>) -> LinearCombination<S>,
                LC: FnOnce(LinearCombination<S>) -> LinearCombination<S>,
        {
            let a = a(LinearCombination::zero()).as_ref();
            for (var, _) in a {
                match var {
                    Variable(Index::Input(p)) => {
                        unreachable!()
                    }

                    Variable(Index::Aux(u)) => {

                    }
                }
            }
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
