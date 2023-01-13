#![cfg_attr(not(any(test, feature = "std")), no_std)]

#[cfg(not(any(test, feature = "std")))]
extern crate alloc;

#[cfg(test)]
#[macro_use]
extern crate std;

use pairing::Engine;

// pub mod prover;
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

#[cfg(any(test, feature = "std"))]
pub mod assignments {
    use bellman::{ConstraintSystem, LinearCombination, SynthesisError, Variable, Index, Circuit};
    use pairing::group::ff::{ Field, PrimeField };
    use super::*;

    #[derive(Default, Debug)]
    pub struct ExtractAssignments<S: PrimeField> {
        input_assignment:  Vec<S>,
        num_inputs: usize,
        aux_assignment: Vec<S>,
        num_aux: usize,
        a_constraints: Vec<Index>,
        b_constraints: Vec<Index>,
    }

    impl<S: PrimeField> ExtractAssignments<S> {
        pub fn get_assignments(&self) -> (Vec<S>, Vec<S>) {
            (self.input_assignment.clone(), self.aux_assignment.clone())
        }

        pub fn get_constraints(&self) -> (Vec<bool>, Vec<bool>) {
            fn eval(input: Vec<Index>, p: usize, a: usize) -> Vec<bool> {
                let mut out = vec![false; p+a];
                
                for idx in input {
                    match idx {
                        Index::Input(i) => {
                            out[i] = true;
                        }

                        Index::Aux(i) => {
                            out[p+i] = true;
                        }
                    }
                }
                out
            }  
            
            (eval(self.a_constraints.clone(), self.num_inputs, self.num_aux), eval(self.b_constraints.clone(),self.num_inputs, self.num_aux))
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
            self.num_aux += 1;
            Ok(Variable::new_unchecked(Index::Aux(self.num_aux - 1)))
        }

        fn alloc_input<F, A, AR>(&mut self, _: A, f: F) -> Result<Variable, SynthesisError>
        where
            F: FnOnce() -> Result<S, SynthesisError>,
            A: FnOnce() -> AR,
            AR: Into<String>,
        {
            self.input_assignment.push(f()?);
            self.num_inputs += 1;
            Ok(Variable::new_unchecked(Index::Aux(self.num_inputs - 1)))
        }

        fn enforce<A, AR, LA, LB, LC>(&mut self, _: A, a: LA, b: LB, _: LC)
            where
                A: FnOnce() -> AR,
                AR: Into<String>,
                LA: FnOnce(LinearCombination<S>) -> LinearCombination<S>,
                LB: FnOnce(LinearCombination<S>) -> LinearCombination<S>,
                LC: FnOnce(LinearCombination<S>) -> LinearCombination<S>,
        {
            fn eval<S: PrimeField>(
                lc: LinearCombination<S>,
                output: &mut Vec<Index>
            ) {
                for (var, _) in lc.as_ref() {
                    output.push(var.get_unchecked());
                }
            }
            
            eval(a(LinearCombination::zero()), &mut self.a_constraints);
            eval(b(LinearCombination::zero()), &mut self.b_constraints);
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

    pub fn extract_assignments<C, E>(circuit: C) -> Result<ExtractAssignments<E::Fr>, SynthesisError>
    where
        E: Engine,
        C: Circuit<E::Fr>
    {
        let mut cs = ExtractAssignments::<E::Fr>::default();
        cs.alloc_input(|| "one", || Ok(E::Fr::one()))?;
        circuit.synthesize(&mut cs)?;
        for i in 0..cs.num_inputs {
            cs.enforce(|| "", |lc| lc + Variable::new_unchecked(Index::Input(i)), |lc| lc, |lc| lc);
        }
        Ok(cs)
    }
}


#[cfg(test)]
mod tests {
    use super::assignments;
    use std::mem;
    use bls12_381::{ Scalar };

    #[test]
    fn test() {
        let a = assignments::ExtractAssignments::<Scalar>::default();
        println!("{}", mem::size_of_val(&a));
    }
}