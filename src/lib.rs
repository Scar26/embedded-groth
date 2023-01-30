#![cfg_attr(not(any(test, feature = "std")), no_std)]

#[cfg(not(any(test, feature = "std")))]
extern crate alloc;
#[cfg(not(any(test, feature = "std")))]
use alloc::vec::Vec;

use ff::PrimeField;
use pairing::Engine;

pub mod prover;
pub mod verifier;
mod poly;

pub enum VerificationError {
    InvalidProof,
    InvalidVerifyingKey,
}


#[derive(Clone, Debug)]
pub struct Proof<E: Engine> {
    pub a: E::G1Affine,
    pub b: E::G2Affine,
    pub c: E::G1Affine,
}
#[derive(Clone)]
pub struct VerificationKey<E: Engine> {
    pub alpha_g1: E::G1Affine,

    pub beta_g1: E::G1Affine,
    pub beta_g2: E::G2Affine,

    pub gamma_g2: E::G2Affine,

    pub delta_g1: E::G1Affine,
    pub delta_g2: E::G2Affine,

    // LP_i = [(beta * A_i(tau) + alpha * B_i(tau) + C_i(tau))/gamma]*G_1
    // for all public inputs.
    pub ic: Vec<E::G1Affine>
}

#[derive(Clone)]
pub struct Parameters<E: Engine> {
    pub vk: VerificationKey<E>,

    // H query
    // h_i = (tau^i*Z_x(tau)/delta)*G1
    pub h: Vec<E::G1Affine>,

    // L query
    // l_i = (L_i(tau)/delta)*G1
    pub l: Vec<E::G1Affine>,

    pub a:  Vec<E::G1Affine>,
    
    pub b_g1: Vec<E::G1Affine>,
    pub b_g2: Vec<E::G2Affine>,
}

/// The first index is not really needed since bellman doesn't allow
/// unconstrained variables but it's included just for cleanliness and 
/// potential multi threaded speedups
#[derive(Default, Debug, Clone)]
pub struct QAP<S: PrimeField> {
    pub a: Vec<(usize, Vec<(S, usize)>)>,
    pub b: Vec<(usize, Vec<(S, usize)>)>,
    pub c: Vec<(usize, Vec<(S, usize)>)>,
}

#[cfg(any(test, feature = "std"))]
pub mod assignments {
    use std::collections::HashMap;
    use bellman::{ConstraintSystem, LinearCombination, SynthesisError, Variable, Index, Circuit};
    use pairing::group::ff::{ Field, PrimeField };
    use super::*;
    #[derive(Default, Debug)]
    pub struct AnalyzeCircuit<S: PrimeField> {
        input_assignment:  Vec<S>,
        num_inputs: usize,
        aux_assignment: Vec<S>,
        num_aux: usize,
        num_constraints: usize,
        extract_assignments: bool,
        at: Vec<(Index, S, usize)>,
        bt: Vec<(Index, S, usize)>,
        ct: Vec<(Index, S, usize)>,
    }

    impl<S: PrimeField> AnalyzeCircuit<S> {
        pub fn get_num_states(&self) -> (usize, usize) {
            (self.num_inputs, self.num_aux)
        }

        pub fn get_assignments(&self) -> (Vec<S>, Vec<S>) {
            (self.input_assignment.clone(), self.aux_assignment.clone())
        }

        // Only call after synthesize
        pub fn qap(self) -> QAP<S> {
            fn collect<S: PrimeField>(v: Vec<(Index, S, usize)>, p: usize) -> Vec<(usize, Vec<(S, usize)>)> {
                let mut map: HashMap<usize, Vec<(S, usize)>> = HashMap::new();
                for (var, coeff, constraint) in v.into_iter() {
                    let i = match var {
                        Index::Input(i) => i,
                        Index::Aux(i) => p + i,
                    };

                    match map.get_mut(&i) {
                        Some(constraints) => {
                            constraints.push((coeff, constraint));
                        },

                        None => { 
                            map.insert(i, vec![(coeff, constraint)]);
                        },
                    }
                } 
                map.into_iter().collect()
            }

            QAP {
                a: collect(self.at, self.num_inputs),
                b: collect(self.bt, self.num_inputs),
                c: collect(self.ct, self.num_inputs),
            }
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

    impl<S: PrimeField> ConstraintSystem<S> for AnalyzeCircuit<S> {
        type Root = AnalyzeCircuit<S>;

        fn alloc<F, A, AR>(&mut self, _: A, f: F) -> Result<Variable, SynthesisError>
        where
            F: FnOnce() -> Result<S, SynthesisError>,
            A: FnOnce() -> AR,
            AR: Into<String>,
        {
            if self.extract_assignments {
                self.aux_assignment.push(f()?);
            }
            self.num_aux += 1;
            Ok(Variable::new_unchecked(Index::Aux(self.num_aux - 1)))
        }

        fn alloc_input<F, A, AR>(&mut self, _: A, f: F) -> Result<Variable, SynthesisError>
        where
            F: FnOnce() -> Result<S, SynthesisError>,
            A: FnOnce() -> AR,
            AR: Into<String>,
        {
            if self.extract_assignments {
                self.input_assignment.push(f()?);
            }
            self.num_inputs += 1;
            Ok(Variable::new_unchecked(Index::Aux(self.num_inputs - 1)))
        }

        fn enforce<A, AR, LA, LB, LC>(&mut self, _: A, a: LA, b: LB, c: LC)
            where
                A: FnOnce() -> AR,
                AR: Into<String>,
                LA: FnOnce(LinearCombination<S>) -> LinearCombination<S>,
                LB: FnOnce(LinearCombination<S>) -> LinearCombination<S>,
                LC: FnOnce(LinearCombination<S>) -> LinearCombination<S>,
        {
            fn eval<S: PrimeField>(
                lc: LinearCombination<S>,
                output: &mut Vec<(Index, S, usize)>,
                current_constraint: usize,
            ) {
                for (var, c) in lc.as_ref() {
                    output.push((var.get_unchecked(), c.clone(), current_constraint))
                }
            }
            
            eval(a(LinearCombination::zero()), &mut self.at, self.num_constraints);
            eval(b(LinearCombination::zero()), &mut self.bt, self.num_constraints);
            self.num_constraints += 1;
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

    pub fn extract_assignments<C, E>(circuit: C) -> Result<AnalyzeCircuit<E::Fr>, SynthesisError>
    where
        E: Engine,
        C: Circuit<E::Fr>
    {
        let mut cs = AnalyzeCircuit::<E::Fr>{
            extract_assignments: true,
            ..Default::default()
        };

        cs.alloc_input(|| "one", || Ok(E::Fr::one()))?;
        circuit.synthesize(&mut cs)?;
        for i in 0..cs.num_inputs {
            cs.enforce(|| "", |lc| lc + Variable::new_unchecked(Index::Input(i)), |lc| lc, |lc| lc);
        }
        Ok(cs)
    }

    pub fn extract_circuit<C, S>(circuit: C) -> Result<QAP<S>, SynthesisError>
    where
        S: PrimeField,
        C: Circuit<S>
    {
        let mut cs = AnalyzeCircuit::<S>{
            extract_assignments: false,
            ..Default::default()
        };
        
        cs.alloc_input(|| "one", || Ok(S::one()))?;
        circuit.synthesize(&mut cs)?;
        for i in 0..cs.num_inputs {
            cs.enforce(|| "", |lc| lc + Variable::new_unchecked(Index::Input(i)), |lc| lc, |lc| lc);
        }

        Ok(QAP { a: vec![], b: vec![], c: vec![] })
    }
}


#[cfg(test)]
mod tests {
    use super::assignments;
    use std::mem;
    use bls12_381::{ Scalar };

    #[test]
    fn test() {
        let a = assignments::AnalyzeCircuit::<Scalar>::default();
        println!("{}", mem::size_of_val(&a));
    }
}