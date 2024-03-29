#![allow(unused_imports)]
#![allow(unused_variables)]
use bellman;

// For randomness (during paramgen and proof generation)
use rand::{thread_rng, Rng};
use groth16::{ prover, assignments, Parameters as GrothParams, VerificationKey as GrothVk };
use std::sync::Arc;

// Bring in some tools for using pairing-friendly curves
use pairing::Engine;

use ff::{ Field, PrimeField };

// We're going to use the BLS12-381 pairing-friendly elliptic curve.
use bls12_381::{ Bls12, Scalar };

// We'll use these interfaces to construct our circuit.
use self::bellman::{
    Circuit,
    ConstraintSystem,
    SynthesisError
};

// We're going to use the Groth16 proving system.
use self::bellman::groth16::{
    Proof as BellmanProof,
    generate_random_parameters,
    prepare_verifying_key,
    create_proof,
    verify_proof,
};

// proving that I know x such that x^3 + x + 5 == 35
// Generalized: x^3 + x + 5 == out
#[derive(Clone, Copy)]
pub struct CubeDemo<S: PrimeField> {
    pub x: Option<S>,
}

impl<S: PrimeField> Circuit<S> for CubeDemo<S> {
    fn synthesize<CS: ConstraintSystem<S>>(
        self, 
        cs: &mut CS
    ) -> Result<(), SynthesisError>
    {
        // Flattened into quadratic equations (x^3 + x + 5 == 35): 
        // x * x = tmp_1
        // tmp_1 * x = y
        // y + x = tmp_2
        // tmp_2 + 5 = out
        // Resulting R1CS with w = [one, x, tmp_1, y, tmp_2, out]
        
        // Allocate the first private "auxiliary" variable
        let x_val = self.x;
        let x = cs.alloc(|| "x", || {
            x_val.ok_or(SynthesisError::AssignmentMissing)
        })?;
        
        // Allocate: x * x = tmp_1
        let tmp_1_val = x_val.map(|mut e| {
            e = e.square();
            e
        });
        let tmp_1 = cs.alloc(|| "tmp_1", || {
            tmp_1_val.ok_or(SynthesisError::AssignmentMissing)
        })?;
        // Enforce: x * x = tmp_1
        cs.enforce(
            || "tmp_1",
            |lc| lc + x,
            |lc| lc + x,
            |lc| lc + tmp_1
        );
        
        // Allocate: tmp_1 * x = y
        let x_cubed_val = tmp_1_val.map(|mut e| {
            e.mul_assign(&x_val.unwrap());
            e
        });
        let x_cubed = cs.alloc(|| "x_cubed", || {
            x_cubed_val.ok_or(SynthesisError::AssignmentMissing)
        })?;
        // Enforce: tmp_1 * x = y
        cs.enforce(
            || "x_cubed",
            |lc| lc + tmp_1,
            |lc| lc + x,
            |lc| lc + x_cubed
        );
        
        // Allocating the public "primary" output uses alloc_input
        let out = cs.alloc_input(|| "out", || {
            let mut tmp = x_cubed_val.unwrap();
            tmp.add_assign(&x_val.unwrap());
            tmp.add_assign(&S::from_str_vartime("5").unwrap());
            Ok(tmp)
        })?;    
        // tmp_2 + 5 = out
        // => (tmp_2 + 5) * 1 = out
        cs.enforce(
            || "out",
            |lc| lc + x_cubed + x + (S::from_str_vartime("5").unwrap(), CS::one()),
            |lc| lc + CS::one(),
            |lc| lc + out
        );
        // lc is an inner product of all variables with some vector of coefficients
        // bunch of variables added together with some coefficients
        
        // usually if mult by 1 can do more efficiently
        // x2 * x = out - x - 5
        
        // mult quadratic constraints 
        // 
        
        Ok(())
    }
}

fn main(){
    // This may not be cryptographically safe, use
    // `OsRng` (for example) in production software.
    let mut rng = thread_rng();
    
    println!("Creating parameters...");
    
    // Create parameters for our circuit
    let params = {
        let c = CubeDemo::<Scalar> {
            x: None
        };

        generate_random_parameters::<Bls12, _, _>(c, &mut rng).unwrap()
    };

    println!("params 1: {}, {}, {}, {}, {}", params.l.len(), params.h.len(), params.a.len(), params.b_g1.len(), params.b_g2.len());
    // Prepare the verification key (for proof verification)
    let pvk = prepare_verifying_key(&params.vk);

    println!("Creating proofs...");
    
    // Create an instance of circuit
    let c = CubeDemo::<Scalar> {
        x: Scalar::from_str_vartime("3")
    };

    let assignments = assignments::extract_assignments::<CubeDemo<Scalar>, Bls12>(c).unwrap();
    let (inputsassign, auxassign) = assignments.get_assignments();
    let m = assignments.num_constraints();
    let cap = assignments.qap();
    let r = Scalar::random(&mut rng);
    let s = Scalar::random(&mut rng);

    let bellproof = create_proof(c, &params, r.clone(), s.clone()).unwrap();
    println!("bellman proof: {:?}", bellproof);

    let grothparams = assignments::create_params(params);
    let grothproof = prover::create_proof::<Bls12>(grothparams, inputsassign.as_ref(), auxassign.as_ref(), r, s, cap, m);
    println!("groth proof: {:?}", grothproof);
    let g2bproof = BellmanProof {
        a: grothproof.a.clone(),
        b: grothproof.b.clone(),
        c: grothproof.c.clone(),
    };

    verify_proof(&pvk, &g2bproof, &inputsassign[1..]).unwrap();
}
