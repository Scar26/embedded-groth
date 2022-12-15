#![allow(unused_imports)]
#![allow(unused_variables)]
extern crate bellman;
extern crate pairing;
extern crate rand;
extern crate ff;

// For randomness (during paramgen and proof generation)
use self::rand::{thread_rng, Rng};
use groth16::assignments;


// Bring in some tools for using pairing-friendly curves
use self::pairing::Engine;

use ff::PrimeField;

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
    Proof,
    generate_random_parameters,
    prepare_verifying_key,
    create_random_proof,
    verify_proof,
};

// proving that I know x such that x^3 + x + 5 == 35
// Generalized: x^3 + x + 5 == out
pub struct CubeDemo<S: PrimeField> {
    pub x: Option<S>,
}

// Algorithm
// proof = proof_generate(params, inputs) 
// 
// 
//
//

impl <S: PrimeField> Circuit<S> for CubeDemo<S> {
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
    let _  = assignments::ExtractAssignments{vec![], vec![]};
    let rng = &mut thread_rng();
    
    println!("Creating parameters...");
    
    // Create parameters for our circuit
    let params = {
        let c = CubeDemo::<Scalar> {
            x: None
        };

        generate_random_parameters::<Bls12, _, _>(c, rng).unwrap()
    };

    println!("params 1: {}, {}, {}, {}, {}", params.l[0], params.h.len(), params.a.len(), params.b_g1.len(), params.b_g2.len());
    // Prepare the verification key (for proof verification)
    let pvk = prepare_verifying_key(&params.vk);

    println!("Creating proofs...");
    
    // Create an instance of circuit
    let c = CubeDemo::<Scalar> {
        x: Scalar::from_str_vartime("3")
    };
    
    // Create a groth16 proof with our parameters.
    let proof = create_random_proof(c, &params, rng).unwrap();
        
    if let Ok(_) = verify_proof(
        &pvk,
        &proof,
        &[Scalar::from_str_vartime("35").unwrap()]
    ) {
        println!("proof verified")
    } else {
        println!("verification failed")
    }


}
