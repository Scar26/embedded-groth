use core::ops::{AddAssign, Mul, Neg};

use pairing::{ Engine };
use ff::Field;


use crate::{ VerificationKey, VerificationError, Proof };

pub fn verify_proof<E: Engine>(
    proof: &Proof<E>,
    public_inputs: &[E::Fr],
    vk: VerificationKey<E>,
) -> Result<(), VerificationError> {
    if (public_inputs.len() + 1) != vk.ic.len() {
        return Err(VerificationError::InvalidVerifyingKey);
    }

    let mut acc: E::G1 = vk.ic[0].into();

    for (i, b) in public_inputs.iter().zip(vk.ic.iter().skip(1)) {
        acc.add_assign(&(*b * i));
    }
    
    let mut rhs = E::pairing(&proof.a, &proof.b.into());
    rhs.add_assign(E::pairing(&acc.into(), &vk.gamma_g2.mul(E::Fr::one().neg()).into()));
    rhs.add_assign(E::pairing(&proof.c, &vk.delta_g2.mul(E::Fr::one().neg()).into()));

    if E::pairing(&vk.alpha_g1, &vk.beta_g2) == rhs {
        Ok(())
    } else {
        Err(VerificationError::InvalidProof)
    }

}