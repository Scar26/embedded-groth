use pairing::Engine;

use crate::{VerificationError, Proof, Parameters};

pub fn verify_proof<E: Engine>(
    proof: &Proof<E>,
    params: &Parameters<E>,
    public_inputs: &[E::Fr],
) -> Result<(), VerificationError> {
    if (public_inputs.len() + 1) != params.vk.ic.len() {
        return Err(VerificationError::InvalidVerifyingKey);
    }

    unimplemented!()
}