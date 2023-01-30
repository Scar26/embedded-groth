use core::ops::AddAssign;

use pairing::Engine;
use pairing::group::{ GroupOps, GroupOpsOwned };
use crate::{ Proof, Parameters, poly };

#[cfg(not(any(test, feature = "std")))]
use alloc::vec::Vec;

pub fn create_proof<E: Engine>(
    params: Parameters<E>,
    inputs: &[E::Fr],
    aux: &[E::Fr],
    input_constraints: &[bool],
    aux_constraints: &[bool],
    r: E::Fr,
    s: E::Fr,
    at: Vec<Vec<E::Fr>>,
    bt: Vec<Vec<E::Fr>>,
    ct: Vec<Vec<E::Fr>>,
) -> Proof<E>
where
    E: Engine,
{
    
    unimplemented!()

}