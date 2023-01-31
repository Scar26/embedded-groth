use core::ops::{SubAssign, MulAssign, AddAssign, Mul};

use pairing::Engine;
use pairing::group::Group;
use crate::{ Proof, Parameters, poly::*, QAP };
use ff::{Field, PrimeField};

#[cfg(not(any(test, feature = "std")))]
use alloc::vec::Vec;

pub struct ProvingError {}

pub fn create_proof<E: Engine>(
    params: Parameters<E>,
    inputs: &[E::Fr],
    aux: &[E::Fr],
    _r: E::Fr,
    _s: E::Fr,
    qap: QAP<E::Fr>,
    num_constraints: usize
) -> Result<Proof<E>, ProvingError>
{

    fn eval<S: PrimeField>(
        input_assignment: &[S],
        aux_assignment: &[S],
        output: &mut [S],
        input: Vec<(usize, Vec<(S, usize)>)>,
        p: usize,
    ) {
        for (i, v) in input.into_iter() {
            if i < p {
                for (mut x, u) in v.into_iter() {
                    x.mul_assign(input_assignment[i]);
                    output[u].add_assign(x);
                }
            } else {
                for (mut x, u) in v.into_iter() {
                    x.mul_assign(aux_assignment[i - p]);
                    output[u].add_assign(x);
                }
            }
        }
    }

    let h = {
        let (omega, m, exp): (E::Fr, usize, u32) = fft_params(num_constraints);
        let mut at = vec![E::Fr::zero(); m];
        let mut bt = vec![E::Fr::zero(); m];
        let mut ct = vec![E::Fr::zero(); m];
    
        eval(inputs, aux, &mut at, qap.a, inputs.len());
        eval(inputs, aux, &mut bt, qap.b, inputs.len());
        eval(inputs, aux, &mut ct, qap.c, inputs.len());
    
        coset_mul_assign(&mut at, bt);
        ifft(&mut ct, &omega, exp);
        coset_fft(&mut ct, &omega, exp);
        sub_eval_domain(&mut at, ct);

        let zinv = {
            let mut t = <E::Fr as PrimeField>::multiplicative_generator();
            t = t.pow_vartime(&[at.len() as u64]);
            t.sub_assign(&E::Fr::one());
            t.invert().unwrap()
        };

        for x in at.iter_mut() {
            x.mul_assign(&zinv);
        }

        icoset_fft(&mut at, &omega, exp);

        let mut acc = E::G1::identity();

        at.truncate(at.len() - 1);
        for (i, x) in at.iter().enumerate() {
            let t = params.h[i].mul(x);
            acc.add_assign(t);
        }

        acc
    };

    unimplemented!()
}