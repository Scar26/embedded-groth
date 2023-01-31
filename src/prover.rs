use core::ops::{SubAssign, MulAssign, AddAssign, Mul, Add};

use pairing::Engine;
use pairing::group::Group;
use crate::{ Proof, Parameters, poly::*, QAP };
use ff::{Field, PrimeField};

#[cfg(not(any(test, feature = "std")))]
use alloc::vec::Vec;

pub fn create_proof<E: Engine>(
    params: Parameters<E>,
    inputs: &[E::Fr],
    aux: &[E::Fr],
    r: E::Fr,
    s: E::Fr,
    qap: QAP<E::Fr>,
    num_constraints: usize
) -> Proof<E>
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

    assert_eq!(aux.len(), params.l.len());
    let l = params.l.iter()
        .zip(aux.iter())
        .fold(E::G1::identity(), |acc, (x, y)| acc.add(x.mul(y)));

    let augmented_inputs: Vec<&E::Fr> = inputs.iter().chain(aux.iter()).collect();

    assert_eq!(params.a_g1.len(), qap.a_constraints.len());
    let at_g1 = params.a_g1.iter()
        .zip(qap.a_constraints.iter())
        .fold(E::G1::identity(), |acc, (x, y)| acc.add(x.mul(augmented_inputs[*y])));

    assert_eq!(params.b_g1.len(), qap.b_constraints.len());
    let bt_g1 = params.b_g1.iter()
        .zip(qap.b_constraints.iter())
        .fold(E::G1::identity(), |acc, (x, y)| acc.add(x.mul(augmented_inputs[*y])));
    
    let bt_g2 = params.b_g2.iter()
        .zip(qap.b_constraints.iter())
        .fold(E::G2::identity(), |acc, (x, y)| acc.add(x.mul(augmented_inputs[*y])));

    let mut a = E::G1::identity();
    a.add_assign(params.vk.alpha_g1);
    a.add_assign(at_g1);
    a.add_assign(params.vk.delta_g1.mul(r));

    let mut b = E::G2::identity();
    b.add_assign(params.vk.beta_g2);
    b.add_assign(bt_g2);
    b.add_assign(params.vk.delta_g2.mul(s));

    let mut c = E::G1::identity();
    c.add_assign(bt_g1);
    c.add_assign(params.vk.delta_g1.mul(s));
    c.add_assign(params.vk.beta_g1);
    c.mul_assign(r);
    c.add_assign(h);
    c.add_assign(l);
    c.add_assign(a.mul(s));
    c.sub_assign(params.vk.delta_g1.mul(r.mul(s)));
    
    
    Proof {
        a: a.into(),
        b: b.into(),
        c: c.into(),
    }
}