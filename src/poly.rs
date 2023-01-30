use ff::PrimeField;
#[cfg(not(any(test, feature = "std")))]
use alloc::vec::Vec;

pub fn fft_params<S: PrimeField>(l: usize) -> (S, usize, u32) {
    let mut m = 1;
    let mut exp = 0;

    while m < l {
        m = m << 1;
        exp += 1;
    }

    let mut omega = S::root_of_unity();
    for _ in exp..S::S {
        omega = omega.square();
    }

    (omega, m, exp)
}

pub fn fft<S: PrimeField>(a: &mut [S], omega: &S, exp: u32) {
    fn bitreverse(mut n: u32, l: u32) -> u32 {
        let mut r = 0;
        for _ in 0..l {
            r = (r << 1) | (n & 1);
            n >>= 1;
        }
        r
    }

    let n = a.len() as u32;
    assert_eq!(n, 1 << exp);

    for k in 0..n {
        let rk = bitreverse(k, exp);
        if k < rk {
            a.swap(rk as usize, k as usize);
        }
    }

    let mut m = 1;
    for _ in 0..exp {
        let w_m = omega.pow_vartime(&[u64::from(n / (2 * m))]);

        let mut k = 0;
        while k < n {
            let mut w = S::one();
            for j in 0..m {
                let mut t = a[(k + j + m) as usize];
                t.mul_assign(&w);
                let mut tmp = a[(k + j) as usize];
                tmp.sub_assign(&t);
                a[(k + j + m) as usize] = tmp;
                a[(k + j) as usize].add_assign(&t);
                w.mul_assign(&w_m);
            }
            k += 2 * m;
        }

        m *= 2;
    }
}

pub fn ifft<S: PrimeField>(a: &mut [S], omega: &S, exp: u32) {
    fft(a, &omega.invert().unwrap(), exp);
    let minv = S::from(a.len() as u64).invert().unwrap();
    for i in a {
        i.mul_assign(&minv);
    }
}

// First argument is modified in place to the multiplication result
// Second argument is modified in place to the evaluation domain
pub fn multiply_coefficient_domain<S: PrimeField>(a: &mut Vec<S>, b: &mut Vec<S>) {
    let (omega, m, exp): (S, usize, u32) = fft_params(a.len() + b.len());
    a.resize(m, S::zero());
    b.resize(m, S::zero());
    fft(a.as_mut_slice(), &omega, exp);
    fft(b.as_mut_slice(), &omega, exp);

    for (x, y) in a.iter_mut().zip(b.iter()) {
        x.mul_assign(y);
    }
}

// First argument is modified in place to the addition result
pub fn add_coefficient_domain<S: PrimeField>(a: &mut Vec<S>, b: &Vec<S>) {
    assert_eq!(a.len(), b.len());
    for (x, y) in a.iter_mut().zip(b.iter()) {
        x.add_assign(y);
    }
}

#[cfg(test)]
mod tests {
    use core::ops::{AddAssign, MulAssign};

    use bls12_381::{ Scalar as BlsScalar};
    use rand::thread_rng;
    use ff::Field;
    use bellman::domain::{ EvaluationDomain, Scalar };
    use bellman::multicore::Worker;


    use super::*;

    #[test]
    fn inverse_correctness() {
        let mut rng = thread_rng();
        let a: Vec<Scalar<BlsScalar>> = (0..32).map(|_| Scalar(BlsScalar::random(&mut rng))).collect();
        let avals: Vec<BlsScalar> = a.iter().map(|t| t.0).collect();
        let mut a2: Vec<BlsScalar> = a.iter().map(|t| t.0).collect();
        
        // fft with bellman
        let mut domain = EvaluationDomain::from_coeffs(a.clone()).unwrap();
        let worker = Worker::new();
        domain.fft(&worker);
        let mut x: Vec<BlsScalar> = domain.as_ref().iter().map(|t| t.0).collect();

        // fft with local functions
        let (omega, m, exp): (BlsScalar, usize, u32) = fft_params(a.len());
        fft(a2.as_mut_slice(), &omega, exp);

        // FFT outputs match up with the bellman evaluation domain
        assert_eq!(x, a2);

        // inverse fft
        domain.ifft(&worker);
        ifft(a2.as_mut_slice(), &omega, exp);
        x = domain.as_ref().iter().map(|t| t.0).collect();

        // iFFT outputs match up with the bellman evaluation domain
        assert_eq!(x, a2);

        // FFT i-FFT cancel out
        assert_eq!(avals, a2);
    }

    #[test]
    fn polynomial_arithmetic() {
        let mut rng = thread_rng();
        let mut a: Vec<BlsScalar> = (0..32).map(|_| BlsScalar::random(&mut rng)).collect();
        let mut b: Vec<BlsScalar> = (0..32).map(|_| BlsScalar::random(&mut rng)).collect();

        let mut naive = vec![BlsScalar::zero(); 64];
        for (i, x) in a.iter().enumerate() {
            for (j, y) in b.iter().enumerate() {
                let mut prod = x.clone();
                prod.mul_assign(y);
                naive[i+j].add_assign(prod)
            }
        }

        multiply_coefficient_domain(&mut a, &mut b);
        let (omega, _, exp): (BlsScalar, usize, u32) = fft_params(a.len());
        ifft(&mut a, &omega, exp);

        assert_eq!(naive, a);
    }
}