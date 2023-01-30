use ff::PrimeField;

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

#[cfg(test)]
mod tests {
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

        assert_eq!(x, a2);

        // inverse fft
        domain.ifft(&worker);
        ifft(a2.as_mut_slice(), &omega, exp);

        x = domain.as_ref().iter().map(|t| t.0).collect();
        assert_eq!(x, a2);
        assert_eq!(avals, a2);
    }
}