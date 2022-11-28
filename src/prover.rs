use pairing::Engine;
use rand_core::RngCore;
pub struct GetAssignment<E: Engine> {
    pub a: E::Fr, 
}

// impl<E: Engine> ConstraintSystem<E::Fr> for GetAssignment<E> {
//     type Root = Self;

//     fn enforce<A, AR, LA, LB, LC>(&mut self, annotation: A, a: LA, b: LB, c: LC)
//         where
//             A: FnOnce() -> AR,
//             AR: Into<alloc::string::String>,
//             LA: FnOnce(bellman::LinearCombination<E::Fr>) -> bellman::LinearCombination<E::Fr>,
//             LB: FnOnce(bellman::LinearCombination<E::Fr>) -> bellman::LinearCombination<E::Fr>,
//             LC: FnOnce(bellman::LinearCombination<E::Fr>) -> bellman::LinearCombination<E::Fr>
//     {
//         // Not needed since constraints are already aggregated in the CRS
//     }
// }


// pub fn generate_random_proof<C, E, R, const P: usize, const A: usize, const M: usize>(
//     circuit: C,
//     parameters: Parameters<E, P, A, M>,
//     mut rng: &mut R,
// ) -> ()
// where
//     E: Engine,
//     C: Circuit<E::Fr>,
//     R: RngCore,
// {
//     let _a = LinearCombination::<E::Fr>::zero();
// }