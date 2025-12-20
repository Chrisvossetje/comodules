// use fp::prime::P2;

// use crate::{matrix::Matrix, rings::finite_fields::F2};

// type M = fp::matrix::Matrix;

// impl Matrix<F2> for M {
//     type UnderlyingRowType = F2;

//     fn zero(domain: usize, codomain: usize) -> Self {
//         M::new(P2, codomain, domain)   
//     }

//     fn get(&self, domain: usize, codomain: usize) -> F2 {
//         self[(codomain, domain)]
//     }

//     fn set(&mut self, domain: usize, codomain: usize, r: F2) {
//         todo!()
//     }

//     fn add_at(&mut self, domain: usize, codomain: usize, r: F2) {
//         todo!()
//     }

//     fn get_row(&self, codomain: usize) -> &[Self::UnderlyingRowType] {
//         todo!()
//     }

//     fn set_row(&mut self, codomain: usize, row: &[Self::UnderlyingRowType]) {
//         todo!()
//     }

//     fn eval_vector(&self, vector: &[F2]) -> Vec<F2> {
//         todo!()
//     }

//     fn compose(&self, rhs: &Self) -> Self {
//         todo!()
//     }

//     fn transpose(&self) -> Self {
//         todo!()
//     }

//     fn domain(&self) -> usize {
//         todo!()
//     }

//     fn codomain(&self) -> usize {
//         todo!()
//     }

//     fn vstack(&mut self, other: &mut Self) {
//         todo!()
//     }

//     fn block_sum(&mut self, other: &Self) {
//         todo!()
//     }

//     fn extend_one_row(&mut self) {
//         todo!()
//     }
// }