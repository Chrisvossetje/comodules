use crate::{field::Field, matrices::flat_matrix::FlatMatrix, rings::univariate_polynomial_ring::UniPolRing};

pub mod morphism;
pub mod cokernel;
pub mod kernel;
pub mod cohomology;
pub mod helper;

type UniPolModule = Vec<Option<u16>>;
#[allow(type_alias_bounds)]
type UniPolMap<F: Field> = FlatMatrix<UniPolRing<F>>;


#[cfg(test)]
mod tests;
