use std::fmt::Debug;
use crate::{matrix::Matrix, ring::CRing};

pub trait Abelian<R: CRing>: Matrix<R> + Clone + Send + Sync + PartialEq + Debug {
    type Module: Default + Send + Sync;

    fn kernel(&self, domain: &Self::Module, codomain: &Self::Module) -> (Self, Self::Module);
    fn cokernel(&self, codomain: &Self::Module) -> (Self, Self, Self::Module);
    
    fn cohomology(f: &Self, g: &Self, n: &Self::Module, q: &Self::Module) -> (Self, Self::Module);

    fn kernel_generators(&self, domain: &Self::Module, codomain: &Self::Module) -> Vec<usize>;

    fn compose(f: &Self, g: &Self, g_codomain: &Self::Module) -> Self;
}
