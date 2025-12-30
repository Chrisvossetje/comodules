use std::fmt::Debug;
use deepsize::DeepSizeOf;
use serde::{Deserialize, Serialize};

use crate::{matrix::Matrix, ring::CRing};

pub trait Abelian<R: CRing>: Matrix<R> + Clone + Send + Sync + PartialEq + Debug {
    /// Default should represents a free generator!
    type Generator: Default + Send + Sync + Copy + Clone + Ord + Debug + DeepSizeOf + Serialize + for<'a> Deserialize<'a>;

    fn kernel(&self, domain: &Vec<Self::Generator>, codomain: &Vec<Self::Generator>) -> (Self, Vec<Self::Generator>);
    fn cokernel(&self, codomain: &Vec<Self::Generator>) -> (Self, Self, Vec<Self::Generator>);
    
    fn cohomology(f: &Self, g: &Self, n: &Vec<Self::Generator>, q: &Vec<Self::Generator>) -> (Self, Vec<Self::Generator>);

    // TODO :
    // ASSUME THIS RETURNS Q_INDEX's WHICH ARE ORDERED WRT DOMAIN!!
    // SIMPLE CHECK HERE IS THAT RETURN VEC IS SORTED
    fn kernel_destroyers(&self, domain: &Vec<Self::Generator>, codomain: &Vec<Self::Generator>) -> Vec<usize>;

    fn compose(f: &Self, g: &Self, g_codomain: &Vec<Self::Generator>) -> Self;
}
