use std::sync::Arc;

use ahash::HashMap;
use algebra::ring::CRing;
use deepsize::DeepSizeOf;

use crate::grading::{Grading, OrderedGrading};

pub trait Coalgebra: DeepSizeOf + Clone {}

pub trait Comodule<G: Grading, C: Coalgebra>: DeepSizeOf + Clone {
    fn fp_comodule(coalgebra: &C, degree: G) -> Self;
}

pub trait CofreeComodule<G: Grading, C: Coalgebra>: DeepSizeOf + Clone {
    type Generator: Default;

    // This is not the correct type yet ???
    fn zero_comodule() -> Self;
    fn get_generators(&self) -> Vec<(usize, G, Option<String>)>;
    fn direct_sum(&mut self, other: &mut Self);
    fn cofree_comodule(
        coalgebra: &C,
        index: usize,
        grade: G,
        limit: G,
        generator: Self::Generator,
    ) -> Self;
}

// pub trait CofreeComodule<G: Grading>: Sized {

// }

pub trait ComoduleMorphism<G: Grading, C: Coalgebra>: Sized + DeepSizeOf {
    type CofreeComodule: CofreeComodule<G, C>;
    type Comodule: Comodule<G, C>;
    type BaseRing: CRing;

    fn cokernel(&self, coalgebra: &C, codomain: &Self::CofreeComodule) -> (Self, Self::Comodule);

    fn inject_codomain_to_cofree(
        coalgebra: &C,
        comodule: &Self::Comodule,
        limit: G,
    ) -> (Self, Self::CofreeComodule)
    where
        G: OrderedGrading;

    fn zero_morphism(comodule: &Self::Comodule) -> Self;

    // domain l == codomain r, l \circ r
    fn compose(l: &Self, r: &Self, codomain: &Self::CofreeComodule) -> Self;

    fn direct_sum(
        a: &mut Self,
        b: &mut Self,
        a_codom: &mut Self::CofreeComodule,
        b_codom: &mut Self::CofreeComodule,
    );

    // fn get_codomain(&self) -> Arc<M>;

    /// (s, gen_index) uniquely defines a generator of Ext
    /// in a specific morphism we only need to know its gen_index
    /// in the resolution we add the s
    /// (from_dot, to_dot, value, line_type)
    fn get_structure_lines(
        &self,
        coalgebra: &C,
        domain: &Self::CofreeComodule,
        codomain: &Self::CofreeComodule,
    ) -> Vec<((usize, G, usize), (usize, G, usize), Self::BaseRing, String)>;
}
