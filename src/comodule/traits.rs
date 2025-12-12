use std::sync::Arc;

use algebra::ring::CRing;
use deepsize::DeepSizeOf;

use crate::grading::{Grading, OrderedGrading};

pub trait Comodule<G: Grading, R: CRing>: DeepSizeOf + Clone {
    type Coalgebra;

    fn fp_comodule(coalgebra: Arc<Self::Coalgebra>, degree: G) -> Self;
}

pub trait CofreeComodule<G: Grading>: DeepSizeOf + Clone {
    type Coalgebra;
    type Generator: Default;

    // This is not the correct type yet ???
    fn zero_comodule(coalgebra: Arc<Self::Coalgebra>) -> Self;
    fn get_generators(&self) -> Vec<(usize, G, Option<String>)>;
    fn direct_sum(&mut self, other: &mut Self);
    fn cofree_comodule(
        coalgebra: Arc<Self::Coalgebra>,
        index: usize,
        grade: G,
        limit: G,
        generator: Self::Generator,
    ) -> Self;
}

// pub trait CofreeComodule<G: Grading>: Sized {

// }

pub trait ComoduleMorphism<G: Grading>: Sized + DeepSizeOf {
    type CofreeComodule: CofreeComodule<G>;
    type Comodule: Comodule<G, Self::BaseRing>;
    type BaseRing: CRing;

    fn cokernel(&self, codomain: &Self::CofreeComodule) -> (Self, Self::Comodule);

    fn inject_codomain_to_cofree(
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
        domain: &Self::CofreeComodule,
        codomain: &Self::CofreeComodule,
    ) -> Vec<((usize, G, usize), (usize, G, usize), Self::BaseRing, String)>;
}
