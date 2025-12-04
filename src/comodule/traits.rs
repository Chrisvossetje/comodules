use std::sync::Arc;

use crate::{
    basiselement::BasisElement,
    grading::{Grading, OrderedGrading}
};

pub trait Comodule<G: Grading>: Sized {
    type Element: BasisElement;
    type Coalgebra;
    type Morphism: ComoduleMorphism<G, Self> + Clone;
    type BaseRing: std::fmt::Debug + Clone;
    type Generator;

    fn zero_comodule(coalgebra: Arc<Self::Coalgebra>) -> Self;

    // This is not the correct type yet
    fn get_generators(&self) -> Vec<(usize, G, Option<String>)>;

    fn fp_comodule(coalgebra: Arc<Self::Coalgebra>, degree: G) -> Self;

    fn direct_sum(&mut self, other: &mut Self);

    fn cofree_comodule(coalgebra: Arc<Self::Coalgebra>, index: usize, grade: G, limit: G, generator: Self::Generator) -> Self;
}

pub trait ComoduleMorphism<G: Grading, M: Comodule<G>> {
    fn cokernel(&self) -> Self;
    fn inject_codomain_to_cofree(&self, limit: G) -> Self 
    where
        G: OrderedGrading; // Question: Shouldn't 'codomain' be 'cokernel'/'comodule'?
    fn zero_morphism(comodule: Arc<M>) -> Self;

    // domain l == codomain r, l \circ r
    fn compose(l: &Self, r: &Self) -> Self;

    fn get_codomain(&self) -> Arc<M>;

    /// (s, gen_index) uniquely defines a generator of Ext
    /// in a specific morphism we only need to know its gen_index
    /// in the resolution we add the s
    /// (from_dot, to_dot, value, line_type)
    fn get_structure_lines(&self) -> Vec<((usize, G, usize), (usize, G, usize), M::BaseRing, String)>;
}

