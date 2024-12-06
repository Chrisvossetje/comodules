use std::{collections::HashMap, sync::Arc};

use crate::{linalg::{field::Field, graded::{BasisElement, BasisIndex, GradedLinearMap, GradedVectorSpace, Grading}, matrix::FieldMatrix}};

use super::comodule::{Comodule, ComoduleMorphism, Tensor};




#[derive(Debug, Clone)]
pub struct kBasisElement {
    name: String,
    generator: bool,
    primitive: Option<usize>,
    generated_index: usize,
}

impl BasisElement for kBasisElement {

}

#[derive(Debug, Clone)]
pub struct kComodule<G: Grading, F: Field> {
    coalgebra: Arc<kComodule<G, F>>,
    space: GradedVectorSpace<G, kBasisElement>,
    coaction: GradedLinearMap<G, F, FieldMatrix<F>>,
    tensor: kTensor<G>,
} 


#[derive(Debug, Clone)]
pub struct kTensor<G: Grading> {
    construct: HashMap<G,G>,
    deconstruct: HashMap<G,G>, 
}

#[derive(Debug, Clone)]
pub struct kComoduleMorphism<G: Grading, F: Field> {
    domain: Arc<kComodule<G, F>>,
    codomain: Arc<kComodule<G, F>>,

    map: GradedLinearMap<G,F, FieldMatrix<F>> // Question: Shouldn't this be a module morphism?
}


impl<G: Grading> Tensor<G> for kTensor<G> {
    fn tensor_to_base() {
        todo!()
    }

    fn base_to_tensor() {
        todo!()
    }
}

impl<G: Grading> Default for kTensor<G> {
    fn default() -> Self {
        Self { construct: Default::default(), deconstruct: Default::default() }
    }
}

impl<G: Grading> kTensor<G> {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn generate<B: BasisElement>(left: &GradedVectorSpace<G, B>, right: &GradedVectorSpace<G, B>) -> Self {
        todo!()
    }
}



impl<G: Grading, F: Field> Comodule<G> for kComodule<G, F> {
    type Element = kBasisElement;

    fn get_generators(&self) -> Vec<(usize, G, Option<String>)> {
        self.space.0.iter().flat_map(|(k,v)| {
            v.iter().filter_map(|b| {
                if b.generator {
                    Some((b.generated_index, *k, Some(b.name.clone())))
                } else {
                    None
                }
            })
        }).collect()
    }
    
    
    fn zero_comodule(comodule: Arc<Self>) -> Self {
        Self {
            coalgebra: comodule.coalgebra.clone(),
            space: GradedVectorSpace::new(),
            coaction: GradedLinearMap::empty(),
            tensor: kTensor::new()
        }
    }
}


impl<G: Grading, F: Field> ComoduleMorphism<G, kComodule<G, F>> for kComoduleMorphism<G, F> {
    fn cokernel(&self) -> Self {
        let cokernel_map = self.map.get_cokernel();
    
        let coker_space = cokernel_map.codomain_space();

        let pivots = cokernel_map.pivots();
        let coalg = self.codomain.coalgebra.as_ref();
        
        let tensor = kTensor::generate(&coalg.space, &coker_space);

        
        todo!()
    }

    fn inject_codomain_to_cofree(&self) -> Self {
        todo!()
    }

    fn zero_morphism(comodule: Arc<kComodule<G, F>>) -> Self {
        let codomain = comodule.clone();
        let zero = Arc::new(kComodule::zero_comodule(comodule));

        // Verify how we want to handle this zero map
        let zero_map = GradedLinearMap::empty();
        Self {
            domain: zero,
            codomain: codomain,
            map: zero_map,
        }
    }

    fn compose(l: Self, r: Self) -> Self {
        let codomain = l.codomain;
        let domain = r.domain;

        let map = l.map.compose(r.map);

        Self {
            domain,
            codomain,
            map,
        }
    }
    
    fn get_structure_lines(&self) -> Vec<(BasisIndex<G>, BasisIndex<G>, usize)> {
        todo!()
    }
    
    fn get_codomain(&self) -> &kComodule<G, F> {
        todo!()
    }
}


