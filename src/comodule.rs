use std::{cell::Ref, ops::Mul, primitive, rc::Rc, sync::Arc};

use crate::{field::Field, graded::{BasisElement, BasisIndex, GradedLinearMap, GradedVectorSpace, Grading, UniGrading}, matrix::{FieldMatrix, Matrix}};

pub trait Comodule<G: Grading> {
    type Element: BasisElement;

    
    fn zero_comodule(comodule: RefType<Self>) -> Self;
    
    // This is not the correct type yet
    fn get_generators(&self) -> Vec<BasisIndex<G>>;
}

pub trait ComoduleMorphism<G: Grading, M: Comodule<G>> {
    fn cokernel(&self) -> Self;
    fn inject_codomain_to_cofree(&self) -> Self; // Question: Shouldn't 'codomain' be 'cokernel'/'comodule'?

    fn zero_morphism(comodule: RefType<M>) -> Self;

    // codomain r == codomain l, l \circ r
    fn compose(l: Self, r: Self) -> Self;

    fn get_structure_lines(&self) -> Vec<(BasisIndex<G>, BasisIndex<G>, usize)>; // Question: Should usize here be Field?
}







pub type RefType<'a, T> = Arc<T>;


#[derive(Debug, Clone)]
pub struct kBasisElement {
    name: String,
    generator: bool,
    primitive: Option<usize>,
    generated_index: usize,
}

impl BasisElement for kBasisElement {

}

pub struct kComodule<'a, G: Grading, F: Field> {
    coalgebra: RefType<'a, GradedVectorSpace<G, kBasisElement>>,
    space: GradedVectorSpace<G, kBasisElement>,
    coaction: GradedLinearMap<G, F, FieldMatrix<F>>,
    // tensor: kTensor,
} 

pub struct kComoduleMorphism<'a, G: Grading, F: Field> {
    domain: RefType<'a, kComodule<'a, G, F>>,
    codomain: RefType<'a, kComodule<'a, G, F>>,

    map: GradedLinearMap<G,F, FieldMatrix<F>> // Question: Shouldn't this be a module morphism?
}

impl<G: Grading, F: Field> Comodule<G> for kComodule<'_, G, F> {
    type Element = kBasisElement;

    fn get_generators(&self) -> Vec<BasisIndex<G>> {
        todo!()
    }
    
    
    fn zero_comodule(comodule: RefType<Self>) -> Self {
        Self {
            coalgebra: comodule.coalgebra.clone(),
            space: GradedVectorSpace::new(),
            coaction: GradedLinearMap::empty(),
        }
    }
}


impl<G: Grading, F: Field> ComoduleMorphism<G, kComodule<'_, G, F>> for kComoduleMorphism<'_, G, F> {
    fn cokernel(&self) -> Self {
        let cokernel = self.map.get_cokernel();
        
        

    }

    fn inject_codomain_to_cofree(&self) -> Self {
        todo!()
    }

    fn zero_morphism(comodule: RefType<kComodule<G, F>>) -> Self {
        let codomain = comodule.clone();
        let zero = RefType::new(kComodule::zero_comodule(comodule));

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
}






// pub struct BasicComodule<M: Module, MM: ModuleMorphism, H: HopfAlgebra> {
//     underlying_module: M,
//     underlying_hopf_algebra: H,

//     coaction: MM,
// }


// impl<M: Module, MM: ModuleMorphism, H: HopfAlgebra> Comodule for BasicComodule<M, MM, H> {

//     fn get_underlying_module(&self) -> &impl Module {
//         &self.underlying_module
//     }

//     fn get_hopf_algebra(&self) -> &impl HopfAlgebra {
//         &self.underlying_hopf_algebra
//     }

//     fn get_cogenerating_module(&self) -> impl ModuleMorphism {
//         todo!()
//     }

//     fn cofree_comodule(hopf: impl HopfAlgebra, module: impl Module) -> impl Comodule {
//         todo!()
//     }

// }


// pub struct BasicComoduleMorphism<M: Module, MM: ModuleMorphism, H: HopfAlgebra> {
//     domain: BasicComodule<M, MM, H>,
//     codomain: BasicComodule<M, MM, H>,

//     morphism: MM
// }

// impl <M: Module, MM: ModuleMorphism, H: HopfAlgebra> ComoduleMorphism for BasicComoduleMorphism<M, MM, H> {
//     fn get_underlying_morphism(&self) -> &impl ModuleMorphism {
//         &self.morphism
//     }

//     fn get_domain(&self) -> &impl Comodule {
//         &self.domain
//     }

//     fn get_codomain(&self) -> &impl Comodule {
//         &self.codomain
//     }

//     fn compute_cokernel(&self) -> impl ComoduleMorphism {
//         let mod_coker = &self.morphism.get_cokernel();

//         //    f      q
//         // C  ->  D  ->  K
//         // |      |      |
//         // v      v      v
//         //CxC -> DxD -> KxK

//         // let k be a basis element in K
//         // q is surjective, thus we get an element b in D which maps to k
//         todo!()
        
//     }
    
//     fn get_zero_morphism_to(comod: BasicComodule<M, MM, H>) -> Self {
//         todo!()
//     }
// }
