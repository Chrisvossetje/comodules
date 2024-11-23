use std::{clone, collections::HashMap, fmt::{Debug, Display}, ops::{Add, AddAssign, Sub, SubAssign}};

use crate::{fp::Field, matrix::Matrix};

pub trait Grading : 'static + Clone + Copy + Debug + Sized + Add<Output=Self> + Sub<Output=Self> + PartialEq + Eq + AddAssign + SubAssign + PartialOrd + Display {}

pub trait GradedBasisElement<G: Grading> : 'static + Clone + Copy { 
    fn get_grading(&self) -> G;
}


pub trait GradedVectorSpace<G: Grading, F: Field> {
    fn get_basis(&self) -> Vec<impl GradedBasisElement<G>>;
}


pub trait GradedLinearMap<G: Grading, F: Field> {
    fn get_domain(&self) -> &impl GradedVectorSpace<G, F>;
    fn get_codomain(&self) -> &impl GradedVectorSpace<G, F>;

    fn get_cokernel(&self) -> Self; //&impl GradedLinearMap<G, F>;
}




impl Grading for isize {}



//
// BasisElement
//

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BasisElement<G: Grading> {
    grading: G,
    index: usize
}

impl<G: Grading> GradedBasisElement<G> for BasisElement<G> {
    fn get_grading(&self) -> G {
        self.grading
    }
}


//
// VectorSpace
//

pub struct VectorSpace<G: Grading, F: Field> {
    basis:  HashMap<G, Vec<BasisElement<G>>>,
    phantom: std::marker::PhantomData<F>, // To make the compiler happy, since it is throwing errors that we don't actually use F
}


impl<G: Grading, F: Field>  GradedVectorSpace<G, F> for VectorSpace<G,F> {
    
    fn get_basis(&self) -> Vec<impl GradedBasisElement<G>> {
        let mut total_basis: Vec<BasisElement<G>> = Vec::new();
        for (_, basis) in self.basis.iter() {
            total_basis.extend(basis.iter());
        }
        return total_basis;
    }
}


//
// LinearMap
//

pub struct LinearMap<G: Grading, F: Field> {
    domain: VectorSpace<G, F>,
    codomain: VectorSpace<G, F>,
    map: HashMap<G,Vec<Vec<F>>> // Grading -> Matrix
}

impl<G: Grading, F: Field> GradedLinearMap<G, F> for LinearMap<G, F> {
    fn get_domain(&self) -> &impl GradedVectorSpace<G, F> {
        &self.domain
    }

    fn get_codomain(&self) -> &impl GradedVectorSpace<G, F> {
        &self.codomain
    }
    
    fn get_cokernel(&self) -> Self {
        unimplemented!()
    }

}