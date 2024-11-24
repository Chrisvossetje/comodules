use std::{hash::Hash, clone, collections::{hash_map, HashMap}, fmt::{Debug, Display}, marker::PhantomData, ops::{Add, AddAssign, Sub, SubAssign}};

use crate::{fp::{Field, Fp, F2}, matrix::{Matrix, FieldMatrix}};

pub trait Grading : 'static + Clone + Hash + Copy + Debug + Sized + Add<Output=Self> + Sub<Output=Self> + PartialEq + Eq + AddAssign + SubAssign + PartialOrd  {}
impl Grading for i32 {}
impl Grading for (i32,i32) {}



pub trait BasisElement : 'static + Debug + Clone {}


pub type BasisIndex<G: Grading> = (G, usize);

// A VectorSpace should be naive / simple
pub type VectorSpace<B: BasisElement> = Vec<B>;


// Maybe make this its own type ???
pub type GradedVectorSpace<G: Grading, B: BasisElement> = HashMap<G, VectorSpace<B>>;


pub struct GradedLinearMap<G: Grading, F: Field, M: Matrix<F>> {
    maps: HashMap<G, M>,
    __: PhantomData<F>,
}

impl<G: Grading, F: Field, M: Matrix<F>> GradedLinearMap<G,F,M> {
    fn get_cokernel(&self) -> Self {
       let kernel: HashMap<G, M> = self.maps.iter().map(|(k,v) | {
           todo!("KERNEL SHOULD BE COKERNEL");
           (*k, v.kernel())
        }).collect();
        GradedLinearMap { maps: kernel, __: PhantomData }
    }
    
    fn get_map(&self, grade: G) -> Option<&M> {
        self.maps.get(&grade)
    }
    
    fn get_kernel(&self) -> Self {
        let kernel: HashMap<G, M> = self.maps.iter().map(|(k,v) | {
            (*k, v.kernel())
         }).collect();
         GradedLinearMap { maps: kernel, __: PhantomData }
    }

    fn transpose() {
        todo!()
    }
}



// EXAMPLE OF A kt MODULE 

#[derive(Debug, Clone)]
pub struct ktBasisElement {

}

pub struct ktModule {
    space: GradedVectorSpace<(i32, i32), ktBasisElement>,
    t: GradedLinearMap<(i32, i32), F2, FieldMatrix<F2>>
}




// //
// // BasisElement
// //

// #[derive(Clone, Copy, Debug, PartialEq, Eq)]
// pub struct BasisElement<G: Grading> {
//     grading: G,
//     index: usize
// }

// impl<G: Grading> GradedBasisElement<G> for BasisElement<G> {
//     fn get_grading(&self) -> G {
//         self.grading
//     }
// }


// //
// // VectorSpace
// //

// pub struct VectorSpace<G: Grading, F: Field> {
//     basis:  HashMap<G, Vec<BasisElement<G>>>,
//     phantom: std::marker::PhantomData<F>, // To make the compiler happy, since it is throwing errors that we don't actually use F
// }


// impl<G: Grading, F: Field> GradedVectorSpace<G, F> for VectorSpace<G,F> {
    
//     fn get_basis(&self) -> Vec<impl GradedBasisElement<G>> {
//         let mut total_basis: Vec<BasisElement<G>> = Vec::new();
//         for (_, basis) in self.basis.iter() {
//             total_basis.extend(basis.iter());
//         }
//         return total_basis;
//     }
// }


// //
// // LinearMap
// //

// pub struct LinearMap<G: Grading, F: Field> {
//     domain: VectorSpace<G, F>,
//     codomain: VectorSpace<G, F>,
//     map: HashMap<G,Vec<Vec<F>>> // Grading -> Matrix
// }

// impl<G: Grading, F: Field> GradedLinearMap<G, F> for LinearMap<G, F> {
//     fn get_domain(&self) -> &impl GradedVectorSpace<G, F> {
//         &self.domain
//     }

//     fn get_codomain(&self) -> &impl GradedVectorSpace<G, F> {
//         &self.codomain
//     }
    
//     fn get_cokernel(&self) -> Self {
//         unimplemented!()
//     }

// }