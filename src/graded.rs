use std::{hash::Hash, clone, collections::{hash_map, HashMap}, fmt::{Debug, Display}, marker::PhantomData, ops::{Add, AddAssign, Sub, SubAssign}};

use crate::{field::{Field, Fp, F2}, matrix::{Matrix, FieldMatrix}};

pub trait Grading : 'static + Clone + Hash + Copy + Debug + Sized + Add<Output=Self> + Sub<Output=Self> + PartialEq + Eq + AddAssign + SubAssign + PartialOrd  {
    fn degree_names() -> Vec<char>;
    fn default_formulas() -> (String, String);
}

impl Grading for i32 {
    fn degree_names() -> Vec<char> {
        vec!['t']
    }
    
    fn default_formulas() -> (String, String) {
        ("t-s".to_string(), "s".to_string())
    }
}

pub type UniGrading = i32;

#[derive(Debug, Clone, Copy, Hash)]
pub struct BiGrading(i32,i32);


pub trait BasisElement : 'static + Debug + Clone {}


pub type BasisIndex<G: Grading> = (G, usize);

// A VectorSpace should be naive / simple, just a list of basis elements!
// Specific modules can implement their own "basis" type which encodes the information they need
pub type VectorSpace<B: BasisElement> = Vec<B>;


// Maybe make this its own type ???
// This is probably fine, as modules will always direct use this type
pub type GradedVectorSpace<G: Grading, B: BasisElement> = HashMap<G, VectorSpace<B>>;

#[derive(Debug, Clone)]
pub struct GradedLinearMap<G: Grading, F: Field, M: Matrix<F>> {
    maps: HashMap<G, M>,
    __: PhantomData<F>,
}

impl<G: Grading, F: Field, M: Matrix<F>> GradedLinearMap<G,F,M> {
    pub fn get_cokernel(&self) -> Self {
       let kernel: HashMap<G, M> = self.maps.iter().map(|(k,v) | {
           todo!("KERNEL SHOULD BE COKERNEL");
           (*k, v.kernel())
        }).collect();
        GradedLinearMap { maps: kernel, __: PhantomData }
    }
    
    pub fn get_kernel(&self) -> Self {
        let kernel: HashMap<G, M> = self.maps.iter().map(|(k,v) | {
            (*k, v.kernel())
         }).collect();
         GradedLinearMap { maps: kernel, __: PhantomData }
    }

    pub fn compose(self, mut rhs: Self) -> Self {
        let compose = self.maps.into_iter().filter_map(|(k,mut v)| {
            match rhs.maps.get_mut(&k) {
                None => { None},
                Some(t) => {
                    Some((k, v.compose(t)))       
                }
            }
        }).collect();
        GradedLinearMap { maps: compose, __: PhantomData }
    }
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