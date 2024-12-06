use std::{clone, collections::{hash_map::{self, IntoIter}, HashMap}, fmt::{Debug, Display}, hash::Hash, marker::PhantomData, ops::{Add, AddAssign, Sub, SubAssign}};

use crate::{field::{Field, Fp, F2}, matrix::{Matrix, FieldMatrix}};

pub trait Grading : 'static + Clone + Hash + Copy + Debug + Sized + Add<Output=Self> + Sub<Output=Self> + PartialEq + Eq + AddAssign + SubAssign + PartialOrd  {
    fn degree_names() -> Vec<char>;
    fn default_formulas() -> (String, String);
    fn export_grade(self) -> Vec<i32>;
}

impl Grading for i32 {
    fn degree_names() -> Vec<char> {
        vec!['t']
    }
    
    fn default_formulas() -> (String, String) {
        ("t-s".to_string(), "s".to_string())
    }

    fn export_grade(self) -> Vec<i32> {
        vec![self]
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
// pub type GradedVectorSpace<G: Grading, B: BasisElement> = HashMap<G, VectorSpace<B>>;

pub struct GradedVectorSpace<G: Grading, B: BasisElement>(
    pub HashMap<G, VectorSpace<B>>,
);

impl<G: Grading, B: BasisElement> IntoIterator for GradedVectorSpace<G,B>{
    type Item = (G,Vec<B>);
    
    type IntoIter = IntoIter<G, Vec<B>>;
    
    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl<G: Grading, B: BasisElement> GradedVectorSpace<G,B> {
    pub fn new() -> Self {
        Self(HashMap::new())
    }


}

#[derive(Debug, Clone)]
pub struct GradedLinearMap<G: Grading, F: Field, M: Matrix<F>> {
    maps: HashMap<G, M>,
    __: PhantomData<F>,
}

impl<G: Grading, F: Field, M: Matrix<F>> GradedLinearMap<G,F,M> {
    pub fn get_cokernel(&self) -> Self {
        todo!("IF DOMAIN DOES NOT HAVE A GRADE THIS IS WRONGGGG");
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
                None => { None },
                Some(t) => {
                    Some((k, v.compose(t)))       
                }
            }
        }).collect();
        GradedLinearMap { maps: compose, __: PhantomData }
    }
    
    pub fn pivots(&self) -> HashMap<G, Vec<(usize,usize)>> {
        self.maps.iter().map(|(k, v)| {
            (*k, v.pivots())
        }).collect()
    }

    pub fn empty() -> Self {
        GradedLinearMap {
            maps: HashMap::new(),
            __: PhantomData,
            
        }
    }


    pub fn zero<B: BasisElement>(domain: GradedVectorSpace<G, B>, codomain: GradedVectorSpace<G, B>) -> Self {
        todo!()
    }

    pub fn codomain_space<B: BasisElement>(&self) -> GradedVectorSpace<G, B> {
        todo!()
    }
}