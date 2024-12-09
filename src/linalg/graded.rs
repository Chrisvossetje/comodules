use std::{clone, collections::{hash_map::{self, IntoIter}, HashMap}, fmt::{Debug, Display}, hash::Hash, marker::PhantomData, ops::{Add, AddAssign, Sub, SubAssign}};

use super::{field::Field, matrix::Matrix};

pub trait Grading : 'static + Clone + Hash + Copy + Debug + Sized + Add<Output=Self> + Sub<Output=Self> + PartialEq + Eq + AddAssign + SubAssign + PartialOrd  {
    fn degree_names() -> Vec<char>;
    fn default_formulas() -> (String, String);
    fn export_grade(self) -> Vec<i32>;

    fn zero() -> Self;
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
    
    fn zero() -> Self {
        0
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

#[derive(Debug, Clone)]
pub struct GradedVectorSpace<G: Grading, B: BasisElement>(
    pub HashMap<G, VectorSpace<B>>,
);

#[derive(Debug, Clone)]
pub struct GradedLinearMap<G: Grading, F: Field, M: Matrix<F>> {
    pub maps: HashMap<G, M>,
    __: PhantomData<F>,
}

impl<G: Grading, B: BasisElement> GradedVectorSpace<G,B> {
    pub fn new() -> Self {
        Self(HashMap::new())
    }
}


impl<G: Grading, B: BasisElement> From<HashMap<G, Vec<B>>> for GradedVectorSpace<G,B> {
    fn from(value: HashMap<G, Vec<B>>) -> Self {
        Self(value)
    }
}

impl<G: Grading, F: Field, M: Matrix<F>> From<HashMap<G, M>> for GradedLinearMap<G,F,M> {
    fn from(value: HashMap<G, M>) -> Self {
        Self {
            maps: value,
            __: PhantomData,
        }
    }
}

impl<G: Grading, F: Field, M: Matrix<F>> GradedLinearMap<G,F,M> {
    pub fn get_cokernel(&self) -> Self {
        let mut kernel: HashMap<G, M> = self.maps.iter().map(|(k,v) | {   
           (*k, v.cokernel())
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

    pub fn codomain_space<B: BasisElement>(&self, b: B) -> GradedVectorSpace<G, B> {
        let space: HashMap<G,Vec<B>> = self.maps.iter().map(|(g,m)| {
            (*g, vec![b.clone(); m.codomain()])
        }).collect();
        GradedVectorSpace(space)
    }
}