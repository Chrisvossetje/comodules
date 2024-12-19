use std::{
    collections::HashMap,
    fmt::Debug,
    hash::Hash,
    marker::PhantomData,
    ops::{Add, AddAssign, Sub, SubAssign},
    str::FromStr,
};

use super::{field::Field, matrix::Matrix};

pub trait Grading:
    'static
    + Clone
    + Hash
    + Copy
    + Debug
    + Sized
    + Add<Output = Self>
    + Sub<Output = Self>
    + PartialEq
    + Eq
    + AddAssign
    + SubAssign
    + PartialOrd
    + Ord
{
    fn degree_names() -> Vec<char>;
    fn default_formulas() -> (String, String);
    fn export_grade(self) -> Vec<i32>;

    fn incr(self) -> Self;
    fn zero() -> Self;

    fn parse(parse: &str) -> Result<Self, ()>;
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

    fn incr(self) -> Self {
        self + 1
    }

    fn parse(parse: &str) -> Result<Self, ()> {
        i32::from_str(parse).map_err(|_| ())
    }
}

pub type UniGrading = i32;

#[derive(Debug, Clone, Copy, Hash)]
pub struct BiGrading(i32, i32);

pub trait BasisElement: 'static + Debug + Clone {}

pub type BasisIndex<G> = (G, usize);


// A VectorSpace should be naive / simple, just a list of basis elements!
// Specific modules can implement their own "basis" type which encodes the information they need
pub type VectorSpace<B> = Vec<B>;

// Maybe make this its own type ???
// This is probably fine, as modules will always direct use this type
// pub type GradedVectorSpace<G: Grading, B: BasisElement> = HashMap<G, VectorSpace<B>>;

#[derive(Debug, Clone, PartialEq)]
pub struct GradedVectorSpace<G: Grading, B: BasisElement>(pub HashMap<G, VectorSpace<B>>);

#[derive(Debug, Clone, PartialEq)]
pub struct GradedLinearMap<G: Grading, F: Field, M: Matrix<F>> {
    pub maps: HashMap<G, M>,
    __: PhantomData<F>,
}

impl<G: Grading, B: BasisElement> GradedVectorSpace<G, B> {
    pub fn new() -> Self {
        Self(HashMap::new())
    }
    pub fn dimension_in_grade(&self, grade: &G) -> usize {
        self.0.get(grade).map(|x| x.len()).unwrap_or(0)
    }
}

impl<G: Grading, B: BasisElement> From<HashMap<G, Vec<B>>> for GradedVectorSpace<G, B> {
    fn from(value: HashMap<G, Vec<B>>) -> Self {
        Self(value)
    }
}

impl<G: Grading, F: Field, M: Matrix<F>> From<HashMap<G, M>> for GradedLinearMap<G, F, M> {
    fn from(value: HashMap<G, M>) -> Self {
        Self {
            maps: value,
            __: PhantomData,
        }
    }
}

impl<G: Grading, F: Field, M: Matrix<F>> GradedLinearMap<G, F, M> {
    pub fn get_cokernel(&self) -> Self {
        let cokernel: HashMap<G, M> = self.maps.iter().map(|(k, v)| (*k, v.cokernel())).collect();
        GradedLinearMap {
            maps: cokernel,
            __: PhantomData,
        }
    }

    pub fn get_kernel(&self) -> Self {
        let kernel: HashMap<G, M> = self.maps.iter().map(|(k, v)| (*k, v.kernel())).collect();
        GradedLinearMap {
            maps: kernel,
            __: PhantomData,
        }
    }

    pub fn vstack(&mut self, other: &mut Self) {
        other.maps.iter_mut().for_each(|(grade, other_mat)| {
            self.maps
                .entry(*grade)
                .and_modify(|self_mat| {
                    self_mat.vstack(other_mat);
                })
                .or_insert(other_mat.clone());
        });
    }

    pub fn block_sum(&mut self, other: &mut Self) {
        other.maps.iter_mut().for_each(|(grade, other_mat)| {
            self.maps
                .entry(*grade)
                .and_modify(|self_mat| {
                    self_mat.block_sum(other_mat);
                })
                .or_insert(other_mat.clone());
        });
    }

    pub fn compose(self, mut rhs: Self) -> Self {
        let mut compose: HashMap<G, M> = self
            .maps
            .iter()
            .filter_map(|(k, v)| match rhs.maps.get_mut(&k) {
                None => None,
                Some(t) => Some((*k, v.compose(t))),
            })
            .collect();

        for (self_gr, val) in self.maps.iter() {
            if !compose.contains_key(self_gr) {
                compose.insert(*self_gr, M::zero(0, val.codomain()));
            }
        }

        for (rhs_gr, val) in rhs.maps.iter() {
            if !compose.contains_key(rhs_gr) {
                compose.insert(*rhs_gr, M::zero(val.domain(), 0));
            }
        }

        GradedLinearMap {
            maps: compose,
            __: PhantomData,
        }
    }

    pub fn pivots(&self) -> HashMap<G, Vec<(usize, usize)>> {
        self.maps.iter().map(|(k, v)| (*k, v.pivots())).collect()
    }

    pub fn empty() -> Self {
        GradedLinearMap {
            maps: HashMap::new(),
            __: PhantomData,
        }
    }

    pub fn zero<B: BasisElement>(
        domain: &GradedVectorSpace<G, B>,
        codomain: &GradedVectorSpace<G, B>,
    ) -> Self {
        let mut maps: HashMap<G, M> = domain
            .0
            .iter()
            .map(|(g, els)| {
                let codom_len = codomain.dimension_in_grade(g);
                (*g, Matrix::zero(els.len(), codom_len))
            })
            .collect();
        codomain.0.iter().for_each(|(g, v)| {
            if !maps.contains_key(g) {
                maps.insert(*g, Matrix::zero(0, v.len()));
            }
        });
        Self {
            maps,
            __: PhantomData,
        }
    }

    pub fn zero_codomain<B: BasisElement>(codomain: &GradedVectorSpace<G, B>) -> Self {
        let maps = codomain
            .0
            .iter()
            .map(|(g, els)| (*g, Matrix::zero(els.len(), 0)))
            .collect();
        Self {
            maps,
            __: PhantomData,
        }
    }

    pub fn codomain_space<B: BasisElement>(&self, b: B) -> GradedVectorSpace<G, B> {
        let space: HashMap<G, Vec<B>> = self
            .maps
            .iter()
            .filter_map(|(g, m)| match m.codomain() {
                0 => None,
                s => Some((*g, vec![b.clone(); s])),
            })
            .collect();
        GradedVectorSpace(space)
    }
}
