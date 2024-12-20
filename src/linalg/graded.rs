use std::{collections::HashMap, fmt::Debug, marker::PhantomData};

use ahash::RandomState;
use rayon::prelude::*;

use super::{field::Field, grading::Grading, matrix::Matrix};
use serde::{Deserialize, Serialize};

pub trait BasisElement: 'static + Debug + Clone {}

pub type BasisIndex<G> = (G, usize);

#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
pub struct GradedVectorSpace<G: Grading, B: BasisElement>(pub HashMap<G, Vec<B>, RandomState>);

#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
pub struct GradedLinearMap<G: Grading, F: Field, M: Matrix<F>> {
    pub maps: HashMap<G, M, RandomState>,
    __: PhantomData<F>,
}

impl<G: Grading, B: BasisElement> GradedVectorSpace<G, B> {
    pub fn new() -> Self {
        Self(HashMap::default())
    }
    pub fn dimension_in_grade(&self, grade: &G) -> usize {
        self.0.get(grade).map(|x| x.len()).unwrap_or(0)
    }
}

impl<G: Grading, B: BasisElement> From<HashMap<G, Vec<B>, RandomState>>
    for GradedVectorSpace<G, B>
{
    fn from(value: HashMap<G, Vec<B>, RandomState>) -> Self {
        Self(value)
    }
}

impl<G: Grading, F: Field, M: Matrix<F>> From<HashMap<G, M, RandomState>>
    for GradedLinearMap<G, F, M>
{
    fn from(value: HashMap<G, M, RandomState>) -> Self {
        Self {
            maps: value,
            __: PhantomData,
        }
    }
}

impl<G: Grading, F: Field, M: Matrix<F>> GradedLinearMap<G, F, M> {
    pub fn get_cokernel(&self) -> Self {
        let cokernel = self
            .maps
            .par_iter()
            .map(|(k, v)| (*k, v.cokernel()))
            .collect();
        GradedLinearMap {
            maps: cokernel,
            __: PhantomData,
        }
    }

    pub fn get_kernel(&self) -> Self {
        let kernel = self
            .maps
            .par_iter()
            .map(|(k, v)| (*k, v.kernel()))
            .collect();
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
        let mut compose: HashMap<G, M, RandomState> = self
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
        self.maps
            .par_iter()
            .map(|(k, v)| (*k, v.pivots()))
            .collect()
    }

    pub fn empty() -> Self {
        GradedLinearMap {
            maps: HashMap::default(),
            __: PhantomData,
        }
    }

    pub fn zero<B: BasisElement>(
        domain: &GradedVectorSpace<G, B>,
        codomain: &GradedVectorSpace<G, B>,
    ) -> Self {
        let mut maps: HashMap<G, M, RandomState> = domain
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
        let space = self
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
