use std::{fmt::Debug, marker::PhantomData};

use ahash::HashMap;
use algebra::{abelian::Abelian, field::Field, matrix::Matrix};
use deepsize::DeepSizeOf;
use rayon::prelude::*;

use serde::{Deserialize, Serialize};

use crate::grading::grading::Grading;

#[derive(Debug, Clone, PartialEq, Deserialize, Serialize, DeepSizeOf)]
pub struct GradedVectorSpace<G: Grading, B>(pub HashMap<G, Vec<B>>);

#[derive(Debug, Clone, PartialEq, Deserialize, Serialize, DeepSizeOf)]
pub struct GradedLinearMap<G: Grading, F: Field, M: Matrix<F>> {
    pub maps: HashMap<G, M>,
    __: PhantomData<F>,
}

impl<G: Grading, B> GradedVectorSpace<G, B> {
    pub fn new() -> Self {
        Self(HashMap::default())
    }
    pub fn dimension_in_grade(&self, grade: &G) -> usize {
        self.0.get(grade).map(|x| x.len()).unwrap_or(0)
    }
}

impl<G: Grading, B> From<HashMap<G, Vec<B>>> for GradedVectorSpace<G, B> {
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

impl<G: Grading, F: Field, M: Abelian<F>> GradedLinearMap<G, F, M> {
    pub fn empty() -> Self {
        GradedLinearMap {
            maps: HashMap::default(),
            __: PhantomData,
        }
    }
}

impl<G: Grading, F: Field, M: Abelian<F>> GradedLinearMap<G, F, M> {
    pub fn get_cokernel(&self) -> Self {
        let cokernel = self
            .maps
            .par_iter()
            .map(|(k, v)| (*k, v.cokernel(&vec![]).0))
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
            .map(|(k, v)| (*k, v.kernel(&vec![], &vec![]).0))
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
        self.maps
            .par_iter_mut()
            .for_each(|(grade, self_mat)| match other.maps.get(grade) {
                Some(other_mat) => {
                    self_mat.block_sum(other_mat);
                }
                None => {}
            });
        other.maps.drain().for_each(|(g, map)| {
            if !self.maps.contains_key(&g) {
                self.maps.insert(g, map);
            }
        });
    }

    pub fn compose(&self, rhs: &Self) -> Self {
        let mut compose: HashMap<G, M> = self
            .maps
            .par_iter()
            .filter_map(|(k, v)| match rhs.maps.get(&k) {
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
            .map(|(k, v)| {
                let mut col = 0;
                let mut pivots = vec![];
                for row in 0..v.codomain() {
                    while col < v.domain() {
                        if !v.get(col, row).is_zero() {
                            pivots.push((col, row));
                            col += 1;
                            break;
                        }
                        col += 1;
                    }
                }
                (*k, pivots)
            })
            .collect()
    }

    pub fn zero_codomain<B>(codomain: &GradedVectorSpace<G, B>) -> Self {
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

    pub fn codomain_space<B: Clone>(&self, b: B) -> GradedVectorSpace<G, B> {
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
