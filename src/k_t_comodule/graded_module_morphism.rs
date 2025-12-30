use std::marker::PhantomData;

use ahash::HashMap;
use algebra::{
    abelian::Abelian, field::Field, rings::univariate_polynomial_ring::UniPolRing
};
use deepsize::DeepSizeOf;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

use crate::{grading::{grading::Grading}, k_comodule::graded_space::GradedVectorSpace, k_t_comodule::k_t_coalgebra::ktCoalgebra, traits::Coalgebra};


#[derive(Debug, Clone, PartialEq, Default, DeepSizeOf)]
pub struct GradedktFieldMap<G: Grading, F: Field, M: Abelian<UniPolRing<F>>> {
    pub maps: HashMap<G, M>,
    pub(crate) _p: PhantomData<F>,
}

impl<G: Grading, F: Field, M: Abelian<UniPolRing<F>>> GradedktFieldMap<G, F, M> {
    pub fn zero_codomain(domain: &GradedVectorSpace<G, M::Generator>) -> Self {
        let maps = domain
            .0
            .iter()
            .map(|(g, v)| (*g, M::zero(v.len(), 0)))
            .collect();

        Self { maps, _p: PhantomData }
    }

    pub fn vstack(&mut self, other: &mut Self) {
        self.maps
            .iter_mut()
            .for_each(|(grade, self_mat)| match other.maps.get_mut(grade) {
                Some(other_mat) => {
                    self_mat.vstack(other_mat);
                }
                None => {}
            });
        other.maps.drain().for_each(|(g, map)| {
            if !self.maps.contains_key(&g) {
                self.maps.insert(g, map);
            }
        });
    }

    /// Note that calling this function invalidates the domains / codomains it references unless these are also summed
    pub fn block_sum(&mut self, other: &mut Self) {
        self.maps
            .iter_mut()
            .for_each(|(grade, self_mat)| match other.maps.get_mut(grade) {
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

    pub fn verify(
        &self,
        domain: &GradedVectorSpace<G, M::Generator>,
        codomain: &GradedVectorSpace<G, M::Generator>,
    ) -> Result<(), String> {
        for (grade, map) in &self.maps {
            let empty = vec![];
            let dom: Vec<_> = domain
                .0
                .get(grade)
                .unwrap_or(&empty)
                .iter()
                .collect();
            let codom: Vec<_> = codomain
                .0
                .get(grade)
                .unwrap_or(&empty)
                .iter()
                .collect();
            // TODO
            // map.verify(&dom, &codom)?;
        }
        Ok(())
    }

    // pub fn reduce(&self) -> Self {
    //     let mut total_map = HashMap::default();
    //     let mut total_expl = HashMap::default();

    //     for (g, map) in &self.maps {
    //         let expl = self.explained.get(g).expect("The codomain should have been explained for the map in this grade");

    //         let mut new_codoms = vec![];
    //         for codom in 0..map.codomain {
    //             let zero = map.get_row(codom).iter().fold(true, |b, r| r.is_zero() && b);
    //             if !zero {
    //                 new_codoms.push(codom);
    //             }
    //         }

    //         let mut new_map = FlatMatrix::zero(map.domain, new_codoms.len());
    //         let mut new_expl = vec![];

    //         for (index, codom) in new_codoms.iter().enumerate() {
    //             let row = map.get_row(*codom);
    //             new_map.set_row(index, row);

    //             let exp = expl[*codom];
    //             new_expl.push(exp);
    //         }

    //         total_map.insert(*g, new_map);
    //         total_expl.insert(*g, new_expl);
    //     }

    //     Self {
    //         maps: total_map,
    //         explained: total_expl
    //     }
    // }

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

        // TODO ! check if powers are now too high

        Self { maps: compose, _p: PhantomData}
    }
}

impl<G: Grading, F: Field, M: Abelian<UniPolRing<F>>> GradedktFieldMap<G, F, M> {
    pub fn zero(
        domain: &GradedVectorSpace<G, M::Generator>,
        codomain: &GradedVectorSpace<G, M::Generator>,
    ) -> Self {
        let mut maps: HashMap<G, M> = domain
            .0
            .iter()
            .map(|(g, els)| {
                let codom_len = codomain.dimension_in_grade(g);
                (*g, M::zero(els.len(), codom_len))
            })
            .collect();
        codomain.0.iter().for_each(|(g, v)| {
            if !maps.contains_key(g) {
                maps.insert(*g, M::zero(0, v.len()));
            }
        });
        Self { maps, _p: PhantomData }
    }

    /// If self is a map from A -> B, let Q be the cokernel.
    /// Then we return the map from (B -> Q, Q -> B, Q).
    /// "The map from B to Q, an 'inverse' map from Q to B and the cokernel object Q.
    pub fn cokernel(
        &self,
        codomain: &GradedVectorSpace<G, M::Generator>,
    ) -> (Self, Self, GradedVectorSpace<G, M::Generator>) {
        if cfg!(debug_assertions) {
            for g in codomain.0.keys() {
                assert!(self.maps.contains_key(g));
            }
        }

        let mut coker = HashMap::default();
        let mut coker_map = HashMap::default();
        let mut coker_inv_map = HashMap::default();

        for (&g, mat) in &self.maps {
            // TODO, parallelization
            let codom = &codomain.0[&g];
            let (g_map, g_inv_map, module) = mat.cokernel(codom);

            if module.len() > 0 {
                coker.insert(g, module);
            }

            coker_map.insert(g, g_map);
            coker_inv_map.insert(g, g_inv_map);
        }

        let to = GradedktFieldMap { maps: coker_map, _p: PhantomData };
        let inv = GradedktFieldMap {
            maps: coker_inv_map, _p: PhantomData
        };
        let coker = GradedVectorSpace(coker);

        if cfg!(debug_assertions) {
            to.verify(codomain, &coker).unwrap();
        }

        (to, inv, coker)
    }
}
