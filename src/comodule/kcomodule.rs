use std::{collections::HashMap, sync::Arc};

use ahash::RandomState;
use serde::{Deserialize, Serialize};

use crate::linalg::{
    field::Field,
    graded::{BasisElement, GradedLinearMap, GradedVectorSpace},
    grading::Grading,
    matrix::Matrix,
};

use super::{kcoalgebra::kCoalgebra, ktensor::kTensor, traits::Comodule};

#[derive(Debug, Clone, PartialEq, Default, Deserialize, Serialize)]
#[allow(non_camel_case_types)]
pub struct kBasisElement {
    pub name: String,
    pub generator: bool,
    pub primitive: Option<usize>,
    pub generated_index: usize,
}

impl BasisElement for kBasisElement {}

#[derive(Clone, PartialEq)]
#[allow(non_camel_case_types)]
pub struct kComodule<G: Grading, F: Field, M: Matrix<F>> {
    pub coalgebra: Arc<kCoalgebra<G, F, M>>,
    pub space: GradedVectorSpace<G, kBasisElement>,
    pub coaction: GradedLinearMap<G, F, M>,
    pub tensor: kTensor<G>,
}

impl<G: Grading, F: Field, M: Matrix<F>> std::fmt::Debug for kComodule<G, F, M> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self.space.0)
    }
}

impl<G: Grading, F: Field, M: Matrix<F>> kComodule<G, F, M> {
    pub fn verify(&self) -> bool {
        for (&(m_gr, m_id), map) in &self.tensor.construct {
            let &(t_gr, t_id) = map.get(&(G::zero(), 0)).unwrap();
            if t_gr != m_gr {
                return false;
            };

            let val = self.coaction.maps.get(&t_gr).unwrap().get(m_id, t_id);
            if val != F::one() {
                return false;
            };
        }
        true
    }

    pub fn new(
        coalgebra: Arc<kCoalgebra<G, F, M>>,
        space: GradedVectorSpace<G, kBasisElement>,
        coaction: GradedLinearMap<G, F, M>,
        tensor: kTensor<G>,
    ) -> Self {
        let com = Self {
            coalgebra,
            space,
            coaction,
            tensor,
        };
        debug_assert!(com.verify());
        com
    }

    pub fn find_cogens(&self, limit: G) -> usize {
        let mut temp_coac = self.coaction.clone();

        self.space.0.iter().for_each(|(g,els)| {
            (0..els.len()).into_iter().for_each(|domain| {
                let (_, codomain) = self.tensor.construct[&(*g,domain)][&(G::zero(), 0)];
                temp_coac.maps.get_mut(&g).unwrap().set(domain, codomain, F::zero());
            })
        });

        temp_coac.maps.iter().filter(|(&gr,_)| {gr <= limit}).map(|(gr, map)| {
            let kernel = map.kernel();
            kernel.codomain()
        }).sum()
    }
}

impl<G: Grading, F: Field, M: Matrix<F>> Comodule<G> for kComodule<G, F, M> {
    type Element = kBasisElement;
    type Coalgebra = kCoalgebra<G, F, M>;

    fn get_generators(&self) -> Vec<(usize, G, Option<String>)> {
        self.space
            .0
            .iter()
            .flat_map(|(k, v)| {
                v.iter().filter_map(|b| {
                    if b.generator {
                        Some((b.generated_index, *k, Some(b.name.clone())))
                    } else {
                        None
                    }
                })
            })
            .collect()
    }

    fn zero_comodule(coalgebra: Arc<Self::Coalgebra>) -> Self {
        Self {
            coalgebra: coalgebra,
            space: GradedVectorSpace::new(),
            coaction: GradedLinearMap::empty(),
            tensor: kTensor::new(),
        }
    }

    fn fp_comodule(coalgebra: Arc<Self::Coalgebra>) -> Self {
        let zero = G::zero();

        // (zero,
        let el = kBasisElement {
            name: "fp".to_string(),
            generator: false,
            primitive: None,
            generated_index: 0,
        };

        let space_map: HashMap<G, Vec<kBasisElement>, RandomState> =
            [(zero, vec![el])].into_iter().collect();
        let space = GradedVectorSpace::from(space_map);

        let coact_map: HashMap<G, M, RandomState> = [(zero, M::identity(1))].into_iter().collect();
        let coaction = GradedLinearMap::from(coact_map);

        assert_eq!(
            coalgebra
                .space
                .0
                .get(&zero)
                .expect("Coalgebra has no element in grade zero")
                .len(),
            1,
            "Coalgebra is not a connected coalgebra"
        );

        let mut dimensions = HashMap::default();
        dimensions.insert(zero, 1);

        let mut construct = HashMap::default();
        let mut first_entry = HashMap::default();
        first_entry.insert((zero, 0), (zero, 0));
        construct.insert((zero, 0), first_entry);

        let mut deconstruct = HashMap::default();
        deconstruct.insert((zero, 0), ((zero, 0), (zero, 0)));

        let tensor = kTensor {
            construct,
            deconstruct,
            dimensions,
        };

        Self {
            coalgebra,
            space,
            coaction,
            tensor,
        }
    }

    fn direct_sum(&mut self, other: &mut Self) {
        self.coaction.block_sum(&mut other.coaction);

        let self_dimensions = self.space.0.iter().map(|(g, v)| (*g, v.len())).collect();

        other.space.0.iter_mut().for_each(|(g, other_els)| {
            self.space
                .0
                .entry(*g)
                .and_modify(|self_els| {
                    self_els.extend(other_els.drain(0..));
                })
                .or_insert(other_els.drain(0..).collect());
        });

        self.tensor.direct_sum(&mut other.tensor, &self_dimensions);

        debug_assert!(self.verify());
    }

    fn cofree_comodule(coalgebra: Arc<Self::Coalgebra>, index: usize, grade: G, limit: G) -> Self {
        let coaction: HashMap<G, M, RandomState> = coalgebra
            .coaction
            .maps
            .iter()
            .filter_map(|(g, v)| {
                let sum = *g + grade;
                if sum <= limit {
                    Some((*g + grade, v.clone()))
                } else {
                    None
                }
            })
            .collect();
        let space = coalgebra
            .space
            .0
            .iter()
            .filter_map(|(g, v)| {
                let sum = *g + grade;
                if sum <= limit {
                    let k_basis: Vec<kBasisElement> = v
                        .iter()
                        .map(|basis| {
                            let mut el = basis.clone();
                            el.generated_index = index;
                            el
                        })
                        .collect();
                    Some((*g + grade, k_basis))
                } else {
                    None
                }
            })
            .collect();
        let tensor = coalgebra.tensor.add_and_restrict(grade, limit);

        kComodule::new(
            coalgebra,
            GradedVectorSpace(space),
            GradedLinearMap::from(coaction),
            tensor,
        )
    }
}
