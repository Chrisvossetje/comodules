use std::sync::Arc;

use ahash::HashMap;
use algebra::{abelian::Abelian, field::Field, matrix::Matrix};
use deepsize::DeepSizeOf;

use crate::{
    basiselement::kBasisElement,
    comodule::traits::CofreeComodule,
    graded_space::{BasisIndex, GradedLinearMap, GradedVectorSpace},
    grading::Grading,
    helper::hashmap_add_restrict,
    tensor::TensorMap,
};

use super::{kcoalgebra::kCoalgebra, kmorphism::kComoduleMorphism, traits::Comodule};

#[derive(Clone, PartialEq, DeepSizeOf)]
#[allow(non_camel_case_types)]
pub struct kComodule<G: Grading, F: Field, M: Matrix<F>> {
    pub coalgebra: Arc<kCoalgebra<G, F, M>>,
    pub space: GradedVectorSpace<G, ()>,
    pub coaction: GradedLinearMap<G, F, M>,
    pub tensor: TensorMap<G>,
}

#[derive(Debug, Clone, PartialEq, DeepSizeOf)]
#[allow(non_camel_case_types)]
pub struct kCofreeComodule<G: Grading, F: Field, M: Matrix<F>> {
    pub coalgebra: Arc<kCoalgebra<G, F, M>>,
    pub space: GradedVectorSpace<G, ((BasisIndex<G>, u16), ())>,
    pub gen_id_gr: Vec<G>, // TODO : Unnecessaery ?
}

pub type Original<G> = (BasisIndex<G>, u16);

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

            let val = self.coaction.maps.get(&t_gr).unwrap().get(m_id as usize, t_id as usize);
            if val != F::one() {
                return false;
            };
        }
        true
    }

    pub fn new(
        coalgebra: Arc<kCoalgebra<G, F, M>>,
        space: GradedVectorSpace<G, ()>,
        coaction: GradedLinearMap<G, F, M>,
        tensor: TensorMap<G>,
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

    // // TODO : Remove
    // pub fn find_cogens(&self, limit: G) -> usize {
    //     let mut temp_coac = self.coaction.clone();

    //     self.space.0.iter().for_each(|(g, els)| {
    //         (0..els.len()).into_iter().for_each(|domain| {
    //             let (_, codomain) = self.tensor.construct[&(*g, domain)][&(G::zero(), 0)];
    //             temp_coac
    //                 .maps
    //                 .get_mut(&g)
    //                 .unwrap()
    //                 .set(domain, codomain, F::zero());
    //         })
    //     });

    //     temp_coac
    //         .maps
    //         .iter()
    //         .filter(|(&gr, _)| gr <= limit)
    //         .map(|(_, map)| {
    //             let kernel = map.kernel();
    //             kernel.codomain()
    //         })
    //         .sum()
    // }
}

impl<G: Grading, F: Field, M: Abelian<F>> CofreeComodule<G> for kCofreeComodule<G, F, M> {
    type Coalgebra = kCoalgebra<G, F, M>;
    type Generator = ();

    fn get_generators(&self) -> Vec<(usize, G, Option<String>)> {
        self.space
            .0
            .iter()
            .flat_map(|(k, v)| {
                v.iter().filter_map(|b| {
                    if b.0.0.0 == G::zero() {
                        Some((b.0.1 as usize, *k, None))
                    } else {
                        None
                    }
                })
            })
            .collect()
    }

    fn direct_sum(&mut self, other: &mut Self) {
        other.space.0.iter_mut().for_each(|(g, other_els)| {
            self.space
                .0
                .entry(*g)
                .and_modify(|self_els| {
                    self_els.extend(other_els.drain(0..));
                })
                .or_insert(other_els.drain(0..).collect());
        });

        self.gen_id_gr.append(&mut other.gen_id_gr);
    }

    fn cofree_comodule(
        coalgebra: Arc<Self::Coalgebra>,
        index: usize,
        grade: G,
        limit: G,
        generator: Self::Generator,
    ) -> Self {
        let space: HashMap<G, Vec<_>> = coalgebra
            .space
            .0
            .iter()
            .filter_map(|(g, v)| {
                let sum = *g + grade;
                if sum <= limit {
                    let k_basis: Vec<_> = (0..v.len())
                        .map(|j| (((*g, j as u32), index as u16), generator))
                        .collect();
                    Some((*g + grade, k_basis))
                } else {
                    None
                }
            })
            .collect();

        // // TODO! SOME UNSTABLE STUFF
        // for ((t_gr, _), (a, m)) in &tensor.deconstruct {
        //     if *t_gr > grade + grade {
        //         tensor.construct.get_mut(m).unwrap().remove(a);
        //     }
        // }
        // tensor.construct.retain(|_, y| {
        //     y.len() > 0
        // });

        // let keys: Vec<G> = space.keys().map(|g| *g).collect();
        // for g in keys {
        //     if g > grade + grade {
        //         space.remove(&g);
        //         coaction.remove(&g);
        //         for t_id in 0..tensor.get_dimension(&g) {
        //             tensor.deconstruct.remove(&(g, t_id));
        //         }
        //         tensor.dimensions.remove(&g);
        //     }
        // }

        // for ((t_gr,t_id), ((a_gr, _), (m_gr, _))) in &tensor.deconstruct {
        //     if a_gr > m_gr {
        //         coaction.get_mut(t_gr).unwrap().set_row_zero(*t_id);
        //     }
        // }

        Self {
            coalgebra,
            space: GradedVectorSpace::from(space),
            gen_id_gr: vec![grade],
        }
    }

    fn zero_comodule(coalgebra: Arc<Self::Coalgebra>) -> Self {
        Self {
            coalgebra: coalgebra,
            space: GradedVectorSpace::new(),
            gen_id_gr: vec![],
        }
    }
}

impl<G: Grading, F: Field, M: Abelian<F>> Comodule<G, F> for kComodule<G, F, M> {
    type Coalgebra = kCoalgebra<G, F, M>;
    // type Morphism = kComoduleMorphism<G, F, M>;
    // type Generator = ();
    // type BaseRing = F;

    fn fp_comodule(coalgebra: Arc<Self::Coalgebra>, degree: G) -> Self {
        let zero = G::zero();

        let space_map: HashMap<G, Vec<_>> = [(degree, vec![()])].into_iter().collect();
        let space = GradedVectorSpace::from(space_map);

        let coact_map: HashMap<G, M> = [(degree, M::identity(1))].into_iter().collect();
        let coaction = GradedLinearMap::from(coact_map);

        // CONNECTED ASSUMPTION
        debug_assert_eq!(
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
        dimensions.insert(degree, 1);

        let mut construct = HashMap::default();
        let mut first_entry = HashMap::default();
        first_entry.insert((zero, 0), (degree, 0));
        construct.insert((degree, 0), first_entry);

        let mut deconstruct = HashMap::default();
        deconstruct.insert((degree, 0), ((zero, 0), (degree, 0)));

        let tensor = TensorMap {
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
}
