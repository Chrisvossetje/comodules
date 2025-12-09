use std::sync::Arc;

use ahash::HashMap;
use algebra::{abelian::Abelian, field::Field, matrix::Matrix};
use serde::de;

use crate::{
    basiselement::kBasisElement, graded_space::{BasisIndex, GradedLinearMap, GradedVectorSpace}, grading::{Grading, UniGrading}, helper::hashmap_add_restrict, tensor::{TensorList, TensorMap, tensor_list_find_tensor_id}
};

use super::{
    kcoalgebra::kCoalgebra, kmorphism::kComoduleMorphism, traits::Comodule,
};



#[derive(Clone, PartialEq)]
#[allow(non_camel_case_types)]
pub struct kComodule<G: Grading, F: Field, M: Matrix<F>> {
    pub coalgebra: Arc<kCoalgebra<G, F, M>>,
    pub space: GradedVectorSpace<G, kBasisElement>,
    pub coaction: GradedLinearMap<G, F, M>,
    pub tensor: TensorList<G>,
}

impl<G: Grading, F: Field, M: Matrix<F>> std::fmt::Debug for kComodule<G, F, M> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self.space.0)
    }
}

impl<G: Grading, F: Field, M: Matrix<F>> kComodule<G, F, M> {
    pub fn verify(&self) -> bool {
        for (&m_gr, m_els) in &self.space.0 {
            for m_id in 0..m_els.len() {
                let t_id = tensor_list_find_tensor_id(&self.tensor, m_gr, ((G::zero(), 0), (m_gr, m_id)));
    
                let val = self.coaction.maps.get(&m_gr).unwrap().get(m_id, t_id);
                if val != F::one() {
                    return false;
                };
            } 
        }
        true
    }

    pub fn new(
        coalgebra: Arc<kCoalgebra<G, F, M>>,
        space: GradedVectorSpace<G, kBasisElement>,
        coaction: GradedLinearMap<G, F, M>,
        tensor: TensorList<G>,
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

impl<G: Grading, F: Field, M: Abelian<F>> Comodule<G> for kComodule<G, F, M> {
    type Element = kBasisElement;
    type Coalgebra = kCoalgebra<G, F, M>;
    type Morphism = kComoduleMorphism<G, F, M>;
    type Generator = ();
    type BaseRing = F;

    fn get_generators(&self) -> Vec<(usize, G, Option<String>)> {
        self.space
            .0
            .iter()
            .flat_map(|(k, v)| {
                v.iter().filter_map(|b| {
                    if b.generator {
                        Some((b.generated_index, *k, None))
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
            tensor: TensorList::default(),
        }
    }

    fn fp_comodule(coalgebra: Arc<Self::Coalgebra>, degree: G) -> Self {
        let zero = G::zero();

        let el = kBasisElement {
            name: "fp".to_string(),
            generator: false,
            primitive: None,
            generated_index: 0,
        };

        let space_map: HashMap<G, Vec<kBasisElement>> = [(degree, vec![el])].into_iter().collect();
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

        let mut tensor = HashMap::default();
        tensor.insert(degree, vec![((G::zero(), 0), (degree, 0))]);

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


        fn tensor_list_direct_sum<G: Grading>(left: &mut TensorList<G>, right: &mut TensorList<G>, dimensions: &HashMap<G,usize>) {
            for (gr, right_els) in right {
                let left_els = left.entry(*gr).or_default();
                for (a, (m_gr, m_id)) in right_els.drain(..) {
                    let module_size = *dimensions.get(&m_gr).unwrap_or(&0);
                    left_els.push((a,(m_gr, m_id + module_size)));
                }
            }
        }

        tensor_list_direct_sum(&mut self.tensor, &mut other.tensor, &self_dimensions);

        debug_assert!(self.verify());
    }

    fn cofree_comodule(coalgebra: Arc<Self::Coalgebra>, index: usize, grade: G, limit: G, _generator: Self::Generator) -> Self {
        let coaction: HashMap<G, M> = hashmap_add_restrict(&coalgebra.coaction.maps, grade, limit);
        let space: HashMap<G, Vec<kBasisElement>> = coalgebra
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

        let tensor: HashMap<_,_> = coalgebra.tensor.iter().filter_map(|(g,els)| {
            let sum = *g + grade;
            if sum <= limit {
                let v: Vec<_> = els.iter().map(|(a,(m_gr, m_id))| (*a, (*m_gr + grade, *m_id))).collect();
                Some((*g + grade,v))
            } else {
                None
            }
        }).collect();
        
        
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
        


        kComodule::new(
            coalgebra,
            GradedVectorSpace(space),
            GradedLinearMap::from(coaction),
            tensor,
        )
    }    
}
