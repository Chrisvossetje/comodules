use std::sync::Arc;

use ahash::HashMap;
use serde::{Deserialize, Serialize};

use crate::{
    comodule::{
        kcomodule::kBasisElement, rmorphism::RComoduleMorphism, tensor::Tensor, traits::Comodule,
    }, helper::{hashmap_add_restrict, hashmap_add_restrict_transform}, linalg::{
        field::Field,
        flat_matrix::FlatMatrix,
        matrix::RModMorphism,
        module::{GradedModule, GradedModuleMap, PolyGrading},
        ring::UniPolRing,
    }
};

#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
pub struct RCoalgebra<G: PolyGrading, F: Field> {
    pub space: GradedModule<G, kBasisElement>,
    pub coaction: GradedModuleMap<G, F>,
    pub tensor: Tensor<G>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct RComodule<G: PolyGrading, F: Field> {
    pub coalgebra: Arc<RCoalgebra<G, F>>,
    pub space: GradedModule<G, kBasisElement>,
    pub coaction: GradedModuleMap<G, F>,
    pub tensor: Tensor<G>,
}

impl<G: PolyGrading, F: Field> Comodule<G> for RComodule<G, F> {
    type Element = kBasisElement;
    type Coalgebra = RCoalgebra<G, F>;
    type Morphism = RComoduleMorphism<G, F>;

    fn zero_comodule(coalgebra: Arc<Self::Coalgebra>) -> Self {
        Self {
            coalgebra: coalgebra,
            space: GradedModule::default(),
            coaction: GradedModuleMap::default(),
            tensor: Tensor::default(),
        }
    }

    fn get_generators(&self) -> Vec<(usize, G, Option<String>)> {
        self.space
            .0
            .iter()
            .flat_map(|(k, v)| {
                v.iter().filter_map(|b| {
                    if b.0.generator {
                        let s = format!("{:?}", b.1);
                        Some((b.0.generated_index, *k, Some(s)))
                    } else {
                        None
                    }
                })
            })
            .collect()
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

        let space_map: HashMap<G, Vec<(kBasisElement, Option<usize>)>> =
            [(zero, vec![(el, None)])].into_iter().collect();
        let space = GradedModule(space_map);

        let coact_map: HashMap<G, FlatMatrix<UniPolRing<F>>> =
            [(zero, FlatMatrix::identity(1))].into_iter().collect();
        let mut domain_explain = HashMap::default();
        domain_explain.insert(zero, vec![((zero, 0), 0)]);
        let mut codomain_explain = HashMap::default();
        codomain_explain.insert(zero, vec![((zero, 0), 0)]);
        let coaction = GradedModuleMap {
            maps_domain: domain_explain,
            maps_codomain: codomain_explain,
            maps: coact_map,
        };

        let mut dimensions = HashMap::default();
        dimensions.insert(zero, 1);

        let mut construct = HashMap::default();
        let mut first_entry = HashMap::default();
        first_entry.insert((zero, 0), (zero, 0));
        construct.insert((zero, 0), first_entry);

        let mut deconstruct = HashMap::default();
        deconstruct.insert((zero, 0), ((zero, 0), (zero, 0)));

        let tensor = Tensor {
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

        // todo!();
        // debug_assert!(self.verify());
    }

    fn cofree_comodule(coalgebra: Arc<Self::Coalgebra>, index: usize, grade: G, limit: G) -> Self {
        let coact_maps = hashmap_add_restrict(&coalgebra.coaction.maps, grade, limit);
        let coact_domain =
            hashmap_add_restrict_transform(&coalgebra.coaction.maps_domain, grade, limit, |x| {
                x.into_iter()
                    .map(|((g, id), power)| ((g + grade, id), power))
                    .collect()
            });
        let coact_codomain =
            hashmap_add_restrict_transform(&coalgebra.coaction.maps_domain, grade, limit, |x| {
                x.into_iter()
                    .map(|((g, id), power)| ((g + grade, id), power))
                    .collect()
            });

        let space = hashmap_add_restrict_transform(&coalgebra.space.0, grade, limit, |v| {
            v.iter()
                .map(|basis| {
                    let mut el = basis.clone();
                    el.0.generated_index = index;
                    el
                })
                .collect()
        });
        let tensor = coalgebra.tensor.add_and_restrict(grade, limit);

        RComodule {
            coalgebra,
            space: GradedModule(space),
            coaction: GradedModuleMap {
                maps_domain: coact_domain,
                maps_codomain: coact_codomain,
                maps: coact_maps,
            },
            tensor,
        }
    }
}
