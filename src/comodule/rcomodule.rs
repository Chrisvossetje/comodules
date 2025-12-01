use std::{sync::Arc};
use ahash::HashMap;
use serde::{Deserialize, Serialize};

use crate::{
    basiselement::kBasisElement, comodule::{
        rmorphism::RComoduleMorphism, traits::Comodule,
    }, grading::{Grading, UniGrading}, helper::{hashmap_add_restrict, hashmap_add_restrict_transform}, linalg::{
        field::Field,
        flat_matrix::FlatMatrix,
        matrix::RModMorphism,
        ring::{CRing, UniPolRing},
    }, module::{module::GradedModule, morphism::GradedModuleMap}, tensor::Tensor
};

#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
pub struct RCoalgebra<G: Grading, F: Field> {
    pub space: GradedModule<G, kBasisElement>,
    pub coaction: GradedModuleMap<G, F>,
    pub tensor: Tensor<G>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct RComodule<G: Grading, F: Field> {
    pub coalgebra: Arc<RCoalgebra<G, F>>,
    pub space: GradedModule<G, kBasisElement>,
    pub coaction: GradedModuleMap<G, F>,
    pub tensor: Tensor<G>,
}


impl<G: Grading, F: Field> RComodule<G, F> {
    pub fn verify(&self) -> Result<(), String> {
        for (g, _) in &self.space.0 {
            if !self.coaction.maps.contains_key(&g) {return Err("Coaction does not contain grade of the space".to_owned()) };
        }

        for (g,m) in &self.coaction.maps {
            if m.domain != self.space.dimension_in_grade(&g) { 
                return Err("Coaction domain does not have correct size".to_owned()) 
            };
            if m.codomain != self.tensor.get_dimension(&g) { 
                return Err("Coaction explanation does not have correct size".to_owned()) 
            };
        }
         
        let tensor_module = self.space.generate_tensor_as_module(&self.coalgebra.as_ref().space, &self.tensor);
        self.coaction.verify(&self.space, &tensor_module)?;

        for (&(_m_gr, m_id), map) in &self.tensor.construct {
            let mut flag = false; 
            
            for (&(a_gr, a_id), &(t_gr, t_id)) in map {
                if !(a_gr == G::zero() && a_id == 0) {
                    continue;
                }

                let val = self.coaction.maps.get(&t_gr).unwrap().get(m_id, t_id);
                if val != UniPolRing::one() {
                    return Err(format!("Value of m -> 1 otimes m is not 1, for basis element {}, {}", _m_gr, m_id ));
                } else {
                    if flag {
                        return Err("There are two elements 1 times m which m maps to".to_owned());
                    }
                    flag = true;
                    break;
                }
            }
            
            if !flag { 
                return Err("m does not map to 1 otimes m".to_owned()) 
            }
        }
        
        Ok(())
    }
}

impl<G: Grading, F: Field> Comodule<G> for RComodule<G, F> {
    type Element = kBasisElement;
    type Coalgebra = RCoalgebra<G, F>;
    type Morphism = RComoduleMorphism<G, F>;
    type BaseRing = UniPolRing<F>;
    
    type Generator = (UniGrading, Option<u16>);

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
                        let s = format!("{:?} | {:?}", b.2, b.1);
                        Some((b.0.generated_index, *k, Some(s)))
                    } else {
                        None
                    }
                })
            })
            .collect()
    }

    fn fp_comodule(coalgebra: Arc<Self::Coalgebra>, degree: G) -> Self {
        let zero = G::zero();

        let el = kBasisElement {
            name: "fp".to_string(),
            generator: false,
            primitive: None,
            generated_index: 0,
        };

        let space_map =
            [(degree, vec![(el, UniGrading(0), None)])].into_iter().collect();
        let space = GradedModule(space_map);

        let coact_map: HashMap<G, FlatMatrix<UniPolRing<F>>> =
            [(degree, FlatMatrix::identity(1))].into_iter().collect();

        let coaction = GradedModuleMap {
            maps: coact_map,
        };

        let mut dimensions = HashMap::default();
        dimensions.insert(degree, 1);

        let mut construct = HashMap::default();
        let mut first_entry = HashMap::default();
        first_entry.insert((zero, 0), (degree, 0));
        construct.insert((degree, 0), first_entry);

        let mut deconstruct = HashMap::default();
        deconstruct.insert((degree, 0), ((zero, 0), (degree, 0)));

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

        if cfg!(debug_assertions) {
            assert!(self.tensor.is_correct());
            self.verify().unwrap();

        }
    }

    fn cofree_comodule(coalgebra: Arc<Self::Coalgebra>, index: usize, grade: G, limit: G, generator: Self::Generator) -> Self {
        let coact_maps = match generator.1 {
            Some(power) => {
                hashmap_add_restrict_transform(&coalgebra.coaction.maps, grade, limit, |mut map| {
                for x in 0..map.domain {
                    for y in 0..map.codomain {
                        let el = map.get(x, y);
                        if el.1 >= power { // WTF ??
                            map.set(x, y, UniPolRing::zero());
                        }
                    }
                }
                map
            })
            },
            None => {
                hashmap_add_restrict(&coalgebra.coaction.maps, grade, limit)
            },
        };
            
        let space = hashmap_add_restrict_transform(&coalgebra.space.0, grade, limit, |v| {
            v.iter()
                .map(|basis| {
                    let mut el = basis.clone();
                    el.0.generated_index = index;
                    el.1 = el.1 + generator.0;
                    el.2 = generator.1;
                    el
                })
                .collect()
        });

        let tensor = coalgebra.tensor.add_and_restrict(grade, limit);

        let module = RComodule {
            coalgebra,
            space: GradedModule(space),
            coaction: GradedModuleMap {
                maps: coact_maps,
            },
            tensor,
        };

        if cfg!(debug_assertions) {
            module.verify().unwrap();
        }

        module
    }
}
