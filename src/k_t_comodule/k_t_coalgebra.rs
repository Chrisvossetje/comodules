use std::marker::PhantomData;

use ahash::HashMap;
use algebra::{abelian::Abelian, field::Field, matrices::flat_matrix::FlatMatrix, matrix::Matrix, ring::CRing, rings::{finite_fields::F2, univariate_polynomial_ring::UniPolRing}};
use deepsize::DeepSizeOf;
use itertools::Itertools;
use serde::{Deserialize, Serialize};

use crate::{grading::{grading::Grading, tensor::TensorMap, unigrading::UniGrading}, k_comodule::{graded_space::GradedVectorSpace, kcoalgebra::kCoalgebra}, k_t_comodule::{graded_module_morphism::GradedktFieldMap, k_t_comodule::{ktCofreeComodule, ktComodule}, k_t_morphism::ktComoduleMorphism}, traits::Coalgebra, types::{AGeneratorIndex, BasisElement, CoalgebraIndex, CoalgebraIndexType}};

#[derive(Debug, Clone, PartialEq, Deserialize, Serialize, DeepSizeOf)]
#[allow(non_camel_case_types)]
pub struct ktCoalgebra<G: Grading, F: Field, M: Abelian<UniPolRing<F>>> {
    pub space: GradedVectorSpace<G, (BasisElement, M::Generator)>,
    pub coaction: HashMap<CoalgebraIndex<G>, Vec<(CoalgebraIndex<G>, CoalgebraIndex<G>, UniPolRing<F>)>>,
}


impl <G: Grading, F: Field, M: Abelian<UniPolRing<F>>> Coalgebra<G> for ktCoalgebra<G,F,M> {
    type BaseField = F;
    type BaseRing = UniPolRing<F>;
    type RingMorph = M;
    type Comod = ktComodule<G, F, M>;
    type CofMod = ktCofreeComodule<G, Self>;
    type ComodMorph = ktComoduleMorphism<G, F, M>;

    fn size_in_degree(&self, g: G) -> usize {
        self.space.dimension_in_grade(&g)
    }

    fn coaction(
        &self,
        i: CoalgebraIndex<G>,
    ) -> &[(CoalgebraIndex<G>, CoalgebraIndex<G>, Self::BaseRing)] {
        self.coaction.get(&i).unwrap().as_slice()
    }


    fn basering_comodule(&self, shift: G) -> Self::Comod {
        let zero = G::zero();
        
        let space_map: HashMap<G, Vec<_>> = [(shift, vec![M::Generator::default()])].into_iter().collect();
        let space = GradedVectorSpace::from(space_map);

        let coact_map: HashMap<G, M> = [(shift, M::identity(1))].into_iter().collect();
        let coaction = GradedktFieldMap { maps: coact_map, _p: PhantomData };

        let mut dimensions = HashMap::default();
        dimensions.insert(shift, 1);

        let mut construct = HashMap::default();
        let mut first_entry = HashMap::default();
        first_entry.insert((zero, 0), (shift, 0));
        construct.insert((shift, 0), first_entry);

        let mut deconstruct = HashMap::default();
        deconstruct.insert((shift, 0), ((zero, 0), (shift, 0)));

        let tensor = TensorMap {
            construct,
            deconstruct,
            dimensions,
        };

        Self::Comod {
            space,
            coaction,
            tensor,
        }
    }

    fn cofree_comodule(&self, index: usize, shift: G, limit: G, generator: <Self::RingMorph as Abelian<Self::BaseRing>>::Generator) -> Self::CofMod {
        let space: HashMap<G, Vec<_>> = self
            .space
            .0
            .iter()
            .filter_map(|(g, v)| {
                let sum = *g + shift;
                if sum <= limit {
                    let k_basis: Vec<_> = (0..v.len())
                        .map(|j| (((*g, j as CoalgebraIndexType), index as AGeneratorIndex), generator))
                        .collect();
                    Some((*g + shift, k_basis))
                } else {
                    None
                }
            })
            .collect();

        Self::CofMod {
            space: GradedVectorSpace(space),
        }
    }    

    // fn cofree_comodule(coalgebra: Arc<Self::Coalgebra>, index: usize, grade: G, limit: G, generator: Self::Generator) -> Self {
    //     let coact_maps = match generator.1 {
    //         Some(power) => {
    //             hashmap_add_restrict_transform(&coalgebra.coaction.maps, grade, limit, |mut map| {
    //             for x in 0..map.domain() {
    //                 for y in 0..map.codomain() {
    //                     let el = map.get(x, y);
    //                     if el.1 >= power { // WTF ??
    //                         map.set(x, y, UniPolRing::zero());
    //                     }
    //                 }
    //             }
    //             map
    //         })
    //         },
    //         None => {
    //             hashmap_add_restrict(&coalgebra.coaction.maps, grade, limit)
    //         },
    //     };
            
    //     let space = hashmap_add_restrict_transform(&coalgebra.space.0, grade, limit, |v| {
    //         v.iter()
    //             .map(|basis| {
    //                 let mut el = basis.clone();
    //                 el.0.generated_index = index;
    //                 el.1 = el.1 + generator.0;
    //                 el.2 = generator.1;
    //                 el
    //             })
    //             .collect()
    //     });

    //     let tensor = coalgebra.tensor.add_and_restrict(grade, limit);

    //     let module = ktComodule {
    //         coalgebra,
    //         space: GradedModule(space),
    //         coaction: GradedModuleMap {
    //             maps: coact_maps,
    //         },
    //         tensor,
    //     };

    //     if cfg!(debug_assertions) {
    //         module.verify().unwrap();
    //     }

    //     module
    // }
}

#[allow(non_snake_case)]
pub fn A1_C() -> ktCoalgebra<UniGrading, F2, FlatMatrix<UniPolRing<F2>>> {
    let input = include_str!("../../examples/direct/A(1)_C.txt");
    let coalg = ktCoalgebra::parse(input, UniGrading::infty()).unwrap().0;

    coalg
}

#[allow(non_snake_case)]
pub fn A0_C() -> ktCoalgebra<UniGrading, F2, FlatMatrix<UniPolRing<F2>>> {
    let mut space = HashMap::default();
    space.insert(
        UniGrading(0),
        vec![(
            BasisElement {
                name: "1".to_owned(),
                excess: 0,
            },
            (0,
            None)
        )],
    );

    space.insert(
        UniGrading(1),
        vec![(
            BasisElement {
                name: "t0".to_owned(),
                excess: 1,
            },
            (0,
            None)
        )],
    );

    let mut coaction = HashMap::default();
    coaction.insert((UniGrading::zero(), 0), vec![((UniGrading::zero(), 0),(UniGrading::zero(), 0), UniPolRing::one())]);
    coaction.insert((UniGrading(1), 0), vec![
        ((UniGrading::zero(), 0),(UniGrading(1), 0), UniPolRing::one()),
        ((UniGrading(1), 0),(UniGrading::zero(), 0), UniPolRing::one()) 
        ]);

    ktCoalgebra {
        space: GradedVectorSpace::from(space),
        coaction,
    }
}

pub fn tensor_k_coalgebra(
    coalgebra: kCoalgebra<UniGrading, F2, FlatMatrix<F2>>,
) -> ktCoalgebra<UniGrading, F2, FlatMatrix<UniPolRing<F2>>> {

    let (space, coaction) = (coalgebra.space, coalgebra.coaction);

    let space: HashMap<_, _> = space
        .0
        .into_iter()
        .map(|x| {
            let gr = UniGrading(x.0 .0);
            let module: Vec<(BasisElement, (i32, Option<u16>))> =
                x.1.into_iter().map(|y| (y, (0, None))).collect();
            (gr, module)
        })
        .collect();

    let coaction: HashMap<(UniGrading, u16), Vec<((UniGrading, u16), (UniGrading, u16), UniPolRing<F2>)>> = coaction.into_iter().map(|(gr, m)| {        
        let new_m: Vec<_> = m.into_iter().map(|x| (x.0, x.1, UniPolRing(x.2,0))).collect();
        (gr, new_m)
    }).collect();

    ktCoalgebra {
        space: GradedVectorSpace::from(space),
        coaction,
    }
}
