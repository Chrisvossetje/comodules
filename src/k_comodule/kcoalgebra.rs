use std::marker::PhantomData;

use ahash::HashMap;
use algebra::{
    abelian::Abelian, field::Field, matrices::flat_matrix::FlatMatrix, matrix::Matrix, ring::CRing,
    rings::finite_fields::F2,
};
use deepsize::DeepSizeOf;

use serde::{Deserialize, Serialize};

use crate::{
    grading::{grading::Grading, tensor::TensorMap, unigrading::UniGrading},
    k_comodule::{
        graded_space::{GradedLinearMap, GradedVectorSpace},
        kcomodule::{kCofreeComodule, kComodule},
        kmorphism::kComoduleMorphism,
    },
    traits::Coalgebra,
    types::{AGeneratorIndex, BasisElement, CoalgebraIndex, CoalgebraIndexType},
};

#[derive(Debug, Clone, PartialEq, Deserialize, Serialize, DeepSizeOf)]
#[allow(non_camel_case_types)]
pub struct kCoalgebra<G: Grading, F: Field, M: Matrix<F>> {
    pub space: GradedVectorSpace<G, BasisElement>,
    pub coaction: HashMap<CoalgebraIndex<G>, Vec<(CoalgebraIndex<G>, CoalgebraIndex<G>, F)>>,
    pub(crate) _p: PhantomData<M>,
}

impl<G: Grading, F: Field, M: Abelian<F>> Coalgebra<G> for kCoalgebra<G, F, M> {
    type BaseField = F;
    type BaseRing = F;
    type RingMorph = M;
    type Comod = kComodule<G, Self>;
    type CofMod = kCofreeComodule<G, Self>;
    type ComodMorph = kComoduleMorphism<G, F, M>;

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

        let space_map: HashMap<G, Vec<_>> = [(shift, vec![M::Generator::default()])]
            .into_iter()
            .collect();
        let space = GradedVectorSpace::from(space_map);

        let coact_map: HashMap<G, M> = [(shift, M::identity(1))].into_iter().collect();
        let coaction = GradedLinearMap::from(coact_map);

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

    fn cofree_comodule(
        &self,
        index: usize,
        shift: G,
        limit: G,
        generator: <Self::RingMorph as Abelian<Self::BaseRing>>::Generator,
    ) -> Self::CofMod {
        let space: HashMap<G, Vec<_>> = self
            .space
            .0
            .iter()
            .filter_map(|(g, v)| {
                let sum = *g + shift;
                if sum <= limit {
                    let k_basis: Vec<_> = (0..v.len())
                        .map(|j| {
                            (
                                ((*g, j as CoalgebraIndexType), index as AGeneratorIndex),
                                generator,
                            )
                        })
                        .collect();
                    Some((*g + shift, k_basis))
                } else {
                    None
                }
            })
            .collect();

        kCofreeComodule {
            space: GradedVectorSpace(space),
        }
    }
}

#[allow(non_snake_case)]
pub fn A0_coalgebra() -> kCoalgebra<UniGrading, F2, FlatMatrix<F2>> {
    let mut space = GradedVectorSpace::new();
    space.0.insert(
        UniGrading(0),
        vec![BasisElement {
            name: "1".to_owned(),
            excess: 0,
        }],
    );

    space.0.insert(
        UniGrading(1),
        vec![BasisElement {
            name: "xi1".to_owned(),
            excess: 0,
        }],
    );

    let mut coaction = HashMap::default();
    coaction.insert(
        (UniGrading::zero(), 0),
        vec![((UniGrading::zero(), 0), (UniGrading::zero(), 0), F2::one())],
    );
    coaction.insert(
        (UniGrading(1), 0),
        vec![
            ((UniGrading::zero(), 0), (UniGrading(1), 0), F2::one()),
            ((UniGrading(1), 0), (UniGrading::zero(), 0), F2::one()),
        ],
    );

    kCoalgebra {
        space,
        coaction,
        _p: PhantomData,
    }
}
