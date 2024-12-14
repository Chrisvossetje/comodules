use std::{collections::HashMap, sync::Arc};

use crate::linalg::{
    field::{Field, F2},
    graded::{BasisElement, GradedLinearMap, GradedVectorSpace, Grading, UniGrading},
    matrix::{FieldMatrix, Matrix},
};

use super::{comodule::Comodule, ktensor::kTensor};

#[derive(Debug, Clone, PartialEq, Default)]
#[allow(non_camel_case_types)]
pub struct kBasisElement {
    pub name: String,
    pub generator: bool,
    pub primitive: Option<usize>,
    pub generated_index: usize,
}

impl BasisElement for kBasisElement {}

#[derive(Debug, Clone, PartialEq)]
#[allow(non_camel_case_types)]
pub struct kCoalgebra<G: Grading, F: Field> {
    pub space: GradedVectorSpace<G, kBasisElement>,
    pub coaction: GradedLinearMap<G, F, FieldMatrix<F>>,
    pub tensor: kTensor<G>,
}

#[derive(Clone, PartialEq)]
#[allow(non_camel_case_types)]
pub struct kComodule<G: Grading, F: Field> {
    pub coalgebra: Arc<kCoalgebra<G, F>>,
    pub space: GradedVectorSpace<G, kBasisElement>,
    pub coaction: GradedLinearMap<G, F, FieldMatrix<F>>,
    pub tensor: kTensor<G>,
}

impl<G: Grading, F: Field>  std::fmt::Debug for kComodule<G,F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self.space.0)
    }
}

impl<G: Grading, F: Field> Comodule<G> for kComodule<G, F> {
    type Element = kBasisElement;
    type Coalgebra = kCoalgebra<G, F>;

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

        let space_map: HashMap<G, Vec<kBasisElement>> = [(zero, vec![el])].into_iter().collect();
        let space = GradedVectorSpace::from(space_map);

        let coact_map: HashMap<G, FieldMatrix<F>> =
            [(zero, FieldMatrix::identity(1))].into_iter().collect();
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

        let mut dimensions = HashMap::new();
        dimensions.insert(zero, 1);

        let mut construct = HashMap::new();
        let mut first_entry = HashMap::new();
        first_entry.insert((zero, 0), (zero, 0));
        construct.insert((zero, 0), first_entry);

        let mut deconstruct = HashMap::new();
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
                    self_els.append(other_els);
                })
                .or_insert(other_els.clone());
        });

        self.tensor.direct_sum(&mut other.tensor, &self_dimensions);
    }

    fn cofree_comodule(coalgebra: Arc<Self::Coalgebra>, index: usize, grade: G, limit: G) -> Self {
        let coaction: HashMap<G, FieldMatrix<F>> = coalgebra
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
        let tensor = coalgebra.tensor.add_and_restrict(grade, limit);

        kComodule {
            coalgebra,
            space: GradedVectorSpace(space),
            coaction: GradedLinearMap::from(coaction),
            tensor,
        }
    }
}

#[allow(non_snake_case)]
pub fn A0_coalgebra() -> kCoalgebra<UniGrading, F2> {
    let mut space = GradedVectorSpace::new();
    space.0.insert(
        0,
        vec![kBasisElement {
            name: "1".to_owned(),
            generator: true,
            primitive: None,
            generated_index: 0,
        }],
    );

    space.0.insert(
        1,
        vec![kBasisElement {
            name: "xi1".to_owned(),
            generator: false,
            primitive: Some(0),
            generated_index: 0,
        }],
    );

    let mut dimensions = HashMap::new();
    dimensions.insert(0, 1);
    dimensions.insert(1, 2);

    let mut construct = HashMap::new();
    let mut first_entry = HashMap::new();
    first_entry.insert((0, 0), (0, 0));
    first_entry.insert((1, 0), (1, 0));

    let mut second_entry = HashMap::new();
    second_entry.insert((0, 0), (1, 1));

    construct.insert((0, 0), first_entry);
    construct.insert((1, 0), second_entry);

    let mut deconstruct = HashMap::new();
    deconstruct.insert((0, 0), ((0, 0), (0, 0)));
    deconstruct.insert((1, 0), ((1, 0), (0, 0)));
    deconstruct.insert((1, 1), ((0, 0), (1, 0)));

    let tensor = kTensor {
        construct,
        deconstruct,
        dimensions,
    };

    let mut coaction = GradedLinearMap::empty();
    coaction.maps.insert(
        0,
        FieldMatrix {
            data: vec![vec![F2::one()]],
            domain: 1,
            codomain: 1,
        },
    );
    coaction.maps.insert(
        1,
        FieldMatrix {
            data: vec![vec![F2::one()], vec![F2::one()]],
            domain: 1,
            codomain: 2,
        },
    );

    kCoalgebra {
        space,
        coaction,
        tensor,
    }
}
