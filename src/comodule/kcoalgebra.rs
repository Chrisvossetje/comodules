use core::panic;
use std::collections::HashMap;

use itertools::Itertools;

use crate::linalg::{
    field::{Field, F2},
    graded::{GradedLinearMap, GradedVectorSpace},
    grading::{Grading, UniGrading},
    matrix::Matrix,
    row_matrix::RowMatrix,
};

use super::{kcomodule::kBasisElement, ktensor::kTensor};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[allow(non_camel_case_types)]
pub struct kCoalgebra<G: Grading, F: Field, M: Matrix<F>> {
    pub space: GradedVectorSpace<G, kBasisElement>,
    pub coaction: GradedLinearMap<G, F, M>,
    pub tensor: kTensor<G>,
}

impl<G: Grading, F: Field, M: Matrix<F>> kCoalgebra<G, F, M> {
    pub fn set_primitives(&mut self) {
        let mut primitive_index = 0;

        for (grade, basis_elements) in self.space.0.iter_mut().sorted_by_key(|(g, _)| *g) {
            let coact_map = &self.coaction.maps[grade];
            for (index, el) in basis_elements.iter_mut().enumerate() {
                let mut non_zero_count = 0;

                // Count non-zero entries
                for t_id in 0..coact_map.codomain() {
                    if !coact_map.get(index, t_id).is_zero() {
                        non_zero_count += 1;
                    }
                }

                // Check if exactly 2 non-zero entries
                if non_zero_count == 2 {
                    el.primitive = Some(primitive_index);
                    primitive_index += 1;
                }
            }
        }
    }

    pub fn set_generator(&mut self) {
        let grade_zero = self.space.0.get_mut(&G::zero());
        if let Some(basis) = grade_zero {
            if basis.len() == 1 {
                if let Some(basis_element) = basis.first_mut() {
                    basis_element.generator = true;
                }
            } else {
                panic!("Coalgebra is not connected, no unique generator found in grade 0");
            }
        } else {
            panic!("Grade 0 not found in space");
        }
    }
}

#[allow(non_snake_case)]
pub fn A0_coalgebra() -> kCoalgebra<UniGrading, F2, RowMatrix<F2>> {
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

    let mut dimensions = HashMap::default();
    dimensions.insert(0, 1);
    dimensions.insert(1, 2);

    let mut construct = HashMap::default();
    let mut first_entry = HashMap::default();
    first_entry.insert((0, 0), (0, 0));
    first_entry.insert((1, 0), (1, 1));

    let mut second_entry = HashMap::default();
    second_entry.insert((0, 0), (1, 0));

    construct.insert((0, 0), first_entry);
    construct.insert((1, 0), second_entry);

    let mut deconstruct = HashMap::default();
    deconstruct.insert((0, 0), ((0, 0), (0, 0)));
    deconstruct.insert((1, 1), ((1, 0), (0, 0)));
    deconstruct.insert((1, 0), ((0, 0), (1, 0)));

    let tensor = kTensor {
        construct,
        deconstruct,
        dimensions,
    };

    let mut coaction = GradedLinearMap::empty();
    coaction.maps.insert(
        0,
        RowMatrix {
            data: vec![vec![F2::one()]],
            domain: 1,
            codomain: 1,
        },
    );
    coaction.maps.insert(
        1,
        RowMatrix {
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
