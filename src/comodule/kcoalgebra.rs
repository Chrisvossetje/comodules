use std::{collections::HashMap, fmt::Error};

use ahash::RandomState;
use itertools::Itertools;

use crate::linalg::{
    field::{Field, F2},
    graded::{BasisIndex, GradedLinearMap, GradedVectorSpace},
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

    pub fn parse(
        input: &str,
    ) -> Result<
        (
            kCoalgebra<G, F, M>,
            HashMap<String, BasisIndex<G>, RandomState>,
        ),
        Error,
    > {
        #[derive(Debug, Clone, PartialEq)]
        enum State {
            None,
            Field,
            Basis,
            Generator,
            Coaction,
        }

        let mut state = State::None;
        let mut field: Option<i32> = None;
        let mut basis: Vec<(String, G)> = vec![];
        let mut generator = vec![];
        let mut coaction_lut = vec![];

        for line in input.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            match line {
                _ if line.starts_with("- FIELD") => {
                    if state != State::None {
                        panic!("Field should be first to parse");
                    }
                    state = State::Field;
                }
                _ if line.starts_with("- BASIS") => {
                    if state != State::Field {
                        panic!("Expected FIELD to be parsed first");
                    }
                    state = State::Basis;
                }
                _ if line.starts_with("- GENERATOR") => {
                    if state != State::Basis {
                        panic!("Expected BASIS to be parsed first");
                    }
                    state = State::Generator;
                }
                _ if line.starts_with("- COACTION") => {
                    if state != State::Generator {
                        panic!("Expected GENERATOR to be parsed first");
                    }
                    state = State::Coaction;
                }
                _ => match state {
                    State::Field => {
                        if field.is_some() {
                            panic!("Field already parsed");
                        }
                        field = Some(line.parse::<i32>().expect("Invalid field value"));
                    }
                    State::Basis => {
                        let (name, grade) = line.split_once(":").expect("Invalid BASIS format");
                        basis.push((name.trim().to_string(), (G::parse(grade.trim()).unwrap())));
                    }
                    State::Generator => {
                        generator.push(line.to_string());
                    }
                    State::Coaction => {
                        let (name, tensors) =
                            line.split_once(":").expect("Invalid COACTION format");
                        let mut ts = vec![];
                        for t in tensors.split('+') {
                            let t = t.trim();
                            if let Some(f) = field {
                                if f == 2 {
                                    let (l, r) = t.split_once('|').expect("Invalid tensor format");
                                    ts.push((
                                        1.to_string(),
                                        l.trim().to_string(),
                                        r.trim().to_string(),
                                    ));
                                } else {
                                    let (v, rest) =
                                        t.split_once('*').expect("Invalid tensor scalar");
                                    let (l, r) =
                                        rest.split_once('|').expect("Invalid tensor format");
                                    ts.push((
                                        v.trim().to_string(),
                                        l.trim().to_string(),
                                        r.trim().to_string(),
                                    ));
                                }
                            }
                        }
                        coaction_lut.push((name.trim().to_string(), ts));
                    }
                    _ => panic!("Unexpected state"),
                },
            }
        }

        // Verify state
        if state != State::Coaction {
            panic!("Coalgebra definition is not complete");
        }

        // Create basis dictionary
        let mut basis_dict: HashMap<String, (kBasisElement, G), RandomState> = HashMap::default();
        for (name, grade) in basis.iter() {
            if basis_dict.contains_key(name) {
                panic!("Name in basis appears twice");
            }
            basis_dict.insert(
                name.clone(),
                (
                    kBasisElement {
                        name: name.clone(),
                        generator: generator.contains(name),
                        primitive: None,
                        generated_index: 0,
                    },
                    *grade,
                ),
            );
        }

        // Transform basis
        let mut transformed = HashMap::default();
        let mut basis_translate = HashMap::default();

        for (name, (el, gr)) in basis_dict.iter().sorted_by_key(|(name, _)| *name) {
            transformed.entry(*gr).or_insert(vec![]).push(el.clone());
            basis_translate.insert(name.clone(), (*gr, transformed[&gr].len() - 1));
        }

        let graded_space = GradedVectorSpace(transformed);
        let tensor = kTensor::generate(&graded_space, &graded_space);

        let mut coaction: HashMap<G, M, RandomState> = HashMap::default();
        for (gr, elements) in &graded_space.0 {
            let domain = elements.len();
            let codomain = *tensor.dimensions.get(gr).unwrap_or(&0);
            coaction.insert(*gr, M::zero(domain, codomain));
        }

        for (b, ls) in coaction_lut {
            let (gr, id) = basis_translate[&b];
            for (scalar, l, r) in ls {
                let l_id = basis_translate[&l];
                let r_id = basis_translate[&r];
                assert_eq!((l_id.0 + r_id.0), gr, "Grades are not homogenous");
                let t_id = tensor.construct[&r_id][&l_id];
                assert_eq!(t_id.0, gr, "Grades are not homogenous");
                coaction
                    .get_mut(&gr)
                    .unwrap()
                    .set(id, t_id.1, F::parse(&scalar).unwrap());
            }
        }

        let mut coalg = kCoalgebra {
            space: graded_space,
            tensor: tensor,
            coaction: GradedLinearMap::from(coaction),
        };

        coalg.set_primitives();

        Ok((coalg, basis_translate))
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
