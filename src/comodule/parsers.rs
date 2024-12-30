use std::{collections::HashMap, fmt::Error};

use ahash::RandomState;
use itertools::Itertools;

use crate::linalg::{field::Field, graded::{BasisIndex, GradedLinearMap, GradedVectorSpace}, grading::Grading, matrix::Matrix};

use super::{kcoalgebra::kCoalgebra, kcomodule::kBasisElement, ktensor::kTensor};



impl<G: Grading, F: Field, M: Matrix<F>> kCoalgebra<G, F, M> {

    pub fn parse_direct(
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

    pub fn parse_polynomial_hopf_algebra(
        input: &str,
        max_grading: G,
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
            Generator,
            Relations,
            Coaction,
        }

        let mut state = State::None;
        let mut field: Option<usize> = None;
        let mut generators: Vec<(String, G)> = vec![];
        let mut relations: Vec<Monomial> = vec![];
        let mut coactions: Vec<Tensor<F>> = vec![];
        let mut generator_translate: HashMap<String, usize> = HashMap::new();

        for line in input.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            match line {
                _ if line.starts_with("- FIELD") => {
                    if state != State::None {
                        panic!("Field should be the first to parse");
                    }
                    state = State::Field;
                }
                _ if line.starts_with("- GENERATOR") => {
                    if state != State::Field {
                        panic!("Expected FIELD to be parsed first");
                    }
                    state = State::Generator;
                }
                _ if line.starts_with("- RELATION") => {
                    if state != State::Generator {
                        panic!("Expected GENERATOR to be parsed first");
                    }
                    state = State::Relations;
                }
                _ if line.starts_with("- COACTION") => {
                    if state != State::Relations {
                        panic!("Expected RELATIONS to be parsed first");
                    }
                    state = State::Coaction;
                }
                _ => match state {
                    State::Field => {
                        if field.is_some() {
                            panic!("Field already parsed");
                        }
                        field = Some(line.parse::<usize>().expect("Invalid field value"));
                    }
                    State::Generator => {
                        let (name, grade) = line.split_once(':').expect(format!("Invalid GENERATOR format, got: {}", line).as_str());
                        let grade = G::parse(grade.trim()).expect("Invalid grading");
                        generator_translate.insert(name.trim().to_string(), generators.len());
                        generators.push((name.trim().to_string(), grade));
                    }
                    State::Relations => {
                        let monomial = parse_monomial(line, &generator_translate, generators.len());
                        relations.push(monomial);
                    }
                    State::Coaction => {
                        let (name, tensors) = line.split_once(':').expect("Invalid COACTION format");
                        let tensors: Vec<_> = tensors.split('+').map(|t| {
                            let (l, r) = t.split_once('|').expect("Invalid tensor format");
                            (F::one(), parse_monomial(l.trim(), &generator_translate, generators.len()), parse_monomial(r.trim(), &generator_translate, generators.len()))
                        }).collect();
                        assert_eq!(name.trim(), generators[coactions.len()].0, "Coaction must match generator order");
                        coactions.push(tensors);
                    }
                    _ => panic!("Unexpected state"),
                },
            }
        }

        if state != State::Coaction {
            panic!("Coalgebra definition is not complete");
        }

        let n = generators.len();
        let one_monomial: Monomial = vec![0; n]; // All exponents zero => 1
        let mut basis_information: HashMap<Monomial, Tensor<F>> = HashMap::new();
        let mut queue: Vec<(Monomial, usize)> = vec![];

        // Initialize basis information for the unit monomial (1)
        basis_information.insert(one_monomial.clone(), vec![(F::one(), one_monomial.clone(), one_monomial.clone())]);

        // Populate the queue with initial values
        for i in 0..n {
            queue.push((one_monomial.clone(), i));
        }

        // BFS loop
        while let Some((current_monomial, generator_index)) = queue.pop() {
            if let Some(next_monomial) = multiply_monomial_by_generator(&current_monomial, generator_index, &relations) {
                // Check if `next_monomial` is already processed
                if !basis_information.contains_key(&next_monomial) {
                    // Check if the grade of the monomial is valid
                    let next_grade = monomial_to_grade(&next_monomial, &generators);
                    if next_grade <= max_grading {
                        // Calculate the coaction for the new monomial
                        let coaction_result = multiply_coaction_elements(
                            basis_information.get(&current_monomial).unwrap(),
                            &coactions[generator_index],
                            &relations,
                        );

                        // Store the result and push new states to the queue
                        basis_information.insert(next_monomial.clone(), coaction_result);
                        for i in 0..n {
                            queue.push((next_monomial.clone(), i));
                        }
                    }
                }
            }
        }
        // Construct the basis structure
        let mut basis: HashMap<G, Vec<kBasisElement>, RandomState> = HashMap::default();
        let mut monomial_to_grade_index: HashMap<Monomial, (G, usize), RandomState> = HashMap::default();

        for monomial in basis_information.keys() {
            let grade = monomial_to_grade(monomial, &generators);
            let label = monomial_to_string(monomial, &generators);

            let element = kBasisElement {
                name: label.clone(),
                generator: generators.contains(&(label, grade)),
                primitive: None,
                generated_index: 0,
            };

            let index = basis.entry(grade).or_insert_with(Vec::new).len();
            monomial_to_grade_index.insert(monomial.clone(), (grade, index));
            basis.entry(grade).or_insert_with(Vec::new).push(element);
        }

        // Create the graded vector space A
        let coalg_vector_space = GradedVectorSpace::from(basis.clone());

        // create A \otimes A
        let tensor = kTensor::generate(&coalg_vector_space, &coalg_vector_space);

        // Create the coaction map
        let mut coaction: HashMap<G, M, RandomState> = HashMap::default();

        for grade in tensor.dimensions.keys() { 
            
            // This loop is really weird, it should be a loop over basis_information, but I am too lazy to rewrite it
            // Todo: Create and store mutable 'map' for each grade, then iterate over basis_information 
            let tensor_rows = tensor.dimensions[grade];
            let basis_cols = &basis[grade].len();
            let mut map = M::zero(*basis_cols, tensor_rows);

            for (monomial, coaction_elements) in &basis_information { 
                let (basis_grade, basis_index) = monomial_to_grade_index.get(monomial).unwrap();
                if basis_grade != grade {
                    continue;
                }

                for (coeff, a, b) in coaction_elements {
                    let a_grade_index = monomial_to_grade_index.get(a).unwrap();
                    let b_grade_index = monomial_to_grade_index.get(b).unwrap();
                    let (_, tensor_index) = tensor.construct[&b_grade_index][&a_grade_index];

                    map.set(*basis_index, tensor_index, *coeff);
                }
            }

            coaction.insert(*grade, map);
        }
        
        let mut coalg = kCoalgebra {
            space: coalg_vector_space,
            tensor: tensor,
            coaction: GradedLinearMap::from(coaction),
        };

        coalg.set_primitives();

        // I think we don't need the generator translate in the rest of the code
        // but i left it here for now ¯\_(ツ)_/¯
        let empty_generator_translate = generator_translate.into_iter().map(|(name, _)| (name, (G::zero(), 0))).collect();

        Ok((coalg, empty_generator_translate))

    }
}

// Helper functions

type Monomial = Vec<usize>;
type Tensor<F> = Vec<(F, Monomial, Monomial)>;

fn parse_monomial(name: &str, generator_translate: &HashMap<String, usize>, size: usize) -> Monomial {
    let mut mon = vec![0; size];
    for el in name.split('*') {
        let parts: Vec<&str> = el.split('^').collect();
        let (name, expo) = match parts.len() {
            1 => (parts[0].trim(), 1),
            2 => (parts[0].trim(), parts[1].parse::<usize>().expect("Invalid exponent")),
            _ => panic!("Invalid monomial format"),
        };
        if name == "1" {
            continue;
        }
        let index = *generator_translate.get(name).expect("Generator not found");
        mon[index] = expo;
    }
    mon
}




// monomial opertations
fn increment_monomial(m: &Monomial, index: usize) -> Monomial {
    let mut new_monomial = m.clone();
    new_monomial[index] += 1;
    new_monomial
}

fn multiply_monomial_by_generator(m: &Monomial, index: usize, relations: &Vec<Monomial>) -> Option<Monomial> {
    let new_monomial = increment_monomial(m, index);
    reduce_monomial(&new_monomial, relations)
}

fn add_monomials(a: &Monomial, b: &Monomial) -> Monomial {
    a.iter().zip(b.iter()).map(|(x, y)| x + y).collect()
}

fn monomial_is_divisible(a: &Monomial, b: &Monomial) -> bool {
    a.iter().zip(b.iter()).all(|(x, y)| x >= y)
}

fn reduce_monomial(m: &Monomial, relations: &Vec<Monomial>) -> Option<Monomial> {
    for relation in relations {
        if monomial_is_divisible(m, relation) {
            return None;
        }
    }
    Some(m.clone())
}

fn monomial_to_grade<G: Grading>(m: &Monomial, generators: &Vec<(String, G)>) -> G {
    m.iter().zip(generators.iter()).map(|(x, (_, g))| g.integer_multiplication(*x as i32)).sum::<G>()
}

fn monomial_to_string<G: Grading>(m: &Monomial, generators: &Vec<(String, G)>) -> String {
    let generators = m.iter().zip(generators.iter()).filter(|(x, _)| **x > 0).map(|(x, (name, _))| {
        if *x == 1 {
            name.clone()
        } else {
            format!("{}^{}", name, x)
        }
    }).collect::<Vec<String>>(); 
    if generators.len() == 0 {
        1.to_string()
    } else {
        generators.join("*")
    }
}

fn multiply_monomials(
    a: &Monomial,
    b: &Monomial,
    relations: &Vec<Monomial>,
) -> Option<Monomial> {
    let product = add_monomials(a, b);
    reduce_monomial(&product, relations)
}

fn multiply_tensor_terms<F: Field>(
    a: &(F, Monomial, Monomial),
    b: &(F, Monomial, Monomial),
    relations: &Vec<Monomial>,
) -> Option<(F, Monomial, Monomial)> {
    let (a_coeff, a_left, a_right) = a;
    let (b_coeff, b_left, b_right) = b;

    let left_product = multiply_monomials(a_left, b_left, relations)?;
    let right_product = multiply_monomials(a_right, b_right, relations)?;

    let result_coeff = *a_coeff * *b_coeff;
    Some((result_coeff, left_product, right_product))
}

fn multiply_coaction_elements<F: Field>(
    a: &Tensor<F>,
    b: &Tensor<F>,
    relations: &Vec<Monomial>,
) -> Tensor<F> {
    let mut result: Vec<(F, Monomial, Monomial)> = vec![];

    for x in a {
        for y in b {
            if let Some(term) = multiply_tensor_terms(x, y, relations) {
                let (coeff, left, right) = term;
                if let Some(existing) = result.iter_mut().find(|(_, l, r)| l == &left && r == &right) {
                    existing.0 += coeff;
                } else {
                    result.push((coeff, left, right));
                }
            }
        }
    }

    result.retain(|(coeff, _, _)| *coeff != F::zero());
    result
}
