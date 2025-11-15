use std::{ops::AddAssign, sync::Arc};

use ahash::HashMap;
use itertools::Itertools;

use crate::{
    basiselement::kBasisElement,
    grading::{Grading, OrderedGrading},
    linalg::{
        field::Field,
        graded::{BasisIndex, GradedLinearMap, GradedVectorSpace},
        matrix::Matrix,
    }, tensor::Tensor,
};

use super::{
    kcoalgebra::kCoalgebra,
    kcomodule::kComodule,
};

impl<G: Grading + OrderedGrading, F: Field, M: Matrix<F>> kCoalgebra<G, F, M> {
    fn check_translator(&self, translate: &HashMap<String, BasisIndex<G>>) -> bool {
        for (gr, space) in &self.space.0 {
            for (index, el) in space.iter().enumerate() {
                let comp = translate.get(&el.name);
                assert!(
                    comp.is_some(),
                    "Name: {} was not found in translator",
                    el.name
                );
                assert_eq!(comp.unwrap().0, *gr, "Grades do not coincide");
                assert_eq!(comp.unwrap().1, index, "Indices do not coincide");
            }
        }
        return true;
    }

    pub fn parse(
        input: &str,
        max_grading: G,
    ) -> Result<
        (
            kCoalgebra<G, F, M>,
            HashMap<String, BasisIndex<G>>,
        ),
        String,
    > {
        if input.contains("- BASIS") {
            Self::parse_direct(input)
        } else {
            Self::parse_polynomial_hopf_algebra(input, max_grading)
        }
    }

    fn parse_direct(
        input: &str,
    ) -> Result<
        (
            kCoalgebra<G, F, M>,
            HashMap<String, BasisIndex<G>>,
        ),
        String,
    > {
        #[derive(Debug, Clone, PartialEq)]
        enum State {
            None,
            Field,
            Basis,
            Coaction,
        }
        let mut state = State::None;
        let mut field: Option<i32> = None;
        let mut basis: Vec<(String, G)> = vec![];
        let mut coaction_lut = vec![];

        for (line_num, line) in input.lines().enumerate() {
            let line_num = line_num + 1;
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            match line {
                _ if line.starts_with("- FIELD") => {
                    if state != State::None {
                        return Err(format!("Line {}: Field should be first to parse", line_num));
                    }
                    state = State::Field;
                }
                _ if line.starts_with("- BASIS") => {
                    if state != State::Field {
                        return Err(format!(
                            "Line {}: Expected FIELD to be parsed first",
                            line_num
                        ));
                    }
                    state = State::Basis;
                }
                _ if line.starts_with("- COACTION") => {
                    if state != State::Basis {
                        return Err(format!(
                            "Line {}: Expected BASIS to be parsed first",
                            line_num
                        ));
                    }
                    state = State::Coaction;
                }
                _ => match state {
                    State::Field => {
                        if field.is_some() {
                            return Err(format!("Line {}: Field already parsed", line_num));
                        }
                        field = Some(line.parse::<i32>().map_err(|_| {
                            format!("Line {}: Invalid field value '{}'", line_num, line)
                        })?);
                        if field.unwrap() as usize != F::get_characteristic() {
                            return Err(format!(
                                "Line {}: Field does not have the expected characteristic",
                                line_num
                            ));
                        }
                    }
                    State::Basis => {
                        let (name, grade) = line.split_once(":").ok_or(format!(
                            "Line {}: Invalid BASIS format '{}' - expected 'name:grade'",
                            line_num, line
                        ))?;
                        basis.push((
                            name.trim().to_string(),
                            G::parse(grade.trim()).map_err(|e| {
                                format!(
                                    "Line {}: Invalid grade '{}' - {}",
                                    line_num,
                                    grade.trim(),
                                    e
                                )
                            })?,
                        ));
                    }
                    State::Coaction => {
                        let (name, tensors) = line.split_once(":").ok_or(format!(
                            "Line {}: Invalid COACTION format '{}' - expected 'name:tensors'",
                            line_num, line
                        ))?;
                        let mut ts = vec![];
                        for t in tensors.split('+') {
                            let t = t.trim();
                            let (s, t) = match t.split_once('.') {
                                Some((s, t)) => (s, t),
                                None => ("1", t),
                            };

                            let (l, r) = t
                                .split_once('|')
                                .ok_or(format!("Line {}: Invalid tensor format '{}' - expected 'left|right' or scalar.left|right", line_num, t))?;
                            ts.push((
                                s.trim().to_string(),
                                l.trim().to_string(),
                                r.trim().to_string(),
                            ));
                        }
                        coaction_lut.push((name.trim().to_string(), ts));
                    }
                    _ => return Err(format!("Line {}: Unexpected state", line_num)),
                },
            }
        }

        // Verify state
        if state != State::Coaction {
            return Err("Coalgebra definition is not complete - missing sections".to_owned());
        }

        // Create basis dictionary
        let mut basis_dict: HashMap<String, (kBasisElement, G)> = HashMap::default();
        for (name, grade) in basis.iter() {
            if basis_dict.contains_key(name) {
                return Err(format!("Basis element '{}' appears twice", name));
            }
            basis_dict.insert(
                name.clone(),
                (
                    kBasisElement {
                        name: name.clone(),
                        generator: false,
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
        let tensor = Tensor::generate(&graded_space, &graded_space);

        let mut coaction: HashMap<G, M> = HashMap::default();
        for (gr, elements) in &graded_space.0 {
            let domain = elements.len();
            let codomain = *tensor.dimensions.get(gr).unwrap_or(&0);
            coaction.insert(*gr, M::zero(domain, codomain));
        }

        for (b, ls) in coaction_lut {
            let (gr, id) = basis_translate.get(&b).ok_or(format!(
                "Basis element '{}' not found in coaction definition",
                b
            ))?;
            for (scalar, l, r) in ls {
                let l_id = basis_translate.get(&l).ok_or(format!(
                    "Left element '{}' not found in basis for coaction of '{}'",
                    l, b
                ))?;
                let r_id = basis_translate.get(&r).ok_or(format!(
                    "Right element '{}' not found in basis for coaction of '{}'",
                    r, b
                ))?;
                if (l_id.0 + r_id.0) != *gr {
                    return Err(format!(
                        "Grades are not homogenous for coaction of '{}': {} + {} != {}",
                        b, l_id.0, r_id.0, gr
                    ));
                };
                let t_id = tensor.construct[&r_id][&l_id];
                if (t_id.0) != *gr {
                    return Err(format!(
                        "Tensor grades are not homogenous for coaction of '{}': {} != {}",
                        b, t_id.0, gr
                    ));
                };
                coaction
                    .get_mut(&gr)
                    .ok_or(format!("Expected coaction to exist in grade {}", gr))?
                    .set(
                        *id,
                        t_id.1,
                        F::parse(&scalar).map_err(|e| {
                            format!("Invalid scalar '{}' for coaction of '{}': {}", scalar, b, e)
                        })?,
                    );
            }
        }

        let mut coalg = kCoalgebra {
            space: graded_space,
            tensor: tensor,
            coaction: GradedLinearMap::from(coaction),
        };

        coalg.set_primitives();
        coalg.set_generator()?;
        coalg.reduce();

        debug_assert!(Self::check_translator(&coalg, &basis_translate));

        Ok((coalg, basis_translate))
    }

    fn parse_polynomial_hopf_algebra(
        input: &str,
        max_grading: G,
    ) -> Result<
        (
            kCoalgebra<G, F, M>,
            HashMap<String, BasisIndex<G>>,
        ),
        String,
    > {
        #[derive(Debug, Clone, PartialEq)]
        enum State {
            None,
            Field,
            Generator,
            Relations,
            Coaction,
        }

        let max_grading = max_grading.incr().incr();
        let mut state = State::None;
        let mut field: Option<usize> = None;
        let mut generators: Vec<(String, G)> = vec![];
        let mut relations: Vec<Monomial> = vec![];
        let mut coactions: Vec<HelperTensor<F>> = vec![];
        let mut generator_translate: HashMap<String, usize> = HashMap::default();
        let mut basis_translate: HashMap<String, BasisIndex<G>> = HashMap::default();

        for (line_num, line) in input.lines().enumerate() {
            let line_num = line_num + 1;
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            match line {
                _ if line.starts_with("- FIELD") => {
                    if state != State::None {
                        return Err(format!(
                            "Line {}: Field should be the first to parse",
                            line_num
                        ));
                    }
                    state = State::Field;
                }
                _ if line.starts_with("- GENERATOR") => {
                    if state != State::Field {
                        return Err(format!(
                            "Line {}: Expected FIELD to be parsed first",
                            line_num
                        ));
                    }
                    state = State::Generator;
                }
                _ if line.starts_with("- RELATION") => {
                    if state != State::Generator {
                        return Err(format!(
                            "Line {}: Expected GENERATOR to be parsed first",
                            line_num
                        ));
                    }
                    state = State::Relations;
                }
                _ if line.starts_with("- COACTION") => {
                    if state != State::Relations {
                        return Err(format!(
                            "Line {}: Expected RELATIONS to be parsed first",
                            line_num
                        ));
                    }
                    state = State::Coaction;
                }
                _ => match state {
                    State::Field => {
                        if field.is_some() {
                            return Err(format!("Line {}: Field already parsed", line_num));
                        }
                        field = Some(line.parse::<usize>().map_err(|_| {
                            format!("Line {}: Invalid field value '{}'", line_num, line)
                        })?);
                        if field.unwrap() as usize != F::get_characteristic() {
                            return Err(format!(
                                "Line {}: Field does not have the expected characteristic",
                                line_num
                            ));
                        }
                    }
                    State::Generator => {
                        let (name, grade) = line.split_once(':').ok_or(format!(
                            "Line {}: Invalid GENERATOR format '{}' - expected 'name:grade'",
                            line_num, line
                        ))?;
                        let grade = G::parse(grade.trim()).map_err(|e| {
                            format!(
                                "Line {}: Invalid grade '{}' - {}",
                                line_num,
                                grade.trim(),
                                e
                            )
                        })?;
                        if grade == G::zero() {
                            return Err(format!(
                                "Line {}: Grade of generator cannot be zero",
                                line_num
                            ));
                        }
                        generator_translate.insert(name.trim().to_string(), generators.len());
                        generators.push((name.trim().to_string(), grade));
                    }
                    State::Relations => {
                        let monomial = parse_monomial(line, &generator_translate, generators.len())
                            .map_err(|e| {
                                format!("Line {}: Invalid relation '{}' - {}", line_num, line, e)
                            })?;
                        relations.push(monomial);
                    }
                    State::Coaction => {
                        let (name, tensors) = line.split_once(':').ok_or(format!(
                            "Line {}: Invalid COACTION format '{}' - expected 'name:tensors'",
                            line_num, line
                        ))?;
                        let tensors: Vec<(F, Vec<usize>, Vec<usize>)> = tensors
                            .split('+')
                            .map::<Result<(F, Vec<usize>, Vec<usize>), String>, _>(|t| {
                                let (s,t) = match t.split_once('.') {
                                    Some((s, t)) => {
                                        (F::parse(s.trim()).map_err(|e| {
                                            format!("Line {}: Scalar {s} with tensor {t} could not be parsed. {e}", line_num)
                                        })?, t)
                                    },
                                    None => {
                                        (F::one(), t)
                                    }
                                };
                                let (l, r) = t
                                    .split_once('|')
                                    .ok_or(format!("Line {}: Invalid tensor format '{}' - expected 'left|right'", line_num, t))?;

                                let left = parse_monomial(
                                    l.trim(),
                                    &generator_translate,
                                    generators.len(),
                                ).map_err(|e| format!("Line {}: Invalid left monomial '{}' - {}", line_num, l.trim(), e))?;
                                // TODO: RHS can only be monomial ?
                                let right = parse_monomial(
                                    r.trim(),
                                    &generator_translate,
                                    generators.len(),
                                ).map_err(|e| format!("Line {}: Invalid right monomial '{}' - {}", line_num, r.trim(), e))?;
                                Ok((s, left, right))
                            })
                            .try_collect()?;
                        if name.trim() != generators[coactions.len()].0 {
                            return Err(format!("Line {}: Coaction for '{}' must match generator order, expected '{}'", line_num, name.trim(), generators[coactions.len()].0));
                        }

                        coactions.push(tensors);
                    }
                    _ => return Err(format!("Line {}: Unexpected state", line_num)),
                },
            }
        }

        if state != State::Coaction {
            return Err("Coalgebra definition is not complete - missing sections".to_owned());
        }

        let n = generators.len();
        let one_monomial: Monomial = vec![0; n]; // All exponents zero => 1
        let mut monomial_coaction: HashMap<Monomial, HelperTensor<F>> = HashMap::default();
        let mut queue: Vec<Monomial> = vec![one_monomial.clone()];

        // Initialize basis information for the unit monomial (1)
        monomial_coaction.insert(
            one_monomial.clone(),
            vec![(F::one(), one_monomial.clone(), one_monomial.clone())],
        );

        basis_translate.insert("1".to_owned(), (G::zero(), 0));

        println!("Entering BFS");
        let mut i = 0;
        // BFS loop
        while let Some(current_monomial) = queue.pop() {
            for generator_index in 0..n {
                if let Some(next_monomial) =
                    multiply_monomial_by_generator(&current_monomial, generator_index, &relations)
                {
                    i += 1;
                    if i % 100 == 0 {
                        println!(
                            "Processed {} elements, current basis_information length: {}",
                            i,
                            monomial_coaction.len()
                        );
                    }
                    // Check if `next_monomial` is already processed
                    if !monomial_coaction.contains_key(&next_monomial) {
                        // Check if the grade of the monomial is valid
                        let next_grade = monomial_to_grade(&next_monomial, &generators);
                        if next_grade <= max_grading {
                            // Calculate the coaction for the new monomial
                            let coaction_result = multiply_coaction_elements(
                                monomial_coaction.get(&current_monomial).ok_or(format!(
                                    "Basis monomial '{}' could not be found in queue",
                                    monomial_to_string(&current_monomial, &generators)
                                ))?,
                                &coactions[generator_index],
                                &relations,
                            );

                            // Store the result and push new state to the queue
                            monomial_coaction.insert(next_monomial.clone(), coaction_result);

                            queue.push(next_monomial.clone());
                        }
                    }
                }
            }
        }
        // Construct the basis structure
        let mut basis: HashMap<G, Vec<kBasisElement>> = HashMap::default();
        let mut monomial_to_grade_index: HashMap<Monomial, (G, usize)> =
            HashMap::default();

        for monomial in monomial_coaction.keys().sorted() {
            let grade = monomial_to_grade(monomial, &generators);
            let label = monomial_to_string(monomial, &generators);

            let element = kBasisElement {
                name: label.clone(),
                generator: false,
                primitive: None,
                generated_index: 0,
            };

            let index = basis.entry(grade).or_insert_with(Vec::new).len();
            monomial_to_grade_index.insert(monomial.clone(), (grade, index));

            let index = basis.get(&grade).map_or(0, |el| el.len());
            basis_translate.insert(label, (grade, index));
            basis.entry(grade).or_insert_with(Vec::new).push(element);
        }

        // Create the graded vector space A
        let coalg_vector_space = GradedVectorSpace::from(basis.clone());

        // create A \otimes A
        let tensor = Tensor::generate(&coalg_vector_space, &coalg_vector_space);

        // Create the coaction map
        let mut coaction: HashMap<G, M> = HashMap::default();

        for (grade, els) in basis.iter() {
            let tensor_rows = tensor.dimensions[grade];
            let basis_cols = els.len();
            let map = M::zero(basis_cols, tensor_rows);
            coaction.insert(*grade, map);
        }

        for (monomial, coaction_elements) in &monomial_coaction {
            let (basis_grade, basis_index) =
                monomial_to_grade_index.get(monomial).ok_or(format!(
                    "Expected monomial '{}' to exist in lookup",
                    monomial_to_string(monomial, &generators)
                ))?;

            let map = coaction.get_mut(basis_grade).ok_or(format!("Expected a coaction to exist in dimension {basis_grade}. For the element {:?}, {:?}", monomial, coaction_elements))?;
            for (coeff, a, b) in coaction_elements {
                let a_grade_index = monomial_to_grade_index.get(a).ok_or(format!(
                    "Expected left monomial '{}' to exist when constructing coaction",
                    monomial_to_string(a, &generators)
                ))?;
                let b_grade_index = monomial_to_grade_index.get(b).ok_or(format!(
                    "Expected right monomial '{}' to exist when constructing coaction",
                    monomial_to_string(b, &generators)
                ))?;
                let (_, tensor_index) = tensor.construct[&b_grade_index][&a_grade_index];

                map.set(*basis_index, tensor_index, *coeff);
            }
        }

        let mut coalg = kCoalgebra {
            space: coalg_vector_space,
            tensor: tensor,
            coaction: GradedLinearMap::from(coaction),
        };

        coalg.set_primitives();
        coalg.set_generator()?;
        coalg.reduce();

        debug_assert!(Self::check_translator(&coalg, &basis_translate));

        Ok((coalg, basis_translate))
    }
}

// Helper functions

type Monomial = Vec<usize>;
type HelperTensor<F> = Vec<(F, Monomial, Monomial)>;

fn parse_name_exponent(el: &str) -> Result<(&str, usize), String> {
    let parts: Vec<&str> = el.split('^').collect();
    Ok(match parts.len() {
        1 => (parts[0].trim(), 1),
        2 => (
            parts[0].trim(),
            parts[1].parse::<usize>().map_err(|_| {
                format!(
                    "Invalid exponent '{}' in monomial element '{}'",
                    parts[1], el
                )
            })?,
        ),
        _ => {
            return Err(format!(
                "Invalid monomial format '{}' - expected 'name' or 'name^exponent'",
                el
            ))
        }
    })
}

fn parse_monomial(
    name: &str,
    generator_translate: &HashMap<String, usize>,
    size: usize,
) -> Result<Monomial, String> {
    let mut mon = vec![0; size];
    for el in name.split(',') {
        let (name, expo) = parse_name_exponent(el)?;
        if name == "1" {
            continue;
        }
        let index = *generator_translate.get(name).ok_or(format!(
            "Generator '{}' not found in monomial '{}'",
            name, name
        ))?;
        mon[index] = expo;
    }
    Ok(mon)
}

// monomial opertations
fn increment_monomial(m: &Monomial, index: usize) -> Monomial {
    let mut new_monomial = m.clone();
    new_monomial[index] += 1;
    new_monomial
}

fn multiply_monomial_by_generator(
    m: &Monomial,
    index: usize,
    relations: &Vec<Monomial>,
) -> Option<Monomial> {
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
    m.iter()
        .zip(generators.iter())
        .map(|(x, (_, g))| g.integer_multiplication(*x as i32))
        .sum::<G>()
}

fn monomial_to_string<G: Grading>(m: &Monomial, generators: &Vec<(String, G)>) -> String {
    let generators = m
        .iter()
        .zip(generators.iter())
        .filter(|(x, _)| **x > 0)
        .map(|(x, (name, _))| {
            if *x == 1 {
                name.clone()
            } else {
                format!("{}^{}", name, x)
            }
        })
        .collect::<Vec<String>>();
    if generators.len() == 0 {
        "1".to_owned()
    } else {
        generators.join(",")
    }
}

fn multiply_monomials(a: &Monomial, b: &Monomial, relations: &Vec<Monomial>) -> Option<Monomial> {
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
    a: &HelperTensor<F>,
    b: &HelperTensor<F>,
    relations: &Vec<Monomial>,
) -> HelperTensor<F> {
    let mut result: Vec<(F, Monomial, Monomial)> = vec![];

    for x in a {
        for y in b {
            if let Some(term) = multiply_tensor_terms(x, y, relations) {
                let (coeff, left, right) = term;
                if let Some(existing) = result
                    .iter_mut()
                    .find(|(_, l, r)| l == &left && r == &right)
                {
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

impl<G: Grading + OrderedGrading, F: Field, M: Matrix<F>> kComodule<G, F, M> {
    pub fn parse(
        input: &str,
        coalgebra: Arc<kCoalgebra<G, F, M>>,
        coalgebra_translate: &HashMap<String, BasisIndex<G>>,
        max_grading: G,
    ) -> Result<kComodule<G, F, M>, String> {
        if input.contains("- BASIS") {
            Self::parse_direct(input, coalgebra, coalgebra_translate)
        } else {
            Self::parse_polynomial(input, coalgebra, coalgebra_translate, max_grading)
        }
    }

    fn parse_direct(
        input: &str,
        coalgebra: Arc<kCoalgebra<G, F, M>>,
        coalgebra_translate: &HashMap<String, BasisIndex<G>>,
    ) -> Result<kComodule<G, F, M>, String> {
        #[derive(Debug, Clone, PartialEq)]
        enum State {
            None,
            Basis,
            Coaction,
        }

        let mut state = State::None;
        let mut basis: Vec<(String, G)> = vec![];
        let mut coaction_lut = vec![];

        for (line_num, line) in input.lines().enumerate() {
            let line_num = line_num + 1;
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            match line {
                _ if line.starts_with("- BASIS") => {
                    if state != State::None {
                        return Err(format!(
                            "Line {}: Expected FIELD to be parsed first",
                            line_num
                        ));
                    }
                    state = State::Basis;
                }
                _ if line.starts_with("- COACTION") => {
                    if state != State::Basis {
                        return Err(format!(
                            "Line {}: Expected BASIS to be parsed first",
                            line_num
                        ));
                    }
                    state = State::Coaction;
                }
                _ => match state {
                    State::Basis => {
                        let (name, grade) = line.split_once(":").ok_or(format!(
                            "Line {}: Invalid BASIS format '{}' - expected 'name:grade'",
                            line_num, line
                        ))?;
                        basis.push((
                            name.trim().to_string(),
                            G::parse(grade.trim()).map_err(|e| {
                                format!(
                                    "Line {}: Invalid grade '{}' - {}",
                                    line_num,
                                    grade.trim(),
                                    e
                                )
                            })?,
                        ));
                    }
                    State::Coaction => {
                        let (name, tensors) = line.split_once(":").ok_or(format!(
                            "Line {}: Invalid COACTION format '{}' - expected 'name:tensors'",
                            line_num, line
                        ))?;
                        let mut ts = vec![];
                        for t in tensors.split('+') {
                            let t = t.trim();
                            let (s, t) = match t.split_once('.') {
                                Some((s, t)) => (s, t),
                                None => ("1", t),
                            };

                            let (l, r) = t
                                .split_once('|')
                                .ok_or(format!("Line {}: Invalid tensor format '{}' - expected 'left|right' or scalar.left|right", line_num, t))?;
                            ts.push((
                                s.trim().to_string(),
                                l.trim().to_string(),
                                r.trim().to_string(),
                            ));
                        }
                        coaction_lut.push((name.trim().to_string(), ts));
                    }
                    _ => return Err(format!("Line {}: Unexpected state", line_num)),
                },
            }
        }

        if state != State::Coaction {
            return Err("Comodule definition is not complete - missing sections".to_owned());
        }

        // Create basis dictionary
        let mut basis_dict: HashMap<String, (kBasisElement, G)> = HashMap::default();
        for (name, grade) in basis.iter() {
            if basis_dict.contains_key(name) {
                return Err(format!("Basis element '{}' appears twice", name));
            }
            basis_dict.insert(
                name.clone(),
                (
                    kBasisElement {
                        name: name.clone(),
                        generator: false,
                        primitive: None,
                        generated_index: 0,
                    },
                    *grade,
                ),
            );
        }

        // Transform basis
        let mut transformed = HashMap::default();
        let mut basis_translate: HashMap<String, BasisIndex<G>> = HashMap::default();

        for (name, (el, gr)) in basis_dict.iter().sorted_by_key(|(name, _)| *name) {
            transformed.entry(*gr).or_insert(vec![]).push(el.clone());
            basis_translate.insert(name.clone(), (*gr, transformed[&gr].len() - 1));
        }

        let graded_space = GradedVectorSpace(transformed);
        let tensor = Tensor::generate(&coalgebra.space, &graded_space);

        let mut coaction: HashMap<G, M> = HashMap::default();
        for (gr, elements) in &graded_space.0 {
            let domain = elements.len();
            let codomain = *tensor.dimensions.get(gr).unwrap_or(&0);
            coaction.insert(*gr, M::zero(domain, codomain));
        }

        for (b, ls) in coaction_lut {
            let (gr, id) = basis_translate.get(&b).ok_or(format!(
                "Basis element '{}' not found in coaction definition",
                b
            ))?;
            for (scalar, l, r) in ls {
                let l_id = coalgebra_translate.get(&l).ok_or(format!(
                    "Left element '{}' not found in coalgebra for coaction of '{}'",
                    l, b
                ))?;
                let r_id = basis_translate.get(&r).ok_or(format!(
                    "Right element '{}' not found in comodule basis for coaction of '{}'",
                    r, b
                ))?;

                if (l_id.0 + r_id.0) != *gr {
                    return Err(format!(
                        "Grades are not homogenous for coaction of '{}': {} + {} != {}",
                        b, l_id.0, r_id.0, gr
                    ));
                }

                let t_id = tensor.construct[&r_id][&l_id];
                if t_id.0 != *gr {
                    return Err(format!(
                        "Tensor grades are not homogenous for coaction of '{}': {} != {}",
                        b, t_id.0, gr
                    ));
                }

                coaction
                    .get_mut(&gr)
                    .ok_or(format!("Expected coaction to exist in grade {}", gr))?
                    .set(
                        *id,
                        t_id.1,
                        F::parse(&scalar).map_err(|e| {
                            format!("Invalid scalar '{}' for coaction of '{}': {}", scalar, b, e)
                        })?,
                    );
            }
        }

        Ok(kComodule::new(
            coalgebra,
            graded_space,
            GradedLinearMap::from(coaction),
            tensor,
        ))
    }

    fn parse_polynomial(
        input: &str,
        coalgebra: Arc<kCoalgebra<G, F, M>>,
        coalgebra_translate: &HashMap<String, BasisIndex<G>>,
        max_grading: G,
    ) -> Result<kComodule<G, F, M>, String> {
        #[derive(Debug, Clone, PartialEq)]
        enum State {
            None,
            Generator,
            Relations,
            Coaction,
        }

        let max_grading = max_grading.incr().incr();
        let mut state = State::None;
        let mut generators: Vec<(String, G)> = vec![];
        let mut relations: Vec<Monomial> = vec![];
        let mut coactions: Vec<Vec<(F, HashMap<String, usize>, Vec<usize>)>> = vec![];
        let mut generator_translate: HashMap<String, usize> = HashMap::default();

        for (line_num, line) in input.lines().enumerate() {
            let line_num = line_num + 1;
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            match line {
                _ if line.starts_with("- GENERATOR") => {
                    if state != State::None {
                        return Err(format!(
                            "Line {}: Expected FIELD to be parsed first",
                            line_num
                        ));
                    }
                    state = State::Generator;
                }
                _ if line.starts_with("- RELATION") => {
                    if state != State::Generator {
                        return Err(format!(
                            "Line {}: Expected GENERATOR to be parsed first",
                            line_num
                        ));
                    }
                    state = State::Relations;
                }
                _ if line.starts_with("- COACTION") => {
                    if state != State::Relations {
                        return Err(format!(
                            "Line {}: Expected RELATION to be parsed first",
                            line_num
                        ));
                    }
                    state = State::Coaction;
                }
                _ => match state {
                    State::Generator => {
                        let (name, grade) = line.split_once(':').ok_or(format!(
                            "Line {}: Invalid GENERATOR format '{}' - expected 'name:grade'",
                            line_num, line
                        ))?;
                        let grade = G::parse(grade.trim()).map_err(|e| {
                            format!(
                                "Line {}: Invalid grade '{}' - {}",
                                line_num,
                                grade.trim(),
                                e
                            )
                        })?;
                        generator_translate.insert(name.trim().to_string(), generators.len());
                        generators.push((name.trim().to_string(), grade));
                    }
                    State::Relations => {
                        let monomial = parse_monomial(line, &generator_translate, generators.len())
                            .map_err(|e| {
                                format!("Line {}: Invalid relation '{}' - {}", line_num, line, e)
                            })?;
                        relations.push(monomial);
                    }
                    State::Coaction => {
                        let (name, tensors) = line.split_once(':').ok_or(format!(
                            "Line {}: Invalid COACTION format '{}' - expected 'name:tensors'",
                            line_num, line
                        ))?;
                        let tensors: Vec<(F, HashMap<String, usize>, Monomial)> = tensors
                            .split('+')
                            .map::<Result<(F, HashMap<String, usize>, Vec<usize>), String>, _>(|t| {
                                let (s,t) = match t.split_once('.') {
                                    Some((s, t)) => {
                                        (F::parse(s.trim()).map_err(|e| {
                                            format!("Line {}: Scalar {s} with tensor {t} could not be parsed. {e}", line_num)
                                        })?, t)
                                    },
                                    None => {
                                        (F::one(), t)
                                    }
                                };
                                let (l, r) = t
                                    .split_once('|')
                                    .ok_or(format!("Line {}: Invalid tensor format '{}' - expected 'left|right'", line_num, t))?;

                                // Parse left side from coalgebra
                                let coalgebra_elements: HashMap<String, usize> = l.trim().
                                    split(',').map(|e| {
                                        let (name, exponent) = parse_name_exponent(e)?;
                                        if name == "1" {
                                            let none: Result<Option<(String, usize)>, String> = Ok(None);
                                            return none;
                                        }
                                        return Ok(Some((name.to_owned(), exponent)));
                                    }).filter_ok(|e| {
                                        e.is_some()
                                    }).map_ok(|e| {e.unwrap()}).try_collect()?;

                                // TODO: RHS can only be monomial ?
                                let right = parse_monomial(
                                    r.trim(),
                                    &generator_translate,
                                    generators.len(),
                                ).map_err(|e| format!("Line {}: Invalid right monomial '{}' - {}", line_num, r.trim(), e))?;                                
                                Ok((s, coalgebra_elements, right))
                            })
                            .try_collect()?;
                        if name.trim() != generators[coactions.len()].0 {
                            return Err(format!("Line {}: Coaction for '{}' must match generator order, expected '{}'", line_num, name.trim(), generators[coactions.len()].0));
                        }
                        coactions.push(tensors);
                    }
                    _ => return Err(format!("Line {}: Unexpected state", line_num)),
                },
            }
        }

        if state != State::Coaction {
            return Err("Comodule definition is not complete - missing sections".to_owned());
        }

        let n = generators.len();
        let one_monomial: Monomial = vec![0; n];
        let mut monomial_coaction: HashMap<Monomial, Vec<(F, HashMap<String, usize>, Monomial)>> =
            HashMap::default();
        let mut queue: Vec<Monomial> = vec![one_monomial.clone()];

        // Initialize basis information for the unit monomial (1)
        monomial_coaction.insert(
            one_monomial.clone(),
            vec![(F::one(), HashMap::default(), one_monomial.clone())],
        );

        let mut i = 0;
        // BFS loop
        while let Some(current_monomial) = queue.pop() {
            for generator_index in 0..n {
                if let Some(next_monomial) =
                    multiply_monomial_by_generator(&current_monomial, generator_index, &relations)
                {
                    i += 1;
                    if i % 100 == 0 {
                        println!(
                            "Processed {} elements, current basis_information length: {}",
                            i,
                            monomial_coaction.len()
                        );
                    }

                    if !monomial_coaction.contains_key(&next_monomial) {
                        let next_grade = monomial_to_grade(&next_monomial, &generators);
                        if next_grade <= max_grading {
                            // Calculate the coaction for the new monomial
                            let coaction_result = multiply_comodule_coaction_elements::<F, G>(
                                monomial_coaction.get(&current_monomial).ok_or(format!(
                                    "Basis monomial '{}' could not be found in queue",
                                    monomial_to_string(&current_monomial, &generators)
                                ))?,
                                &coactions[generator_index],
                                &relations,
                            );

                            monomial_coaction.insert(next_monomial.clone(), coaction_result);
                            queue.push(next_monomial.clone());
                        }
                    }
                }
            }
        }

        // Construct the basis structure
        let mut basis: HashMap<G, Vec<kBasisElement>> = HashMap::default();
        let mut monomial_to_grade_index: HashMap<Monomial, (G, usize)> =
            HashMap::default();

        for monomial in monomial_coaction.keys().sorted() {
            let grade = monomial_to_grade(monomial, &generators);
            let label = monomial_to_string(monomial, &generators);

            let element = kBasisElement {
                name: label.clone(),
                generator: false,
                primitive: None,
                generated_index: 0,
            };

            let index = basis.entry(grade).or_insert_with(Vec::new).len();
            monomial_to_grade_index.insert(monomial.clone(), (grade, index));
            basis.entry(grade).or_insert_with(Vec::new).push(element);
        }

        // Create the graded vector space
        let comodule_vector_space = GradedVectorSpace::from(basis.clone());

        // create C ⊗ M (coalgebra tensor comodule)
        let tensor = Tensor::generate(&coalgebra.space, &comodule_vector_space);

        // Create the coaction map
        let mut coaction: HashMap<G, M> = HashMap::default();

        for (grade, els) in basis.iter() {
            let tensor_rows = tensor.dimensions[grade];
            let basis_cols = els.len();
            let map = M::zero(basis_cols, tensor_rows);
            coaction.insert(*grade, map);
        }

        println!("{:?}", coalgebra_translate);

        for (monomial, coaction_elements) in &monomial_coaction {
            println!("{:?} :\n {:?} \n\n", monomial, coaction_elements);
            let (basis_grade, basis_index) =
                monomial_to_grade_index.get(monomial).ok_or(format!(
                    "Expected monomial '{}' to exist in lookup",
                    monomial_to_string(monomial, &generators)
                ))?;

            let map = coaction.get_mut(basis_grade).ok_or(format!("Expected a coaction to exist in dimension {basis_grade}. For the element {:?}, {:?}", monomial, coaction_elements))?;

            for (coeff, coalg_indices, mod_monomial) in coaction_elements {
                let mod_grade_index = monomial_to_grade_index.get(mod_monomial).ok_or(format!(
                    "Expected comodule monomial '{}' to exist when constructing coaction",
                    monomial_to_string(mod_monomial, &generators)
                ))?;

                let coalg_name = match coalg_indices.len() {
                    0 => "1",
                    _ => &coalg_indices
                        .iter()
                        .sorted()
                        .filter_map(|(name, ind)| match ind {
                            0 => None,
                            1 => Some(name.to_owned()),
                            _ => Some(name.to_owned() + "^" + &ind.to_string()),
                        })
                        .join(","),
                };

                match coalgebra_translate.get(coalg_name) {
                    Some(coalg_index) => {
                        let (_, tensor_index) = tensor
                            .construct
                            .get(mod_grade_index)
                            .ok_or("Module BasisIndex  not found in tensor construction.")?
                            .get(coalg_index)
                            .ok_or("Algebra BasisIndex not found in tensor construction")?;

                        map.set(*basis_index, *tensor_index, *coeff);
                    }
                    None => {
                        continue;
                    }
                }
            }
        }

        Ok(kComodule::new(
            coalgebra,
            comodule_vector_space,
            GradedLinearMap::from(coaction),
            tensor,
        ))
    }
}

fn multiply_comodule_coaction_elements<F: Field, G: Grading>(
    a: &Vec<(F, HashMap<String, usize>, Monomial)>,
    b: &Vec<(F, HashMap<String, usize>, Monomial)>,
    relations: &Vec<Monomial>,
) -> Vec<(F, HashMap<String, usize>, Monomial)> {
    let mut result: Vec<(F, HashMap<String, usize>, Monomial)> = vec![];

    for x in a {
        for y in b {
            let (a_coeff, a_coalg, a_mod) = x;
            let (b_coeff, b_coalg, b_mod) = y;

            if let Some(mod_product) = multiply_monomials(a_mod, b_mod, relations) {
                let result_coeff = *a_coeff * *b_coeff;

                let mut coalg_product = a_coalg.clone();
                for (name, exponent) in b_coalg {
                    coalg_product
                        .entry(name.to_owned())
                        .and_modify(|e| e.add_assign(exponent))
                        .or_insert(*exponent);
                }

                if let Some(existing) = result
                    .iter_mut()
                    .find(|(_, coalg, monom)| coalg == &coalg_product && monom == &mod_product)
                {
                    existing.0 += result_coeff;
                } else {
                    result.push((result_coeff, coalg_product, mod_product));
                }
            }
        }
    }

    result.retain(|(coeff, _, _)| *coeff != F::zero());
    result
}
