#[cfg(test)]
mod tests {
    use ahash::HashMap;
    use algebra::{
        abelian::Abelian,
        field::Field,
        matrices::flat_matrix::FlatMatrix,
        matrix::Matrix,
        rings::{finite_fields::F2, univariate_polynomial_ring::UniPolRing},
    };
    use comodules::{
        export::{Page, SSeq},
        grading::{grading::Grading, unigrading::UniGrading},
        k_comodule::kcoalgebra::{A0_coalgebra, kCoalgebra},
        k_t_comodule::k_t_coalgebra::{A0_C, ktCoalgebra, tensor_k_coalgebra},
        resolution::resolve_by_s::Resolution,
        traits::Coalgebra,
    };
    use itertools::Itertools;

    #[test]
    fn test_a0_c_resolution() {
        let coalgebra = A0_C();

        let kt = ktCoalgebra::basering_comodule(&coalgebra, UniGrading(0));

        let mut res = Resolution::new(coalgebra, kt);

        res.resolve_to_s(20, UniGrading(40));

        let sseq = res.generate_sseq("");

        sseq.save_to_json("LOL.json").unwrap();
    }

    #[test]
    fn test_a0_tensor_resolution() {
        let coalgebra = A0_coalgebra();

        let tensor_coalgebra = tensor_k_coalgebra(coalgebra);

        let kt = ktCoalgebra::basering_comodule(&tensor_coalgebra, UniGrading(0));

        let mut res = Resolution::new(tensor_coalgebra, kt);

        res.resolve_to_s(20, UniGrading(40));

        let _sseq = res.generate_sseq("A0");
        let _comp_sseq: SSeq = serde_json::from_str(include_str!("./A(0).json")).unwrap();

        // TODO :
        // assert_eq!(sseq.pages[0].generators, comp_sseq.pages[0].generators);
    }

    #[test]
    fn test_a1_tensor_resolution() {
        let input = include_str!("../examples/direct/A(1).txt");
        let coalgebra = kCoalgebra::parse(input, UniGrading::infty()).unwrap().0;

        let tensor_coalgebra = tensor_k_coalgebra(coalgebra);

        let kt = ktCoalgebra::basering_comodule(&tensor_coalgebra, UniGrading(0));

        let mut res = Resolution::new(tensor_coalgebra, kt);

        res.resolve_to_s(20, UniGrading(40));

        let sseq = res.generate_sseq("");

        sseq.save_to_json("LOL.json").unwrap();
    }

    #[test]
    fn test_a1_tensor_parser_resolution() {
        let input = include_str!("../examples/direct/A(1)_dual_graded.txt");
        let coalgebra = ktCoalgebra::<UniGrading, F2, FlatMatrix<UniPolRing<F2>>>::parse(
            input,
            UniGrading::infty(),
        )
        .unwrap()
        .0;

        let kt = ktCoalgebra::basering_comodule(&coalgebra, UniGrading(0));

        let mut res = Resolution::new(coalgebra, kt);

        res.resolve_to_s(20, UniGrading(40));

        let sseq = res.generate_sseq("");

        sseq.save_to_json("LOL.json").unwrap();
    }

    #[test]
    fn test_a_tensor_resolution() {
        let input = include_str!("../examples/polynomial/A.txt");
        let coalgebra = kCoalgebra::parse(input, UniGrading(20)).unwrap().0;

        let tensor_coalgebra = tensor_k_coalgebra(coalgebra);

        let kt = ktCoalgebra::basering_comodule(&tensor_coalgebra, UniGrading(0));

        let mut res = Resolution::new(tensor_coalgebra, kt);

        res.resolve_to_s(20, UniGrading(15));

        let sseq = res.generate_sseq("");

        sseq.save_to_json("LOL.json").unwrap();
    }

    #[test]
    fn test_a1_c_resolution() {
        let input = include_str!("../examples/direct/A(1)_C.txt");
        let coalgebra = ktCoalgebra::<UniGrading, F2, FlatMatrix<UniPolRing<F2>>>::parse(
            input,
            UniGrading::infty(),
        )
        .unwrap()
        .0;

        let kt = ktCoalgebra::basering_comodule(&coalgebra, UniGrading(0));

        let mut res = Resolution::new(coalgebra, kt);

        res.resolve_to_s(25, UniGrading(80));

        let ccpx = to_cochain_cpx(res);

        let mut gens = vec![];

        for (s, ((f_full, n_full), (g_full, q_full))) in
            ccpx.0.iter().zip(ccpx.1).tuple_windows().enumerate()
        {
            let mut count = 0;
            for (gr, f) in f_full {
                let zero_map = FlatMatrix::zero(0, 0);
                let g = match g_full.get(gr) {
                    Some(g) => g,
                    None => &zero_map,
                };
                let empty = vec![];

                let n = n_full.get(gr).unwrap_or(&empty);
                let q = q_full.get(gr).unwrap_or(&empty);

                let (_, cohom) = FlatMatrix::cohomology(f, g, &n, &q);

                for b in cohom {
                    let l = format!("{:?}", b);

                    println!("{:?}", l);
                    gens.push((s, count, gr.export_grade(), Some(l)));
                    count += 1;
                }
            }
        }

        println!("{:?}", gens);

        let page = Page {
            id: 2,
            generators: gens,
            structure_lines: vec![],
        };

        let (x_formula, y_formula) = UniGrading::default_formulas();

        let sseq = SSeq {
            name: "A(1)_C".to_owned(),
            degrees: UniGrading::degree_names(),
            x_formula,
            y_formula,
            pages: vec![page],
            differentials: vec![],
        };

        // let sseq = res.generate_sseq("");

        sseq.save_to_json("A(1)_C.json").unwrap();
    }

    #[test]
    fn test_a1_c_parser() {
        let input = include_str!("../examples/polynomial/A(1)_C.txt");
        let coalgebra =
            ktCoalgebra::<UniGrading, F2, FlatMatrix<UniPolRing<F2>>>::parse(input, UniGrading(80))
                .unwrap()
                .0;

        let kt = ktCoalgebra::basering_comodule(&coalgebra, UniGrading(0));

        let mut res = Resolution::new(coalgebra, kt);

        res.resolve_to_s(25, UniGrading(80));

        let ccpx = to_cochain_cpx(res);

        let mut gens = vec![];

        for (s, ((f_full, n_full), (g_full, q_full))) in
            ccpx.0.iter().zip(ccpx.1).tuple_windows().enumerate()
        {
            let mut count = 0;
            for (gr, f) in f_full {
                let zero_map = FlatMatrix::zero(0, 0);
                let g = match g_full.get(gr) {
                    Some(g) => g,
                    None => &zero_map,
                };
                let empty = vec![];

                let n = n_full.get(gr).unwrap_or(&empty);
                let q = q_full.get(gr).unwrap_or(&empty);

                let (_, cohom) = FlatMatrix::cohomology(f, g, &n, &q);

                for b in cohom {
                    let l = format!("{:?}", b);

                    println!("{:?}", l);
                    gens.push((s, count, gr.export_grade(), Some(l)));
                    count += 1;
                }
            }
        }

        println!("{:?}", gens);

        let page = Page {
            id: 2,
            generators: gens,
            structure_lines: vec![],
        };

        let (x_formula, y_formula) = UniGrading::default_formulas();

        let sseq = SSeq {
            name: "A(1)_C".to_owned(),
            degrees: UniGrading::degree_names(),
            x_formula,
            y_formula,
            pages: vec![page],
            differentials: vec![],
        };

        // let sseq = res.generate_sseq("");

        sseq.save_to_json("A(1)_C.json").unwrap();
    }

    #[test]
    fn test_a2_c_parser() {
        let input = include_str!("../examples/polynomial/A(2)_C.txt");
        let coalgebra =
            ktCoalgebra::<UniGrading, F2, FlatMatrix<UniPolRing<F2>>>::parse(input, UniGrading(30))
                .unwrap()
                .0;

        let kt = ktCoalgebra::basering_comodule(&coalgebra, UniGrading(0));

        let mut res = Resolution::new(coalgebra, kt);

        res.resolve_to_s_with_print(10, UniGrading(30));

        let ccpx = to_cochain_cpx(res);

        let mut gens = vec![];

        for (s, ((f_full, n_full), (g_full, q_full))) in
            ccpx.0.iter().zip(ccpx.1).tuple_windows().enumerate()
        {
            let mut count = 0;
            for (gr, f) in f_full {
                let zero_map = FlatMatrix::zero(0, 0);
                let g = match g_full.get(gr) {
                    Some(g) => g,
                    None => &zero_map,
                };
                let empty = vec![];

                let n = n_full.get(gr).unwrap_or(&empty);
                let q = q_full.get(gr).unwrap_or(&empty);

                let (_, cohom) = FlatMatrix::cohomology(f, g, &n, &q);

                for b in cohom {
                    let l = format!("{:?}", b);

                    if b.1.is_some() {
                        continue;
                    }
                    println!("{:?}", l);
                    gens.push((s, count, gr.export_grade(), Some(l)));
                    count += 1;
                }
            }
        }

        println!("{:?}", gens);

        let page = Page {
            id: 2,
            generators: gens,
            structure_lines: vec![],
        };

        let (x_formula, y_formula) = UniGrading::default_formulas();

        let sseq = SSeq {
            name: "A(2)_C".to_owned(),
            degrees: UniGrading::degree_names(),
            x_formula,
            y_formula,
            pages: vec![page],
            differentials: vec![],
        };

        // let sseq = res.generate_sseq("");

        sseq.save_to_json("A(2)_C.json").unwrap();
    }

    // fn generate_sseq(lol: ?, name: &str) -> SSeq {
    //     let (x_formula, y_formula) = G::default_formulas();

    //     let gens = self
    //         .resolution
    //         .iter()
    //         .enumerate()
    //         .flat_map(|(s, x)| {
    //             let g = x.get_codomain().get_generators();
    //             g.into_iter()
    //                 .map(move |(id, g, name)| (s, id, g.export_grade(), name))
    //         })
    //         .sorted_by_key(|f| (f.0, f.1))
    //         .collect();

    //     let lines = self
    //         .resolution
    //         .iter()
    //         .enumerate()
    //         .flat_map(|(s, x)| {
    //             let g = x.get_structure_lines();
    //             g.into_iter()
    //                 .map(move |(from_gen, to_gen, value, prim_type)| {
    //                     ((s - 1, from_gen), (s, to_gen), value, prim_type)
    //                 })
    //         })
    //         .sorted_by_key(|f| (f.0, f.1))
    //         .collect();

    //     let page = Page {
    //         id: 2,
    //         generators: gens,
    //         structure_lines: lines,
    //     };

    //     SSeq {
    //         name: name.to_owned(),
    //         degrees: G::degree_names(),
    //         x_formula,
    //         y_formula,
    //         pages: vec![page],
    //         differentials: vec![],
    //     }
    // }

    fn to_cochain_cpx<F: Field, M: Abelian<UniPolRing<F>>>(
        res: Resolution<UniGrading, ktCoalgebra<UniGrading, F, M>>,
    ) -> (
        Vec<HashMap<UniGrading, FlatMatrix<UniPolRing<F>>>>,
        Vec<HashMap<UniGrading, Vec<M::Generator>>>,
    ) {
        let mut maps = vec![];
        let mut codomains = vec![];

        // res.resolution[0].domain;
        for (s, a) in res.resolution.iter().enumerate() {
            let mut gr_map = HashMap::default();
            let mut gr_codom = HashMap::default();

            for (gr, map) in &a.0.map.maps {
                let empty = vec![];
                let domain = if s == 0 {
                    &empty
                } else {
                    res.resolution[s - 1].1.space.0.get(gr).unwrap_or(&empty)
                };
                let codomain = a.1.space.0.get(gr).unwrap_or(&empty);

                let mut new_domain = vec![];
                for (id, el) in domain.iter().enumerate() {
                    if el.0.0.0 == UniGrading::zero() {
                        new_domain.push((id, el.clone()));
                    }
                }

                let mut new_codomain = vec![];
                for (id, el) in codomain.iter().enumerate() {
                    if el.0.0.0 == UniGrading::zero() {
                        new_codomain.push((id, el.clone()));
                    }
                }

                let mut new_map = FlatMatrix::zero(new_domain.len(), new_codomain.len());

                for x in 0..new_map.domain() {
                    for y in 0..new_map.codomain() {
                        let original_x = new_domain[x].0;
                        let original_y = new_codomain[y].0;
                        let val = map.get(original_x, original_y);
                        new_map.set(x, y, val);
                    }
                }

                gr_map.insert(*gr, new_map);
                let real_codom = new_codomain.iter().map(|x| x.1.1.clone()).collect();
                gr_codom.insert(*gr, real_codom);
            }
            maps.push(gr_map);
            codomains.push(gr_codom);
        }

        (maps, codomains)
    }
}
