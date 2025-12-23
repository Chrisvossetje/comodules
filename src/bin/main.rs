use algebra::matrices::f2_matrix::F2Matrix;
use algebra::matrices::flat_matrix::FlatMatrix;
use algebra::rings::finite_fields::F2;
use comodules::export::{Page, SSeq};
use comodules::grading::grading::Grading;
use comodules::grading::unigrading::UniGrading;
use comodules::k_comodule::kcoalgebra::kCoalgebra;
use comodules::resolution::superparallel::ParallelResolution;
use comodules::traits::Coalgebra;
use std::time::Instant;

// fn to_cochain_cpx<F: Field>(res: Resolution<UniGrading, RComodule<UniGrading, F>>) -> (Vec<HashMap<UniGrading, FlatMatrix<UniPolRing<F>>>>, Vec<HashMap<UniGrading, Vec<((kBasisElement, UniGrading, Option<u16>), usize)>>>) {
//     let mut maps = vec![];
//     let mut codomains = vec![];

//     for (_, a) in res.resolution.iter().enumerate() {

//         let mut gr_map = HashMap::default();
//         let mut gr_codom = HashMap::default();

//         for (gr,map) in &a.map.maps {
//             let empty = vec![];
//             let domain = a.domain.as_ref().space.0.get(gr).unwrap_or(&empty);
//             let codomain = a.codomain.as_ref().space.0.get(gr).unwrap_or(&empty);

//             let mut new_domain = vec![];
//             for (id, el) in domain.iter().enumerate() {
//                 if el.0.generator {
//                     new_domain.push((id, el.clone()));
//                 }
//             }

//             let mut new_codomain = vec![];
//             for (id, el) in codomain.iter().enumerate() {
//                 if el.0.generator {
//                     new_codomain.push((id, el.clone()));
//                 }
//             }

//             let mut new_map = FlatMatrix::zero(new_domain.len(), new_codomain.len());

//             for x in 0..new_map.domain() {
//                 for y in 0..new_map.codomain() {
//                     let original_x = new_domain[x].0;
//                     let original_y = new_codomain[y].0;
//                     let val = map.get(original_x, original_y);
//                     new_map.set(x, y, val);
//                 }
//             }

//             gr_map.insert(*gr,new_map);
//             let real_codom = new_codomain.iter().map(|x| {
//                 (x.1.clone(), x.0)
//             }).collect();
//             gr_codom.insert(*gr, real_codom);
//         }
//         maps.push(gr_map);
//         codomains.push(gr_codom);
//     }

//     (maps, codomains)
// }

// fn main() {
//     let start = Instant::now();

//     let input = include_str!("../../examples/polynomial/A_C.txt");
//     let coalgebra = RCoalgebra::<UniGrading, F2>::parse(input, UniGrading(20)).unwrap().0;

//     let coalgebra = Arc::new(coalgebra);

//     let kt = RComodule::fp_comodule(coalgebra, UniGrading(0));

//     let mut res = Resolution::new(kt);

//     res.resolve_to_s_with_print(10, UniGrading(20));

//     let lines: Vec<_> = res.resolution[..(res.resolution.len())-1]
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

//     let ccpx = to_cochain_cpx(res);

//     let mut gens = vec![];

//     let mut map: HashMap<(usize, UniGrading), _> = HashMap::default();

//     for (s, ((f_full, n_full), (g_full, q_full))) in ccpx.0.iter().zip(&ccpx.1).tuple_windows().enumerate() {
//         let mut count = 0;
//         for (gr, f) in f_full {
//             let zero_map = FlatMatrix::zero(0, 0);
//             let g = match g_full.get(gr) {
//                 Some(g) => { g
//                 },
//                 None => {
//                     &zero_map
//                 },
//             };
//             let empty = vec![];

//             let n = n_full.get(gr).unwrap_or(&empty);
//             let q = q_full.get(gr).unwrap_or(&empty);

//             let (cohom_to_n, cohom) = FlatMatrix::cohomology(f, g, &n.iter().map(|x| x.0.2).collect(),  &q.iter().map(|x| x.0.2).collect());

//             let mut cohom_plus = vec![];

//             for b in &cohom {
//                 // TODO: deduce second grading somehow ???
//                 let l = format!("{:?}", b);
//                 // let thing = Some((b.0.generated_index, *k, Some(s)));

//                 cohom_plus.push((b.clone(), count));

//                 gens.push((s, count, gr.export_grade(), Some(l)));
//                 count += 1;
//             }

//             map.insert((s,*gr), (cohom_plus, cohom_to_n));
//         }
//     }

//     println!("{:?}", map);
//     println!("{:?}", lines);

//     let mut new_lines = vec![];

//     // We are gonna do the MOST basic implementation which is probably too coarse.
//     for ((s, gr), (cohom_plus, cohom_to_n)) in &map {
//         let n_full = &ccpx.1[*s][gr];

//         // For each basis in cohom, find a representable vector in n.
//         for (cohom_id, cohom_el) in cohom_plus.iter().enumerate() {
//             let v = cohom_to_n.get_column(cohom_id);

//             // Then we check if this maps to some element in some s+1 thing.
//             for n_id in 0..n_full.len() {
//                 let val = v[n_id];
//                 if !val.is_zero() {
//                     let res_n_id = n_full[n_id].1;
//                     for (source, target, map_value, thing) in &lines {
//                         if &source.0 == s && source.1.2 == res_n_id && source.1.1 == *gr{

//                             if (*s == 3) && (*gr == UniGrading(6)) {
//                                 println!("OPLETTEN!");
//                             }

//                             // Then we check if target is represented by some element in the ker of that s+1 thing
//                             // This check will be super rough as "inclusion" in a vector is enoguh for us.
//                             // So if we map to the second index and the ker is represented by (1 1) then we still map to it
//                             let (target_gr, target_orig_id) = (target.1.1, target.1.2);
//                             let target_full = &ccpx.1[target.0][&target_gr];
//                             for (target_id, target_el) in target_full.iter().enumerate() {
//                                 if target_el.1 == target_orig_id {
//                                     // Now we found an element in the target n
//                                     let (cohom_target_n, target_cohom_to_n) = &map[&(target.0, target_gr)];
//                                     assert_eq!(target_full.len(), target_cohom_to_n.codomain());
//                                     assert_eq!(cohom_target_n.len(), target_cohom_to_n.domain());

//                                     println!("target_cohom_to_n:\n{:?}\n\n cohom_target:{:?}", target_cohom_to_n, cohom_target_n);

//                                     // Here we look for a candidate in the cohom of target n
//                                     for target_cohom_id in 0..cohom_target_n.len() {
//                                         let val2 = target_cohom_to_n.get(target_cohom_id, target_id);
//                                         if !val2.is_zero() {
//                                             let total_val = val * val2 * *map_value;
//                                             if let Some(power) = cohom_target_n[target_cohom_id].0 {
//                                                 if total_val.1 >= power {
//                                                     continue;
//                                                 }
//                                             }

//                                             new_lines.push(((*s, cohom_el.1), (target.0, cohom_target_n[target_cohom_id].1), format!("{:?}", total_val), thing.clone()));
//                                         }

//                                     }

//                                     // Now we have found the element (find_id) in the "n" of the target
//                                     // We will try to find a vector which has this

//                                 }
//                             }

//                         }
//                     }
//                 }
//             }

//             // Then we check if that gets mapped to some element of the cokernel.
//         }
//     }

//     println!("{:?}", new_lines);

//     let page = Page {
//         id: 2,
//         generators: gens,
//         structure_lines: new_lines,
//     };

//     let (x_formula, y_formula) = UniGrading::default_formulas();

//     let sseq = SSeq {
//         name: "A_C".to_owned(),
//         degrees: UniGrading::degree_names(),
//         x_formula,
//         y_formula,
//         pages: vec![page],
//         differentials: vec![],
//     };

//     // let sseq = res.generate_sseq("");

//     sseq.save_to_json("A_C.json").unwrap();

//     println!("\nProgram took: {:.2?}", start.elapsed());
// }

// fn main() {

//     let input = include_str!("../../examples/polynomial/A.txt");
//     let coalgebra = kCoalgebra::<UniGrading, F2, FlatMatrix<F2>>::parse(input, UniGrading(60))
//     .unwrap()
//     .0;

//     println!("Size of coalgebra: {:?}\nSize of space:{:?}\nSize of tensor:{:?}\nSize of coaction:{:?}\n",
//                     deepsize::DeepSizeOf::deep_size_of(&coalgebra),
//                     deepsize::DeepSizeOf::deep_size_of(&coalgebra.space),
//                     deepsize::DeepSizeOf::deep_size_of(&coalgebra.tensor),
//                     deepsize::DeepSizeOf::deep_size_of(&coalgebra.coaction),);
//     let coalgebra = coalgebra;

//     let fp = kComodule::fp_comodule(&coalgebra, UniGrading::zero());

//     let mut res: Resolution<UniGrading, kCoalgebra<UniGrading, F2, FlatMatrix<F2>>, kComoduleMorphism<UniGrading, F2, FlatMatrix<F2>>> =
//         Resolution::new(coalgebra, fp);

//     let start = Instant::now();
//     res.resolve_to_s_with_print(40, UniGrading(60));

//         for s in &res.resolution {

//             println!("Size of map: {:?}\nSize of codomain:{:?}\nSize of space:{:?}\nSize of vec:{:?}\n",         deepsize::DeepSizeOf::deep_size_of(&s.0),
//                     deepsize::DeepSizeOf::deep_size_of(&s.1),
//                     deepsize::DeepSizeOf::deep_size_of(&s.1.space),
//                     deepsize::DeepSizeOf::deep_size_of(&s.1.gen_id_gr),);
//         }

//     println!(
//         "Size of resolution: {:?}",
//         deepsize::DeepSizeOf::deep_size_of(&res)
//     );

//     let sseq = res.generate_sseq("");
//     sseq.save_to_json("A.json").unwrap();

//     println!("\nProgram took: {:.2?}", start.elapsed());
// }

fn main() {
    const MAX_GRADING: UniGrading = UniGrading(90);
    const S: usize = 30;

    let start = Instant::now();

    println!("Started processing coalgebra");
    let input = include_str!("../../examples/polynomial/A.txt");
    let coalgebra = kCoalgebra::<UniGrading, F2, F2Matrix>::parse(input, MAX_GRADING)
        .unwrap()
        .0;
    println!("Ended processing coalgebra");

    println!(
        "Elements in coalgebra: {:?}\nSize of coalgebra: {:?}\nSize of space:{:?}\nSize of coaction:{:?}\n",
        coalgebra
            .space
            .0
            .iter()
            .fold(0, |count, g| count + g.1.len()),
        deepsize::DeepSizeOf::deep_size_of(&coalgebra),
        deepsize::DeepSizeOf::deep_size_of(&coalgebra.space),
        deepsize::DeepSizeOf::deep_size_of(&coalgebra.coaction),
    );
    let coalgebra = coalgebra;

    let fp = coalgebra.basering_comodule(UniGrading::zero());
    let res = ParallelResolution::init(coalgebra, fp, S, MAX_GRADING);
    res.populate_with_basering();

    let resolution_time = Instant::now();

    // // Parallal Executor
    // rayon::scope(|i| {
    //     res.recursion_solve(i, 1, UniGrading::zero());
    // });


    // Non Parallal executor
    for g in MAX_GRADING.iterator_from_zero(true) {
        for s in 1..=S {
            if s==2 && g == UniGrading(5) {
                println!("S:{}, G:{}", s, g);
            }
            res.resolve_at_s_g(s, g);
        }
    }


    for s in 0..=S {
        for g in MAX_GRADING.iterator_from_zero(true) {
            debug_assert!(&res.data[s][g.to_index()].1.get().is_some());
        }
    }

    let mut gens = vec![];

    for s in 0..=S {
        for g in MAX_GRADING.iterator_from_zero(true) {
            let mut a_gens = res
                .get_data_cell(s, g)
                .a_gens
                .iter()
                .map(|x| (s, x.1 as usize, g.export_grade(), None))
                .collect();
            gens.append(&mut a_gens);
        }
    }

    let page = Page {
        id: 2,
        generators: gens,
        structure_lines: vec![],
    };

    let (x_formula, y_formula) = UniGrading::default_formulas();

    let sseq = SSeq {
        name: "A".to_owned(),
        degrees: UniGrading::degree_names(),
        x_formula,
        y_formula,
        pages: vec![page],
        differentials: vec![],
    };

    sseq.save_to_json("A_par.json").unwrap();

    println!(
        "\nCoalgebra generation took: {:.2?}",
        (resolution_time - start)
    );
    println!("\nResolution took: {:.2?}", resolution_time.elapsed());

    println!("\nProgram took: {:.2?}", start.elapsed());
}
