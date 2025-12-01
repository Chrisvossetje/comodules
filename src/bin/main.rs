use std::sync::Arc;

use std::time::Instant;
use ahash::HashMap;
use comodules::basiselement::kBasisElement;
use comodules::comodule::traits::ComoduleMorphism;
use comodules::grading::Grading;
use comodules::linalg::field::Field;
use comodules::linalg::matrix::RModMorphism;
use comodules::linalg::ring::UniPolRing;
use comodules::module::module::cohomology;
use comodules::{
    comodule::{kcoalgebra::kCoalgebra, kcomodule::kComodule, rcomodule::{RCoalgebra, RComodule}, traits::Comodule}, export::{Page, SSeq}, grading::UniGrading, linalg::{field::F2, flat_matrix::FlatMatrix}, resolution::Resolution
};
use itertools::Itertools;


fn to_cochain_cpx<F: Field>(res: Resolution<UniGrading, RComodule<UniGrading, F>>) -> (Vec<HashMap<UniGrading, FlatMatrix<UniPolRing<F>>>>, Vec<HashMap<UniGrading, Vec<(kBasisElement, UniGrading, Option<u16>)>>>) {
    let mut maps = vec![];
    let mut codomains = vec![];

    // res.resolution[0].domain;
    for (id,a) in res.resolution.iter().enumerate() {

        let mut gr_map = HashMap::default();
        let mut gr_codom = HashMap::default();

        for (gr,map) in &a.map.maps {
            let empty = vec![];
            let domain = a.domain.as_ref().space.0.get(gr).unwrap_or(&empty);
            let codomain = a.codomain.as_ref().space.0.get(gr).unwrap_or(&empty);

            let mut new_domain = vec![];
            for (id, el) in domain.iter().enumerate() {
                if el.0.generator {
                    new_domain.push((id, el.clone()));
                }
            }
            
            let mut new_codomain = vec![];
            for (id, el) in codomain.iter().enumerate() {
                if el.0.generator {
                    new_codomain.push((id, el.clone()));
                }
            }
            

            let mut new_map = FlatMatrix::zero(new_domain.len(), new_codomain.len());

            for x in 0..new_map.domain {
                for y in 0..new_map.codomain {
                    let original_x = new_domain[x].0;
                    let original_y = new_codomain[y].0;
                    let val = map.get(original_x, original_y);
                    new_map.set(x, y, val);
                }
            }

            gr_map.insert(*gr,new_map);
            let real_codom = new_codomain.iter().map(|x| {
                x.1.clone()
            }).collect();
            gr_codom.insert(*gr, real_codom);
        }
        maps.push(gr_map);
        codomains.push(gr_codom);
    }

    (maps, codomains)
}

// fn to_maps<F: Field>(lines: )

fn main() {
    let start = Instant::now();

    let input = include_str!("../../examples/polynomial/A_C.txt");
    let coalgebra = RCoalgebra::<UniGrading, F2>::parse(input, UniGrading(10)).unwrap().0; 

    let coalgebra = Arc::new(coalgebra);

    let kt = RComodule::fp_comodule(coalgebra, UniGrading(0));

    let mut res = Resolution::new(kt);

    res.resolve_to_s_with_print(5, UniGrading(10));

    let lines: Vec<_> = res.resolution
        .iter()
        .enumerate()
        .flat_map(|(s, x)| {
            let g = x.get_structure_lines();
            g.into_iter()
                .map(move |(from_gen, to_gen, value, prim_type)| {
                    ((s - 1, from_gen), (s, to_gen), value, prim_type)
                })
        })
        .sorted_by_key(|f| (f.0, f.1))
        .collect();
    
    let ccpx = to_cochain_cpx(res);
    

    let mut gens = vec![];
    
    // let map = 


    for (s, ((f_full, n_full), (g_full, q_full))) in ccpx.0.iter().zip(ccpx.1).tuple_windows().enumerate() {
        let mut count = 0;
        for (gr, f) in f_full {
            let zero_map = FlatMatrix::zero(0, 0);
            let g = match g_full.get(gr) {
                Some(g) => { g
                },
                None => {
                    &zero_map
                },
            };
            let empty = vec![];
            
            let n = n_full.get(gr).unwrap_or(&empty);
            let q = q_full.get(gr).unwrap_or(&empty); 

            let (cohom, original) = cohomology(f, g, n, q);



            for b in cohom {
                let l = format!("{:?} | {:?}", b.2, b.1);
                // let thing = Some((b.0.generated_index, *k, Some(s)));

                gens.push((s, count, gr.export_grade(), Some(l)));
                count += 1;
            }
        }            
    }

    let page = Page {
        id: 2,
        generators: gens,
        structure_lines: vec![],
    };

    let (x_formula, y_formula) = UniGrading::default_formulas();

    let sseq = SSeq {
        name: "A_C".to_owned(),
        degrees: UniGrading::degree_names(),
        x_formula,
        y_formula,
        pages: vec![page],
        differentials: vec![],
    };


    // let sseq = res.generate_sseq("");

    sseq.save_to_json("A_C.json").unwrap();

    println!("\nProgram took: {:.2?}", start.elapsed());
}
