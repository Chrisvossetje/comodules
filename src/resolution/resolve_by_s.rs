use std::{
    io::{self, Write},
    marker::PhantomData,
};

use deepsize::DeepSizeOf;
use itertools::Itertools;

use crate::{
    export::{Page, SSeq},
    grading::grading::Grading,
    traits::*,
};

#[derive(Debug, Clone, PartialEq, DeepSizeOf)]
pub struct Resolution<G: Grading, C: Coalgebra<G>> {
    pub coalgebra: C,
    pub comodule: C::Comod,
    pub resolution: Vec<(C::ComodMorph, C::CofMod)>,
    _grading: PhantomData<G>,
}

impl<G: Grading, C: Coalgebra<G>> Resolution<G, C> {
    pub fn new(coalgebra: C, comodule: C::Comod) -> Self {
        Resolution {
            coalgebra,
            comodule,
            resolution: vec![],
            _grading: PhantomData,
        }
    }

    pub fn resolve_to_s(&mut self, s: usize, mut limit: G) {
        let mut last_morph = C::ComodMorph::zero_morphism(&self.comodule);

        if self.resolution.len() == 0 {
            let zero_morph = C::ComodMorph::zero_morphism(&self.comodule);

            let (initial_inject, a_0) =
                C::ComodMorph::inject_codomain_to_cofree(&self.coalgebra, &self.comodule, limit);

            let compose = C::ComodMorph::compose(&initial_inject, &zero_morph, &a_0);

            self.resolution.push((compose, a_0));
            last_morph = initial_inject;
        }

        for _ in self.resolution.len()..=s {
            // Increment limit and get last morphism
            limit = limit.incr();

            let a_i = &self.resolution.last().unwrap().1;
            let (coker_map, coker_comodule) = last_morph.cokernel(&self.coalgebra, &a_i);
            let (inject, a_i) =
                C::ComodMorph::inject_codomain_to_cofree(&self.coalgebra, &coker_comodule, limit);

            let compose = C::ComodMorph::compose(&inject, &coker_map, &a_i);

            last_morph = inject;

            self.resolution.push((compose, a_i));
        }
    }

    /// This crashes on WASM
    pub fn resolve_to_s_with_print(&mut self, s: usize, mut limit: G) {
        let mut last_morph = C::ComodMorph::zero_morphism(&self.comodule);

        println!("Resolving to filtration index: {} \n", s);

        if self.resolution.len() == 0 {
            println!("Resolving for 0",);
            print!("Injecting to 0              ",);
            io::stdout().flush().unwrap();
            let inject_time = std::time::Instant::now();

            let zero_morph = C::ComodMorph::zero_morphism(&self.comodule);

            let (initial_inject, a_0) =
                C::ComodMorph::inject_codomain_to_cofree(&self.coalgebra, &self.comodule, limit);
            let compose = C::ComodMorph::compose(&initial_inject, &zero_morph, &a_0);

            println!("took: {:.2?}\n", inject_time.elapsed());

            self.resolution.push((compose, a_0));
            last_morph = initial_inject;
        }

        for i in self.resolution.len()..=s {
            // Increment limit and get last morphism
            limit = limit.incr();

            println!("Resolving for {}", i);
            let coker_time = std::time::Instant::now();
            print!("Finding cokernel            ");
            // Cokernel
            io::stdout().flush().unwrap();
            let a_i = &self.resolution.last().unwrap().1;
            let (coker_map, coker_comodule) = last_morph.cokernel(&self.coalgebra, &a_i);
            println!("took: {:.2?}", coker_time.elapsed());

            let inject_time = std::time::Instant::now();
            print!("Injecting to cofree         ");
            io::stdout().flush().unwrap();
            // Inject to cofree
            let (inject, a_i) =
                C::ComodMorph::inject_codomain_to_cofree(&self.coalgebra, &coker_comodule, limit);
            println!("took: {:.2?}", inject_time.elapsed());

            let compose_time = std::time::Instant::now();
            print!("Composing morphisms         ");
            // Compose
            io::stdout().flush().unwrap();
            let compose = C::ComodMorph::compose(&inject, &coker_map, &a_i);
            println!("took: {:.2?}", compose_time.elapsed());

            println!(
                "Total time:                 took: {:.2?}\n",
                coker_time.elapsed()
            );

            last_morph = inject;

            self.resolution.push((compose, a_i));
        }
    }

    pub fn generate_sseq(&self, name: &str) -> SSeq {
        let (x_formula, y_formula) = G::default_formulas();

        let gens = self
            .resolution
            .iter()
            .enumerate()
            .flat_map(|(s, x)| {
                let g = x.1.get_generators();
                g.into_iter()
                    .map(move |(id, g, name)| (s, id, g.export_grade(), name))
            })
            .sorted_by_key(|f| (f.0, f.1))
            .collect();

        let lines = self
            .resolution
            .iter()
            .skip(1)
            .enumerate()
            .flat_map(|(s, x)| {
                let s = s + 1; // Note Skip above
                let domain = &self.resolution[s - 1].1;
                let g = x.0.get_structure_lines(&self.coalgebra, domain, &x.1);
                g.into_iter()
                    .map(move |(from_gen, to_gen, value, prim_type)| {
                        (
                            (s - 1, from_gen.0),
                            (s, to_gen.0),
                            format!("{:?}", value),
                            prim_type,
                        )
                    })
            })
            .sorted_by_key(|f| (f.0, f.1))
            .collect();

        let page = Page {
            id: 2,
            generators: gens,
            structure_lines: lines,
        };

        SSeq {
            name: name.to_owned(),
            degrees: G::degree_names(),
            x_formula,
            y_formula,
            pages: vec![page],
            differentials: vec![],
        }
    }
}
