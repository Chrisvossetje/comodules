use std::{
    io::{self, Write}, marker::PhantomData, sync::Arc
};

use itertools::Itertools;

use crate::{
    comodule::traits::{Comodule, ComoduleMorphism},
    export::{Page, SSeq}, grading::OrderedGrading,
};

#[derive(Debug, Clone, PartialEq)]
pub struct Resolution<G: OrderedGrading, M: Comodule<G>> {
    pub comodule: Arc<M>,
    pub resolution: Vec<M::Morphism>,
    _grading: PhantomData<G>,
}

impl<G: OrderedGrading, M: Comodule<G>> Resolution<G, M> {
    pub fn new(comodule: M) -> Self {
        Resolution {
            comodule: Arc::new(comodule),
            resolution: vec![],
            _grading: PhantomData,
        }
    }

    pub fn resolve_to_s(&mut self, s: usize, mut limit: G) {
        let mut last_morph: M::Morphism = M::Morphism::zero_morphism(self.comodule.clone());
        
        if self.resolution.len() == 0 {
            let zero_morph = M::Morphism::zero_morphism(self.comodule.clone());
            
            let initial_inject = zero_morph.inject_codomain_to_cofree(limit);
            
            self.resolution.push(M::Morphism::compose(&initial_inject, &zero_morph));
            last_morph = initial_inject;
        }
        

        for _ in self.resolution.len()..=s {
            // Increment limit and get last morphism
            limit = limit.incr();
            
            let coker = &last_morph.cokernel();
            let inject = coker.inject_codomain_to_cofree(limit);
            
            let combine = M::Morphism::compose(&inject, &coker);
            
            last_morph = inject;

            self.resolution.push(combine);
        }
    }

    /// This crashes on WASM
    pub fn resolve_to_s_with_print(&mut self, s: usize, mut limit: G) {
        let mut last_morph: M::Morphism = M::Morphism::zero_morphism(self.comodule.clone());

        
        println!("Resolving to filtration index: {} \n", s);


        if self.resolution.len() == 0 {
            println!("Resolving for 0",);
            print!("Injecting to 0              ",);
            io::stdout().flush().unwrap();
            let inject_time = std::time::Instant::now();

            let zero_morph = M::Morphism::zero_morphism(self.comodule.clone());

            println!("took: {:.2?}\n", inject_time.elapsed());

            let initial_inject = zero_morph.inject_codomain_to_cofree(limit);
        
            self.resolution.push(M::Morphism::compose(&initial_inject, &zero_morph));
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
            let coker = &last_morph.cokernel();
            println!("took: {:.2?}", coker_time.elapsed());

            let inject_time = std::time::Instant::now();
            print!("Injecting to cofree         ");
            io::stdout().flush().unwrap();
            // Inject to cofree
            let inject = coker.inject_codomain_to_cofree(limit);
            println!("took: {:.2?}", inject_time.elapsed());

            let compose_time = std::time::Instant::now();
            print!("Composing morphisms         ");
            // Compose
            io::stdout().flush().unwrap();
            let combine = M::Morphism::compose(&inject, &coker);
            println!("took: {:.2?}", compose_time.elapsed());

            println!(
                "Total time:                 took: {:.2?}\n",
                coker_time.elapsed()
            );

            last_morph = inject;
            self.resolution.push(combine);
        }
    }

    pub fn generate_sseq(&self, name: &str) -> SSeq {
        let (x_formula, y_formula) = G::default_formulas();

        let gens = self
            .resolution
            .iter()
            .enumerate()
            .flat_map(|(s, x)| {
                let g = x.get_codomain().get_generators();
                g.into_iter()
                    .map(move |(id, g, name)| (s, id, g.export_grade(), name))
            })
            .sorted_by_key(|f| (f.0, f.1))
            .collect();

        let lines = self
            .resolution
            .iter()
            .enumerate()
            .flat_map(|(s, x)| {
                let g = x.get_structure_lines();
                g.into_iter()
                    .map(move |(from_gen, to_gen, value, prim_type)| {
                        ((s - 1, from_gen), (s, to_gen), format!("{:?}", value), prim_type)
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
