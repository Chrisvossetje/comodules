use std::marker::PhantomData;

use crate::{comodule::{Comodule, ComoduleMorphism}, graded::Grading, page::Page, utils::RefType};


pub struct Resolution<G: Grading, M: Comodule<G>, Morph: ComoduleMorphism<G, M>> {
    comodule: M,
    resolution: Vec<Morph>,
    _grading: PhantomData<G>,
}


impl<G: Grading, M: Comodule<G>, Morph: ComoduleMorphism<G, M>> Resolution<G, M, Morph>{
    pub fn new(comodule: M) -> Self {
        Resolution {
            comodule,
            resolution: vec![],
            _grading: PhantomData,
        }
    }

    pub fn resolve_to_s(mut self, s: usize) {
        let comod = RefType::new(self.comodule);
        let zero_morph = Morph::zero_morphism(comod);

        let initial_inject = zero_morph.inject_codomain_to_cofree();
        self.resolution.push(initial_inject);
        
        for i in (0..s) {
            let last_morph = self.resolution.last().unwrap();

            let coker = last_morph.cokernel();

            let inject = last_morph.inject_codomain_to_cofree();

            let combine = Morph::compose(coker, inject);

            self.resolution.push(combine);
        }
    } 

    pub fn generate_page(&self) -> Page {
        let (x_formula, y_formula) = G::default_formulas();
        
        let gens = self.resolution.iter().enumerate().flat_map(|(s, x)| {
            let g = x.get_codomain().get_generators();
            g.into_iter().map(move |(id, g, name)| { 
                (s, id, g.export_grade(), name)
            })
        }).collect();

        Page {
            name: " ?? ".to_string(),
            id: 2,
            degrees: G::degree_names(),
            x_formula,
            y_formula,
            generators: gens,
            structure_lines: todo!(),
            differentials: vec![],
        }
    }
}