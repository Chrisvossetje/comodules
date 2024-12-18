use std::{marker::PhantomData, sync::Arc};

use crate::{
    comodule::comodule::{Comodule, ComoduleMorphism},
    linalg::graded::Grading,
    page::Page,
};

#[derive(Debug, Clone, PartialEq)]
pub struct Resolution<G: Grading, M: Comodule<G>, Morph: ComoduleMorphism<G, M>> {
    comodule: Arc<M>,
    resolution: Vec<Morph>,
    _grading: PhantomData<G>,
}

impl<G: Grading, M: Comodule<G>, Morph: ComoduleMorphism<G, M>> Resolution<G, M, Morph> {
    pub fn new(comodule: M) -> Self {
        Resolution {
            comodule: Arc::new(comodule),
            resolution: vec![],
            _grading: PhantomData,
        }
    }

    pub fn resolve_to_s(&mut self, s: usize, mut limit: G) {
        let zero_morph = Morph::zero_morphism(self.comodule.clone());
        let fixed_limit = limit.incr().incr();

        let initial_inject = zero_morph.inject_codomain_to_cofree(limit, fixed_limit);
        self.resolution.push(initial_inject);

        for _ in 0..s {
            limit = limit.incr();
            let fixed_limit = limit.incr().incr();
            let last_morph = self.resolution.last().unwrap();

            let coker = last_morph.cokernel();

            let inject = coker.inject_codomain_to_cofree(limit, fixed_limit);

            let combine = Morph::compose(inject, coker);

            self.resolution.push(combine);
        }
    }

    pub fn generate_page(&self) -> Page {
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
            .collect();

        let lines = self
            .resolution
            .iter()
            .enumerate()
            .flat_map(|(s, x)| {
                let g = x.get_structure_lines();
                g.into_iter()
                    .map(move |(from_gen, to_gen, value, prim_type)| {
                        ((s - 1, from_gen), (s, to_gen), value, prim_type)
                    })
            })
            .collect();

        // TODO: Structure Lines
        Page {
            name: " ?? ".to_string(),
            id: 2,
            degrees: G::degree_names(),
            x_formula,
            y_formula,
            generators: gens,
            structure_lines: lines,
            differentials: vec![],
        }
    }
}
