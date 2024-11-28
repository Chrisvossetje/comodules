use crate::comodule::{Comodule, ComoduleMorphism};


pub struct Resolution<M: Comodule, F: ComoduleMorphism<M>> {
    comodule: M,
    resolution: Vec<F>,
}


impl<M: Comodule, F: ComoduleMorphism<M>> Resolution<M, F>{
    pub fn new(comodule: M) -> Self {
        Resolution {
            comodule,
            resolution: vec![],
        }
    }

    pub fn resolve_to_s(mut self, s: usize) {
        let zero_morph = F::zero_morphism(self.comodule);

        let initial_inject = zero_morph.injection_codomain_to_cofree();
        self.resolution.push(initial_inject);
        
        for i in (0..s) {
            let last_morph = self.resolution.last().unwrap();

            let coker = last_morph.cokernel();

            let inject = last_morph.injection_codomain_to_cofree();

            let combine = F::compose(coker, inject);

            self.resolution.push(combine);
        }
    } 
}