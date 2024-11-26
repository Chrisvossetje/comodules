




// pub struct Resolution<M: Morphism, A: CoAlgebra, C: CoModule<A>> {
//     comodule: C,
//     coalgebra: A,
//     resolution: Vec<M>,
// }

use std::marker::PhantomData;

use crate::{coalgebra::Coalgebra, comodule::{Comodule, ComoduleMorphism}, field::Field, graded::Grading, matrix::Matrix};
pub struct Resolution<G: Grading, F: Field, M: Matrix<F>, A: Coalgebra<G, F, M>, C: Comodule<G, F, M>, CM: ComoduleMorphism<G, F, M>> {
    comodule: C,
    coalgebra: A,
    resolution: Vec<CM>,

    grading: PhantomData<G>,
    field: PhantomData<F>,
    matrix: PhantomData<M>,
}


impl <G: Grading, F: Field, M: Matrix<F>, A: Coalgebra<G, F, M>, C: Comodule<G, F, M>, CM: ComoduleMorphism<G, F, M>> Resolution<G, F, M, A, C, CM>{
    pub fn new(comodule: C) {
        let coalgebra = comodule.get_coalgebra();
        Resolution {
            comodule,
            coalgebra,
            resolution: vec![],
        }
    }

    pub fn resolve(self, s: usize) {
        let zero_morph = self.comodule.zero_morph();
        for i in (0..s) {
            let m = self.resolution.last().unwrap_or(zero_morph);

            
        
        }
    } 
}