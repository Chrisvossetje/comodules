


pub struct Resolution<M: Morphism, A: CoAlgebra, C: CoModule<A>> {
    comodule: C,
    coalgebra: A,
    resolution: Vec<M>,
}


impl<M, A, C> Resolution<M,A,C> 
where 
A: CoAlgebra, 
C: CoModule<A>,
M: Morphism<C> {
    pub fn new(comodule: C) {

    }

    pub fn resolve(self, s: usize) {
        let zero_morph = self.comodule.zero_morph();
        for i in (0..s) {
            let m = self.resolution.last().unwrap_or(zero_morph);

            
        
        }
    } 
}