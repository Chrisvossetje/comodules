use algebra::abelian::Abelian;
use deepsize::DeepSizeOf;

use crate::{grading::grading::Grading, resolution::datacell::{DataCell, Module}, traits::Coalgebra, types::ComoduleIndexType};

/// Everything below is wrt. some degree g
///
/// to_cokernel: morphism from A\otimes V_i-1 to Q_i
///
/// cokernel: module structure of the cokernel Q_i
/// 
/// repr_vecs: inverse image of some mapping
/// 
#[derive(Debug, DeepSizeOf)]
pub struct CokerCell<G: Grading, C: Coalgebra<G>> {
    pub to_cokernel: C::RingMorph,
    pub cokernel: Vec<<C::RingMorph as Abelian<C::BaseRing>>::Generator>,
    pub repr_vecs: C::RingMorph,
}


impl<G: Grading, C: Coalgebra<G>> CokerCell<G, C> {
    pub fn cokernel(prev_s: &DataCell<G, C>) -> CokerCell<G, C> {
        let codom_module = prev_s.r_gens.iter().map(|x| x.1).collect();
        let (to_cokernel, repr_vecs, cokernel) = prev_s.to_cofree.cokernel(&codom_module);
        CokerCell { to_cokernel, cokernel, repr_vecs }
    }
}
