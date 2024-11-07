
use comodule::{Comodule, ComoduleMorphism};
use hopfalgebra::HopfAlgebra;


#[allow(unused)]
mod linearalgebra;
#[allow(unused)]
mod module;
#[allow(unused)]
mod hopfalgebra;
#[allow(unused)]
mod comodule;
#[allow(unused)]
mod resolution;


fn main() {
    println!("Hello, world!");
}

// #[allow(unused)]
// #[allow()]
// fn resolve(hopf: impl HopfAlgebra, comod: impl Comodule) {
//     let zero = comod.get_zero_morphism_to();

//     let mut morphism = zero;
//     loop {
//         let coker_mor = morphism.get_cokernel();
//         let coker = coker_mor.get_domain();
//         let map_to_cogen_mod = coker.get_cogenerating_module();
//         let injection_to_free = 

//     }
// }
