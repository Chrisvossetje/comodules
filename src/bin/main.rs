use std::sync::Arc;

use std::time::Instant;

use comodules::{
    comodule::{kcoalgebra::kCoalgebra, kcomodule::kComodule, traits::Comodule},
    grading::UniGrading,
    linalg::{field::F2, flat_matrix::FlatMatrix},
    resolution::Resolution,
};

fn main() {
    let start = Instant::now();

    let input = include_str!("../../examples/polynomial/A(2).txt");
    const MAX_GRADING: UniGrading = UniGrading(20);
    let (coalgebra, _translate) = kCoalgebra::parse(input, MAX_GRADING).unwrap();
    // let coalgebra = Arc::new(kCoalgebra::parse(input, MAX_GRADING).unwrap().0);

    let coalgebra = Arc::new(coalgebra);
    
    // let input = include_str!("../../examples/comodule/A(1).txt");
    // let comod = kComodule::parse(input, coalgebra, &translate, MAX_GRADING).unwrap();
    
    let comod = kComodule::fp_comodule(coalgebra, UniGrading(1));
    let mut res: Resolution<UniGrading, kComodule<UniGrading, F2, FlatMatrix<F2>>> =
        Resolution::new(comod);

    res.resolve_to_s_with_print(8, MAX_GRADING);

    let _page = res.generate_sseq("?");

    let _ = _page.save_to_json("page.json");

    println!("\nProgram took: {:.2?}", start.elapsed());
}
