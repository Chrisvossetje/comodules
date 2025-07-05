use std::sync::Arc;

use std::time::Instant;

use comodules::{
    comodule::{kcoalgebra::kCoalgebra, kcomodule::kComodule},
    linalg::{field::F2, flat_matrix::FlatMatrix, grading::UniGrading},
    resolution::Resolution,
};

fn main() {
    let start = Instant::now();

    let input = include_str!("../../examples/polynomial/A(1).txt");
    const MAX_GRADING: i32 = 60;
    let (coalgebra, translate) = kCoalgebra::parse(input, MAX_GRADING).unwrap();
    // let coalgebra = Arc::new(kCoalgebra::parse(input, MAX_GRADING).unwrap().0);

    let coalgebra = Arc::new(coalgebra);
    let input = include_str!("../../examples/comodule/A(1).txt");

    let comod = kComodule::parse(input, coalgebra, &translate, MAX_GRADING).unwrap();
    let mut res: Resolution<UniGrading, kComodule<UniGrading, F2, FlatMatrix<F2>>> =
        Resolution::new(comod);

    res.resolve_to_s_with_print(3, MAX_GRADING);

    let page = res.generate_sseq("?");

    let _ = page.save_to_json("page.json");

    println!("\nProgram took: {:.2?}", start.elapsed());
}
