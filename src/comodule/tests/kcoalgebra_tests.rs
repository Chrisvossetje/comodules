use crate::linalg::graded::BasisElement;

#[derive(Debug, Clone)]
struct MockBasisElement;

impl BasisElement for MockBasisElement {}

#[cfg(test)]
mod tests {
    use itertools::Itertools;

    use crate::{
        comodule::kcoalgebra::{kCoalgebra, A0_coalgebra},
        linalg::field::F2,
    };

    #[test]
    fn test_a0() {
        let input = include_str!("../../../examples/kcoalgebras/A(0).txt");

        let (kcoalg, _) = kCoalgebra::<i32, F2>::parse(input).unwrap();

        assert_eq!(kcoalg.coaction, A0_coalgebra().coaction);

        // HASHMAPS are NOT deterministic SO DON'T COMPARE
        // as construct and deconstruct are dependent on insertion order
        assert_eq!(kcoalg.tensor.dimensions, A0_coalgebra().tensor.dimensions);
        assert_eq!(kcoalg.space, A0_coalgebra().space);
    }

    #[test]
    fn test_a2_consistency() {
        let mut comps = vec![];
        for  _ in 0..100 {
            let input = include_str!("../../../examples/kcoalgebras/A(2).txt");
            
            let (kcoalg, _) = kCoalgebra::<i32, F2>::parse(input).unwrap();
            comps.push(kcoalg);
        }
        assert!(comps.iter().all_equal())
    }
}
