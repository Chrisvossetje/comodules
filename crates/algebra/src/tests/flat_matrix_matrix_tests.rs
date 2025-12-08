#[cfg(test)]
mod tests {
    use crate::{matrices::flat_matrix::FlatMatrix, matrix::Matrix, ring::CRing, rings::finite_fields::{F2, Fp}};

    type TestField = Fp<23>;

    #[test]
    fn test_vstack_1() {
        let mut matrix1 = FlatMatrix {
            data: vec![TestField { 0: 1 }, TestField { 0: 2 }],
            domain: 2,
            codomain: 1,
        };
        let mut matrix2 = FlatMatrix {
            data: vec![TestField { 0: 3 }, TestField { 0: 4 }],
            domain: 2,
            codomain: 1,
        };
        matrix1.vstack(&mut matrix2);
        let expected = FlatMatrix {
            data: vec![
                TestField { 0: 1 },
                TestField { 0: 2 },
                TestField { 0: 3 },
                TestField { 0: 4 },
            ],
            domain: 2,
            codomain: 2,
        };
        assert_eq!(matrix1, expected);
    }

    #[test]
    fn test_vstack_2() {
        let mut matrix1 = FlatMatrix {
            data: vec![F2::one(), F2::one(), F2::zero(), F2::zero()],
            domain: 2,
            codomain: 2,
        };
        let mut matrix2 = FlatMatrix {
            data: vec![F2::zero(), F2::one()],
            domain: 2,
            codomain: 1,
        };
        matrix1.vstack(&mut matrix2);
        let expected = FlatMatrix {
            data: vec![
                F2::one(),
                F2::one(),
                F2::zero(),
                F2::zero(),
                F2::zero(),
                F2::one(),
            ],
            domain: 2,
            codomain: 3,
        };
        assert_eq!(matrix1, expected);
    }

    #[test]
    fn test_block_sum_1() {
        let mut matrix1 = FlatMatrix {
            data: vec![TestField { 0: 1 }, TestField { 0: 2 }],
            domain: 1,
            codomain: 2,
        };
        let mut matrix2 = FlatMatrix {
            data: vec![TestField { 0: 3 }, TestField { 0: 4 }],
            domain: 1,
            codomain: 2,
        };
        matrix1.block_sum(&mut matrix2);
        let expected = FlatMatrix {
            data: vec![
                TestField { 0: 1 },
                TestField { 0: 0 },
                TestField { 0: 2 },
                TestField { 0: 0 },
                TestField { 0: 0 },
                TestField { 0: 3 },
                TestField { 0: 0 },
                TestField { 0: 4 },
            ],
            domain: 2,
            codomain: 4,
        };
        assert_eq!(matrix1, expected);
    }

    #[test]
    fn test_block_sum_2() {
        let mut matrix1 = FlatMatrix {
            data: vec![TestField { 0: 1 }, TestField { 0: 2 }],
            domain: 2,
            codomain: 1,
        };
        let mut matrix2 = FlatMatrix {
            data: vec![TestField { 0: 3 }, TestField { 0: 4 }],
            domain: 1,
            codomain: 2,
        };
        matrix1.block_sum(&mut matrix2);
        let expected = FlatMatrix {
            data: vec![
                TestField { 0: 1 },
                TestField { 0: 2 },
                TestField { 0: 0 },
                TestField { 0: 0 },
                TestField { 0: 0 },
                TestField { 0: 3 },
                TestField { 0: 0 },
                TestField { 0: 0 },
                TestField { 0: 4 },
            ],
            domain: 3,
            codomain: 3,
        };
        assert_eq!(matrix1, expected);
    }

    #[test]
    fn test_get() {
        let matrix = FlatMatrix {
            data: vec![
                TestField { 0: 1 },
                TestField { 0: 2 },
                TestField { 0: 3 },
                TestField { 0: 4 },
            ],
            domain: 2,
            codomain: 2,
        };

        assert_eq!(matrix.get(0, 0), TestField { 0: 1 });
        assert_eq!(matrix.get(1, 0), TestField { 0: 2 });
        assert_eq!(matrix.get(0, 1), TestField { 0: 3 });
        assert_eq!(matrix.get(1, 1), TestField { 0: 4 });
    }

    #[test]
    fn test_set() {
        let mut matrix = FlatMatrix {
            data: vec![
                TestField { 0: 1 },
                TestField { 0: 2 },
                TestField { 0: 3 },
                TestField { 0: 4 },
            ],
            domain: 2,
            codomain: 2,
        };

        matrix.set(0, 0, TestField { 0: 10 });
        matrix.set(1, 1, TestField { 0: 20 });

        assert_eq!(matrix.get(0, 0), TestField { 0: 10 });
        assert_eq!(matrix.get(1, 1), TestField { 0: 20 });
    }

    #[test]
    fn test_add_at() {
        let mut matrix = FlatMatrix {
            data: vec![
                TestField { 0: 1 },
                TestField { 0: 2 },
                TestField { 0: 3 },
                TestField { 0: 4 },
            ],
            domain: 2,
            codomain: 2,
        };

        matrix.add_at(0, 0, TestField { 0: 5 });
        matrix.add_at(1, 1, TestField { 0: 3 });

        assert_eq!(matrix.get(0, 0), TestField { 0: 6 }); // 1 + 5
        assert_eq!(matrix.get(1, 1), TestField { 0: 7 }); // 4 + 3
    }

    #[test]
    fn test_get_row() {
        let matrix = FlatMatrix {
            data: vec![
                TestField { 0: 1 },
                TestField { 0: 2 },
                TestField { 0: 3 },
                TestField { 0: 4 },
            ],
            domain: 2,
            codomain: 2,
        };

        assert_eq!(matrix.get_row(0), &[TestField { 0: 1 }, TestField { 0: 2 }]);
        assert_eq!(matrix.get_row(1), &[TestField { 0: 3 }, TestField { 0: 4 }]);
    }

    #[test]
    fn test_set_row() {
        let mut matrix = FlatMatrix {
            data: vec![
                TestField { 0: 1 },
                TestField { 0: 2 },
                TestField { 0: 3 },
                TestField { 0: 4 },
            ],
            domain: 2,
            codomain: 2,
        };

        let new_row = vec![TestField { 0: 10 }, TestField { 0: 20 }];
        matrix.set_row(0, &new_row);

        assert_eq!(
            matrix.get_row(0),
            &[TestField { 0: 10 }, TestField { 0: 20 }]
        );
    }

    #[test]
    fn test_domain_codomain() {
        let matrix = FlatMatrix {
            data: vec![
                TestField { 0: 1 },
                TestField { 0: 2 },
                TestField { 0: 3 },
                TestField { 0: 4 },
            ],
            domain: 2,
            codomain: 2,
        };

        assert_eq!(matrix.domain(), 2);
        assert_eq!(matrix.codomain(), 2);
    }
}