#[cfg(test)]
mod tests {
    use crate::linalg::{
        field::{Field, Fp, F2},
        flat_matrix::FlatMatrix, // Update the import to use FlatMatrix
        matrix::Matrix,
    };

    type TestField = Fp<23>;

    #[test]
    pub fn kernel_simple() {
        let matrix = FlatMatrix {
            data: vec![
                F2::one(),
                F2::zero(),
                F2::one(),
                F2::one(),
                F2::zero(),
                F2::one(),
                F2::zero(),
                F2::zero(),
                F2::zero(),
                F2::zero(),
                F2::one(),
                F2::one(),
            ],
            domain: 4,
            codomain: 3,
        };

        let kern = matrix.kernel();

        let compare = FlatMatrix {
            data: vec![F2::zero(), F2::zero(), F2::one(), F2::one()],
            domain: 4,
            codomain: 1,
        };
        assert_eq!(kern, compare);
    }

    #[test]
    pub fn kernel_simple_2() {
        let matrix = FlatMatrix {
            data: vec![
                F2::zero(),
                F2::one(),
                F2::zero(),
                F2::one(),
                F2::one(),
                F2::zero(),
                F2::one(),
                F2::one(),
            ],
            domain: 4,
            codomain: 2,
        };

        let kern = matrix.kernel();

        let compare = FlatMatrix {
            data: vec![
                F2::one(),
                F2::zero(),
                F2::one(),
                F2::zero(),
                F2::zero(),
                F2::one(),
                F2::one(),
                F2::one(),
            ],
            domain: 4,
            codomain: 2,
        };
        assert_eq!(kern, compare);
    }

    #[test]
    fn kernel_identity() {
        let matrix = FlatMatrix::<Fp<23>>::identity(1);
        let kernel = matrix.kernel();

        let comp = FlatMatrix::zero(1, 0);
        assert_eq!(kernel, comp)
    }

    #[test]
    fn cokernel_identity() {
        let matrix = FlatMatrix::<Fp<23>>::identity(3);
        let kernel = matrix.cokernel();

        let comp = FlatMatrix::zero(3, 0);
        assert_eq!(kernel, comp)
    }

    #[test]
    fn test_new() {
        let matrix: FlatMatrix<TestField> = FlatMatrix::zero(3, 3);
        assert_eq!(matrix.data.len(), 9); // Data is flat
        assert_eq!(matrix.domain, 3);
        assert_eq!(matrix.codomain, 3);
        assert!(matrix.data.iter().all(|&v| v == TestField::zero()));
    }

    #[test]
    fn test_identity() {
        let identity = FlatMatrix::<TestField>::identity(3);
        for i in 0..3 {
            for j in 0..3 {
                if i == j {
                    assert_eq!(identity.get(i, j), TestField::one());
                } else {
                    assert_eq!(identity.get(i, j), TestField::zero());
                }
            }
        }
    }

    #[test]
    fn test_set_and_get() {
        let mut matrix = FlatMatrix::zero(3, 3);
        matrix.set(1, 2, TestField { 0: 3 });
        assert_eq!(matrix.get(1, 2), TestField { 0: 3 });
    }

    #[test]
    fn test_transpose() {
        let mut matrix = FlatMatrix::zero(2, 3);
        matrix.set(0, 1, TestField { 0: 2 });
        matrix.set(1, 0, TestField { 0: 5 });
        let transpose = matrix.transpose();
        assert_eq!(transpose.get(1, 0), TestField { 0: 2 });
        assert_eq!(transpose.get(0, 1), TestField { 0: 5 });
    }

    #[test]
    fn test_rref_1() {
        let mut matrix = FlatMatrix {
            data: vec![
                TestField { 0: 2 },
                TestField { 0: 4 },
                TestField { 0: 8 },
                TestField { 0: 1 },
                TestField { 0: 2 },
                TestField { 0: 4 },
            ],
            domain: 3,
            codomain: 2,
        };
        matrix.rref();
        let expected = FlatMatrix {
            data: vec![
                TestField { 0: 1 },
                TestField { 0: 2 },
                TestField { 0: 4 },
                TestField { 0: 0 },
                TestField { 0: 0 },
                TestField { 0: 0 },
            ],
            domain: 3,
            codomain: 2,
        };
        assert_eq!(matrix, expected);
    }

    #[test]
    fn test_rref_2() {
        let mut matrix = FlatMatrix {
            data: vec![
                TestField { 0: 0 },
                TestField { 0: 0 },
                TestField { 0: 8 },
                TestField { 0: 0 },
                TestField { 0: 2 },
                TestField { 0: 4 },
            ],
            domain: 3,
            codomain: 2,
        };
        matrix.rref();
        let expected = FlatMatrix {
            data: vec![
                TestField { 0: 0 },
                TestField { 0: 1 },
                TestField { 0: 0 },
                TestField { 0: 0 },
                TestField { 0: 0 },
                TestField { 0: 1 },
            ],
            domain: 3,
            codomain: 2,
        };
        assert_eq!(matrix, expected);
    }

    #[test]
    fn test_rref_3() {
        let mut matrix = FlatMatrix {
            data: vec![
                TestField { 0: 1 },
                TestField { 0: 1 },
                TestField { 0: 0 },
                TestField { 0: 0 },
                TestField { 0: 1 },
                TestField { 0: 1 },
            ],
            domain: 3,
            codomain: 2,
        };
        matrix.rref();
        let expected = FlatMatrix {
            data: vec![
                TestField { 0: 1 },
                TestField { 0: 0 },
                TestField { 0: 22 },
                TestField { 0: 0 },
                TestField { 0: 1 },
                TestField { 0: 1 },
            ],
            domain: 3,
            codomain: 2,
        };
        assert_eq!(matrix, expected);
    }

    #[test]
    fn test_rref_4() {
        let mut matrix = FlatMatrix {
            data: vec![TestField { 0: 1 }, TestField { 0: 0 }],
            domain: 1,
            codomain: 2,
        };
        matrix.rref();
        let expected = FlatMatrix {
            data: vec![TestField { 0: 1 }, TestField { 0: 0 }],
            domain: 1,
            codomain: 2,
        };
        assert_eq!(matrix, expected);
    }

    #[test]
    fn test_rref_5() {
        let mut matrix = FlatMatrix {
            data: vec![
                // F2::zero(), F2::zero(), F2::zero(), F2::zero(),
                F2::zero(),
                F2::one(),
                F2::zero(),
                F2::one(),
                F2::one(),
                F2::zero(),
                F2::one(),
                F2::one(),
            ],
            domain: 4,
            codomain: 2,
        };
        matrix.rref();
        let expected = FlatMatrix {
            data: vec![
                F2::one(),
                F2::zero(),
                F2::one(),
                F2::one(),
                F2::zero(),
                F2::one(),
                F2::zero(),
                F2::one(),
            ],
            domain: 4,
            codomain: 2,
        };
        assert_eq!(matrix, expected);
    }

    #[test]
    fn test_compose_1() {
        let matrix1 = FlatMatrix {
            data: vec![
                TestField { 0: 1 },
                TestField { 0: 2 },
                TestField { 0: 3 },
                TestField { 0: 4 },
            ],
            domain: 2,
            codomain: 2,
        };
        let mut matrix2 = FlatMatrix {
            data: vec![
                TestField { 0: 2 },
                TestField { 0: 0 },
                TestField { 0: 1 },
                TestField { 0: 3 },
            ],
            domain: 2,
            codomain: 2,
        };
        let result = matrix1.compose(&mut matrix2);
        let expected = FlatMatrix {
            data: vec![
                TestField { 0: 4 },
                TestField { 0: 6 },
                TestField { 0: 10 },
                TestField { 0: 12 },
            ],
            domain: 2,
            codomain: 2,
        };
        assert_eq!(result, expected);
    }

    #[test]
    fn test_compose_2() {
        let mut matrix1: FlatMatrix<F2> = FlatMatrix::zero(4, 3);
        let matrix2 = FlatMatrix::zero(3, 5);
        let result = matrix2.compose(&mut matrix1);
        assert_eq!(result.codomain, 5);
        assert_eq!(result.domain, 4);
    }

    #[test]
    fn test_pivots_1() {
        let matrix = FlatMatrix {
            data: vec![
                TestField { 0: 1 },
                TestField { 0: 0 },
                TestField { 0: 3 },
                TestField { 0: 0 },
                TestField { 0: 1 },
                TestField { 0: 4 },
            ],
            domain: 3,
            codomain: 2,
        };
        let pivots = matrix.pivots();
        assert_eq!(pivots, vec![(0, 0), (1, 1)]);
    }

    #[test]
    fn test_pivots_2() {
        let matrix = FlatMatrix {
            data: vec![
                TestField { 0: 0 },
                TestField { 0: 1 },
                TestField { 0: 0 },
                TestField { 0: 0 },
                TestField { 0: 0 },
                TestField { 0: 1 },
            ],
            domain: 3,
            codomain: 2,
        };
        let pivots = matrix.pivots();
        assert_eq!(pivots, vec![(1, 0), (2, 1)]);
    }

    #[test]
    fn test_pivots_3() {
        let matrix = FlatMatrix {
            data: vec![
                F2::one(),
                F2::zero(),
                F2::one(),
                F2::one(),
                F2::zero(),
                F2::one(),
                F2::zero(),
                F2::one(),
            ],
            domain: 4,
            codomain: 2,
        };
        let pivots = matrix.pivots();
        assert_eq!(pivots, vec![(0, 0), (1, 1)]);
    }

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
        let stack = matrix1.vstack(&mut matrix2);
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
        assert_eq!(stack, expected);
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
        let stack = matrix1.vstack(&mut matrix2);
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
        assert_eq!(stack, expected);
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
        let block = matrix1.block_sum(&mut matrix2);
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
        assert_eq!(block, expected);
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
        let block = matrix1.block_sum(&mut matrix2);
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
        assert_eq!(block, expected);
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

    #[test]
    fn test_first_non_zero_entry() {
        let matrix = FlatMatrix {
            data: vec![
                TestField { 0: 0 },
                TestField { 0: 0 },
                TestField { 0: 0 },
                TestField { 0: 4 },
            ],
            domain: 2,
            codomain: 2,
        };

        assert_eq!(matrix.first_non_zero_entry(), Some((1, 1)));

        let matrix = FlatMatrix {
            data: vec![
                TestField { 0: 0 },
                TestField { 0: 0 },
                TestField { 0: 0 },
                TestField { 0: 0 },
            ],
            domain: 2,
            codomain: 2,
        };

        assert_eq!(matrix.first_non_zero_entry(), None);
    }

    #[test]
    fn test_first_non_zero_entry_edge_case() {
        let matrix = FlatMatrix {
            data: vec![
                TestField { 0: 5 },
                TestField { 0: 0 },
                TestField { 0: 0 },
                TestField { 0: 4 },
            ],
            domain: 2,
            codomain: 2,
        };

        assert_eq!(matrix.first_non_zero_entry(), Some((0, 0)));
    }
}
