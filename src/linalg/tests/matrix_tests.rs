#[cfg(test)]
mod tests {
    use crate::linalg::{
        field::{Field, Fp, F2},
        matrix::Matrix,
        row_matrix::RowMatrix,
    };

    type TestField = Fp<23>;

    #[test]
    pub fn kernel_simple() {
        let matrix = RowMatrix {
            data: vec![
                vec![F2::one(), F2::zero(), F2::one(), F2::one()],
                vec![F2::zero(), F2::one(), F2::zero(), F2::zero()],
                vec![F2::zero(), F2::zero(), F2::one(), F2::one()],
            ],
            domain: 4,
            codomain: 3,
        };

        let kern = matrix.kernel();

        let compare = RowMatrix {
            data: vec![vec![F2::zero(), F2::zero(), F2::one(), F2::one()]],
            domain: 4,
            codomain: 1,
        };
        assert_eq!(kern, compare);
    }

    #[test]
    pub fn kernel_simple_2() {
        let matrix = RowMatrix {
            data: vec![
                // vec![F2::zero(), F2::zero(), F2::zero(), F2::zero()],
                vec![F2::zero(), F2::one(), F2::zero(), F2::one()],
                vec![F2::one(), F2::zero(), F2::one(), F2::one()],
            ],
            domain: 4,
            codomain: 2,
        };

        let kern = matrix.kernel();

        let compare = RowMatrix {
            data: vec![
                vec![F2::one(), F2::zero(), F2::one(), F2::zero()],
                vec![F2::zero(), F2::one(), F2::one(), F2::one()],
            ],
            domain: 4,
            codomain: 2,
        };
        assert_eq!(kern, compare);
    }

    #[test]
    fn kernel_identity() {
        let matrix = RowMatrix::<Fp<23>>::identity(1);
        let kernel = matrix.kernel();

        let comp = RowMatrix::zero(1, 0);
        assert_eq!(kernel, comp)
    }

    #[test]
    fn cokernel_identity() {
        let matrix = RowMatrix::<Fp<23>>::identity(3);
        let kernel = matrix.cokernel();

        let comp = RowMatrix::zero(3, 0);
        assert_eq!(kernel, comp)
    }

    #[test]
    fn test_new() {
        let matrix: RowMatrix<TestField> = RowMatrix::zero(3, 3);
        assert_eq!(matrix.data.len(), 3);
        assert_eq!(matrix.data[0].len(), 3);
        assert_eq!(matrix.domain, 3);
        assert_eq!(matrix.codomain, 3);
        assert!(matrix
            .data
            .iter()
            .all(|row| row.iter().all(|&v| v == TestField::zero())));
    }

    #[test]
    fn test_identity() {
        let identity = RowMatrix::<TestField>::identity(3);
        for i in 0..3 {
            for j in 0..3 {
                if i == j {
                    assert_eq!(identity.data[i][j], TestField::one());
                } else {
                    assert_eq!(identity.data[i][j], TestField::zero());
                }
            }
        }
    }

    #[test]
    fn test_set_and_get() {
        let mut matrix = RowMatrix::zero(3, 3);
        matrix.set(1, 2, TestField { 0: 3 });
        assert_eq!(matrix.get(1, 2), TestField { 0: 3 });
    }

    #[test]
    fn test_transpose() {
        let mut matrix = RowMatrix::zero(2, 3);
        matrix.set(0, 1, TestField { 0: 2 });
        matrix.set(1, 0, TestField { 0: 5 });
        let transpose = matrix.transpose();
        assert_eq!(transpose.get(1, 0), TestField { 0: 2 });
        assert_eq!(transpose.get(0, 1), TestField { 0: 5 });
    }

    #[test]
    fn test_rref_1() {
        let mut matrix = RowMatrix {
            data: vec![
                vec![TestField { 0: 2 }, TestField { 0: 4 }, TestField { 0: 8 }],
                vec![TestField { 0: 1 }, TestField { 0: 2 }, TestField { 0: 4 }],
            ],
            domain: 3,
            codomain: 2,
        };
        matrix.rref();
        let expected = RowMatrix {
            data: vec![
                vec![TestField { 0: 1 }, TestField { 0: 2 }, TestField { 0: 4 }],
                vec![TestField { 0: 0 }, TestField { 0: 0 }, TestField { 0: 0 }],
            ],
            domain: 3,
            codomain: 2,
        };
        assert_eq!(matrix, expected);
    }

    #[test]
    fn test_rref_2() {
        let mut matrix = RowMatrix {
            data: vec![
                vec![TestField { 0: 0 }, TestField { 0: 0 }, TestField { 0: 8 }],
                vec![TestField { 0: 0 }, TestField { 0: 2 }, TestField { 0: 4 }],
            ],
            domain: 3,
            codomain: 2,
        };
        matrix.rref();
        let expected = RowMatrix {
            data: vec![
                vec![TestField { 0: 0 }, TestField { 0: 1 }, TestField { 0: 0 }],
                vec![TestField { 0: 0 }, TestField { 0: 0 }, TestField { 0: 1 }],
            ],
            domain: 3,
            codomain: 2,
        };
        assert_eq!(matrix, expected);
    }

    #[test]
    fn test_rref_3() {
        let mut matrix = RowMatrix {
            data: vec![
                vec![TestField { 0: 1 }, TestField { 0: 1 }, TestField { 0: 0 }],
                vec![TestField { 0: 0 }, TestField { 0: 1 }, TestField { 0: 1 }],
            ],
            domain: 3,
            codomain: 2,
        };
        matrix.rref();
        let expected = RowMatrix {
            data: vec![
                vec![TestField { 0: 1 }, TestField { 0: 0 }, TestField { 0: 22 }],
                vec![TestField { 0: 0 }, TestField { 0: 1 }, TestField { 0: 1 }],
            ],
            domain: 3,
            codomain: 2,
        };
        assert_eq!(matrix, expected);
    }

    #[test]
    fn test_rref_4() {
        let mut matrix = RowMatrix {
            data: vec![vec![TestField { 0: 1 }], vec![TestField { 0: 0 }]],
            domain: 1,
            codomain: 2,
        };
        matrix.rref();
        let expected = RowMatrix {
            data: vec![vec![TestField { 0: 1 }], vec![TestField { 0: 0 }]],
            domain: 1,
            codomain: 2,
        };
        assert_eq!(matrix, expected);
    }

    #[test]
    fn test_rref_7() {
        let mut matrix = RowMatrix {
            data: vec![
                // vec![F2::zero(), F2::zero(), F2::zero(), F2::zero()],
                vec![F2::zero(), F2::one(), F2::zero(), F2::one()],
                vec![F2::one(), F2::zero(), F2::one(), F2::one()],
            ],
            domain: 4,
            codomain: 2,
        };
        matrix.rref();
        let expected = RowMatrix {
            data: vec![
                vec![F2::one(), F2::zero(), F2::one(), F2::one()],
                vec![F2::zero(), F2::one(), F2::zero(), F2::one()],
            ],
            domain: 4,
            codomain: 2,
        };
        assert_eq!(matrix, expected);
    }

    #[test]
    fn test_compose() {
        let matrix1 = RowMatrix {
            data: vec![
                vec![TestField { 0: 1 }, TestField { 0: 2 }],
                vec![TestField { 0: 3 }, TestField { 0: 4 }],
            ],
            domain: 2,
            codomain: 2,
        };
        let mut matrix2 = RowMatrix {
            data: vec![
                vec![TestField { 0: 2 }, TestField { 0: 0 }],
                vec![TestField { 0: 1 }, TestField { 0: 3 }],
            ],
            domain: 2,
            codomain: 2,
        };
        let result = matrix1.compose(&mut matrix2);
        let expected = RowMatrix {
            data: vec![
                vec![TestField { 0: 4 }, TestField { 0: 6 }],
                vec![TestField { 0: 10 }, TestField { 0: 12 }],
            ],
            domain: 2,
            codomain: 2,
        };
        assert_eq!(result, expected);
    }

    #[test]
    fn test_compose_2() {
        let mut matrix1: RowMatrix<F2> = RowMatrix::zero(4, 3);
        let matrix2 = RowMatrix::zero(3, 5);
        let result = matrix2.compose(&mut matrix1);
        assert_eq!(result.codomain, 5);
        assert_eq!(result.domain, 4);
    }

    #[test]
    fn test_pivots_1() {
        let matrix = RowMatrix {
            data: vec![
                vec![TestField { 0: 1 }, TestField { 0: 0 }, TestField { 0: 3 }],
                vec![TestField { 0: 0 }, TestField { 0: 1 }, TestField { 0: 4 }],
            ],
            domain: 3,
            codomain: 2,
        };
        let pivots = matrix.pivots();
        assert_eq!(pivots, vec![(0, 0), (1, 1)]);
    }

    #[test]
    fn test_pivots_2() {
        let matrix = RowMatrix {
            data: vec![
                vec![TestField { 0: 0 }, TestField { 0: 1 }, TestField { 0: 0 }],
                vec![TestField { 0: 0 }, TestField { 0: 0 }, TestField { 0: 1 }],
            ],
            domain: 3,
            codomain: 2,
        };
        let pivots = matrix.pivots();
        assert_eq!(pivots, vec![(1, 0), (2, 1)]);
    }

    #[test]
    fn test_pivots_3() {
        let matrix = RowMatrix {
            data: vec![
                vec![F2::one(), F2::zero(), F2::one(), F2::one()],
                vec![F2::zero(), F2::one(), F2::zero(), F2::one()],
            ],
            domain: 4,
            codomain: 2,
        };
        let pivots = matrix.pivots();
        assert_eq!(pivots, vec![(0, 0), (1, 1)]);
    }

    #[test]
    fn test_vstack_1() {
        let mut matrix1 = RowMatrix {
            data: vec![vec![TestField { 0: 1 }], vec![TestField { 0: 2 }]],
            domain: 1,
            codomain: 2,
        };
        let mut matrix2 = RowMatrix {
            data: vec![vec![TestField { 0: 3 }], vec![TestField { 0: 4 }]],
            domain: 1,
            codomain: 2,
        };
        matrix1.vstack(&mut matrix2);
        let expected = RowMatrix {
            data: vec![
                vec![TestField { 0: 1 }],
                vec![TestField { 0: 2 }],
                vec![TestField { 0: 3 }],
                vec![TestField { 0: 4 }],
            ],
            domain: 1,
            codomain: 4,
        };
        assert_eq!(matrix1, expected);
    }

    #[test]
    fn test_vstack_2() {
        let mut matrix1 = RowMatrix {
            data: vec![vec![F2::one(), F2::one()], vec![F2::zero(), F2::zero()]],
            domain: 2,
            codomain: 2,
        };
        let mut matrix2 = RowMatrix {
            data: vec![vec![F2::zero(), F2::one()]],
            domain: 2,
            codomain: 1,
        };
        matrix1.vstack(&mut matrix2);
        let expected = RowMatrix {
            data: vec![
                vec![F2::one(), F2::one()],
                vec![F2::zero(), F2::zero()],
                vec![F2::zero(), F2::one()],
            ],
            domain: 2,
            codomain: 3,
        };
        assert_eq!(matrix1, expected);
    }

    #[test]
    fn test_block_sum_1() {
        let mut matrix1 = RowMatrix {
            data: vec![vec![TestField { 0: 1 }], vec![TestField { 0: 2 }]],
            domain: 1,
            codomain: 2,
        };
        let mut matrix2 = RowMatrix {
            data: vec![vec![TestField { 0: 3 }], vec![TestField { 0: 4 }]],
            domain: 1,
            codomain: 2,
        };
        matrix1.block_sum(&mut matrix2);
        let expected = RowMatrix {
            data: vec![
                vec![TestField { 0: 1 }, TestField { 0: 0 }],
                vec![TestField { 0: 2 }, TestField { 0: 0 }],
                vec![TestField { 0: 0 }, TestField { 0: 3 }],
                vec![TestField { 0: 0 }, TestField { 0: 4 }],
            ],
            domain: 2,
            codomain: 4,
        };
        assert_eq!(matrix1, expected);
    }

    #[test]
    fn test_block_sum_2() {
        let mut matrix1 = RowMatrix {
            data: vec![vec![TestField { 0: 1 }, TestField { 0: 2 }]],
            domain: 2,
            codomain: 1,
        };
        let mut matrix2 = RowMatrix {
            data: vec![vec![TestField { 0: 3 }], vec![TestField { 0: 4 }]],
            domain: 1,
            codomain: 2,
        };
        matrix1.block_sum(&mut matrix2);
        let expected = RowMatrix {
            data: vec![
                vec![TestField { 0: 1 }, TestField { 0: 2 }, TestField { 0: 0 }],
                vec![TestField { 0: 0 }, TestField { 0: 0 }, TestField { 0: 3 }],
                vec![TestField { 0: 0 }, TestField { 0: 0 }, TestField { 0: 4 }],
            ],
            domain: 3,
            codomain: 3,
        };
        assert_eq!(matrix1, expected);
    }

    #[test]
    fn test_get() {
        let data = vec![
            vec![TestField { 0: 1 }, TestField { 0: 2 }],
            vec![TestField { 0: 3 }, TestField { 0: 4 }],
        ];
        let matrix = RowMatrix {
            data,
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
        let data = vec![
            vec![TestField { 0: 1 }, TestField { 0: 2 }],
            vec![TestField { 0: 3 }, TestField { 0: 4 }],
        ];
        let mut matrix = RowMatrix {
            data,
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
        let data = vec![
            vec![TestField { 0: 1 }, TestField { 0: 2 }],
            vec![TestField { 0: 3 }, TestField { 0: 4 }],
        ];
        let mut matrix = RowMatrix {
            data,
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
        let data = vec![
            vec![TestField { 0: 1 }, TestField { 0: 2 }],
            vec![TestField { 0: 3 }, TestField { 0: 4 }],
        ];
        let matrix = RowMatrix {
            data,
            domain: 2,
            codomain: 2,
        };

        assert_eq!(matrix.get_row(0), &[TestField { 0: 1 }, TestField { 0: 2 }]);
        assert_eq!(matrix.get_row(1), &[TestField { 0: 3 }, TestField { 0: 4 }]);
    }

    #[test]
    fn test_set_row() {
        let data = vec![
            vec![TestField { 0: 1 }, TestField { 0: 2 }],
            vec![TestField { 0: 3 }, TestField { 0: 4 }],
        ];
        let mut matrix = RowMatrix {
            data,
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
        let data = vec![
            vec![TestField { 0: 1 }, TestField { 0: 2 }],
            vec![TestField { 0: 3 }, TestField { 0: 4 }],
        ];
        let matrix = RowMatrix {
            data,
            domain: 2,
            codomain: 2,
        };

        assert_eq!(matrix.domain(), 2);
        assert_eq!(matrix.codomain(), 2);
    }

    #[test]
    fn test_first_non_zero_entry() {
        let data = vec![
            vec![TestField { 0: 0 }, TestField { 0: 0 }],
            vec![TestField { 0: 0 }, TestField { 0: 4 }],
        ];
        let matrix = RowMatrix {
            data,
            domain: 2,
            codomain: 2,
        };

        assert_eq!(matrix.first_non_zero_entry(), Some((1, 1)));

        let data = vec![
            vec![TestField { 0: 0 }, TestField { 0: 0 }],
            vec![TestField { 0: 0 }, TestField { 0: 0 }],
        ];
        let matrix = RowMatrix {
            data,
            domain: 2,
            codomain: 2,
        };

        assert_eq!(matrix.first_non_zero_entry(), None);
    }

    #[test]
    fn test_first_non_zero_entry_edge_case() {
        let data = vec![
            vec![TestField { 0: 5 }, TestField { 0: 0 }],
            vec![TestField { 0: 0 }, TestField { 0: 4 }],
        ];
        let matrix = RowMatrix {
            data,
            domain: 2,
            codomain: 2,
        };

        assert_eq!(matrix.first_non_zero_entry(), Some((0, 0)));
    }
}
