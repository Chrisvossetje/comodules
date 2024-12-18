#[cfg(test)]
mod tests {
    use crate::linalg::{
        field::{Field, Fp, F2},
        matrix::{FieldMatrix, Matrix},
    };

    type TestField = Fp<23>;

    #[test]
    pub fn kernel_simple() {
        let matrix = FieldMatrix {
            data: vec![
                vec![F2::one(), F2::zero(), F2::one(), F2::one()],
                vec![F2::zero(), F2::one(), F2::zero(), F2::zero()],
                vec![F2::zero(), F2::zero(), F2::one(), F2::one()],
            ],
            domain: 4,
            codomain: 3,
        };

        let kern = matrix.kernel();

        let compare = FieldMatrix {
            data: vec![vec![F2::zero(), F2::zero(), F2::one(), F2::one()]],
            domain: 4,
            codomain: 1,
        };
        assert_eq!(kern, compare);
    }

    #[test]
    fn kernel_identity() {
        let matrix = FieldMatrix::<Fp<23>>::identity(1);
        let kernel = matrix.kernel();

        let comp = FieldMatrix::zero(1, 0);
        assert_eq!(kernel, comp)
    }

    #[test]
    fn cokernel_identity() {
        let matrix = FieldMatrix::<Fp<23>>::identity(3);
        let kernel = matrix.cokernel();

        let comp = FieldMatrix::zero(3, 0);
        assert_eq!(kernel, comp)
    }

    #[test]
    fn test_new() {
        let matrix: FieldMatrix<TestField> = FieldMatrix::zero(3, 3);
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
        let identity = FieldMatrix::<TestField>::identity(3);
        for i in 0..3 {
            for j in 0..3 {
                if i == j {
                    assert_eq!(identity.data[i][j], TestField::one());
                } else {
                    assert_eq!(identity.data[i][j], TestField::zero());
                }
            }
        }
        let a = TestField { 0: 2 };
    }

    #[test]
    fn test_set_and_get() {
        let mut matrix = FieldMatrix::zero(3, 3);
        matrix.set(1, 2, TestField { 0: 3 });
        assert_eq!(matrix.get(1, 2), TestField { 0: 3 });
    }

    #[test]
    fn test_transpose() {
        let mut matrix = FieldMatrix::zero(2, 3);
        matrix.set(0, 1, TestField { 0: 2 });
        matrix.set(1, 0, TestField { 0: 5 });
        let transpose = matrix.transpose();
        assert_eq!(transpose.get(1, 0), TestField { 0: 2 });
        assert_eq!(transpose.get(0, 1), TestField { 0: 5 });
    }

    #[test]
    fn test_rref_1() {
        let mut matrix = FieldMatrix {
            data: vec![
                vec![TestField { 0: 2 }, TestField { 0: 4 }, TestField { 0: 8 }],
                vec![TestField { 0: 1 }, TestField { 0: 2 }, TestField { 0: 4 }],
            ],
            domain: 3,
            codomain: 2,
        };
        matrix.rref();
        let expected = FieldMatrix {
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
        let mut matrix = FieldMatrix {
            data: vec![
                vec![TestField { 0: 0 }, TestField { 0: 0 }, TestField { 0: 8 }],
                vec![TestField { 0: 0 }, TestField { 0: 2 }, TestField { 0: 4 }],
            ],
            domain: 3,
            codomain: 2,
        };
        matrix.rref();
        let expected = FieldMatrix {
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
        let mut matrix = FieldMatrix {
            data: vec![
                vec![TestField { 0: 1 }, TestField { 0: 1 }, TestField { 0: 0 }],
                vec![TestField { 0: 0 }, TestField { 0: 1 }, TestField { 0: 1 }],
            ],
            domain: 3,
            codomain: 2,
        };
        matrix.rref();
        let expected = FieldMatrix {
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
        let mut matrix = FieldMatrix {
            data: vec![vec![TestField { 0: 1 }], vec![TestField { 0: 0 }]],
            domain: 1,
            codomain: 2,
        };
        matrix.rref();
        let expected = FieldMatrix {
            data: vec![vec![TestField { 0: 1 }], vec![TestField { 0: 0 }]],
            domain: 1,
            codomain: 2,
        };
        assert_eq!(matrix, expected);
    }

    #[test]
    fn test_compose() {
        let matrix1 = FieldMatrix {
            data: vec![
                vec![TestField { 0: 1 }, TestField { 0: 2 }],
                vec![TestField { 0: 3 }, TestField { 0: 4 }],
            ],
            domain: 2,
            codomain: 2,
        };
        let mut matrix2 = FieldMatrix {
            data: vec![
                vec![TestField { 0: 2 }, TestField { 0: 0 }],
                vec![TestField { 0: 1 }, TestField { 0: 3 }],
            ],
            domain: 2,
            codomain: 2,
        };
        let result = matrix1.compose(&mut matrix2);
        let expected = FieldMatrix {
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
        let mut matrix1: FieldMatrix<F2> = FieldMatrix::zero(4, 3);
        let matrix2 = FieldMatrix::zero(3, 5);
        let result = matrix2.compose(&mut matrix1);
        assert_eq!(result.codomain, 5);
        assert_eq!(result.domain, 4);
    }

    #[test]
    fn test_pivots_1() {
        let matrix = FieldMatrix {
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
        let matrix = FieldMatrix {
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
    fn test_vstack_1() {
        let mut matrix1 = FieldMatrix {
            data: vec![vec![TestField { 0: 1 }], vec![TestField { 0: 2 }]],
            domain: 1,
            codomain: 2,
        };
        let mut matrix2 = FieldMatrix {
            data: vec![vec![TestField { 0: 3 }], vec![TestField { 0: 4 }]],
            domain: 1,
            codomain: 2,
        };
        matrix1.vstack(&mut matrix2);
        let expected = FieldMatrix {
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
        let mut matrix1 = FieldMatrix {
            data: vec![vec![F2::one(), F2::one()], vec![F2::zero(), F2::zero()]],
            domain: 2,
            codomain: 2,
        };
        let mut matrix2 = FieldMatrix {
            data: vec![vec![F2::zero(), F2::one()]],
            domain: 2,
            codomain: 1,
        };
        matrix1.vstack(&mut matrix2);
        let expected = FieldMatrix {
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
        let mut matrix1 = FieldMatrix {
            data: vec![vec![TestField { 0: 1 }], vec![TestField { 0: 2 }]],
            domain: 1,
            codomain: 2,
        };
        let mut matrix2 = FieldMatrix {
            data: vec![vec![TestField { 0: 3 }], vec![TestField { 0: 4 }]],
            domain: 1,
            codomain: 2,
        };
        matrix1.block_sum(&mut matrix2);
        let expected = FieldMatrix {
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
        let mut matrix1 = FieldMatrix {
            data: vec![vec![TestField { 0: 1 }, TestField { 0: 2 }]],
            domain: 2,
            codomain: 1,
        };
        let mut matrix2 = FieldMatrix {
            data: vec![vec![TestField { 0: 3 }], vec![TestField { 0: 4 }]],
            domain: 1,
            codomain: 2,
        };
        matrix1.block_sum(&mut matrix2);
        let expected = FieldMatrix {
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
}
