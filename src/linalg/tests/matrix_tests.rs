#[cfg(test)]
mod tests {
    use crate::linalg::{field::{Field, Fp, F2}, matrix::{FieldMatrix, Matrix}};
    use super::*;
    type TestField = Fp::<23>;

    #[test]
    pub fn kernel_simple() {
        let mut matrix = vec![
            vec![F2::one(), F2::zero(), F2::one(), F2::one()],
            vec![F2::zero(), F2::one(), F2::zero(), F2::zero()],
            vec![F2::zero(), F2::zero(), F2::one(), F2::one()],
        ];
        
        println!("Original Matrix:");
        for row in &matrix {
            println!("{:?}", row);
        }
        matrix.kernel();
        println!("Rref Matrix:");
        for row in matrix.clone() {
            println!("{:?}", row);
        }
        let kern = matrix.rref_kernel();

        println!("Kernel Matrix:");
        for row in &kern {
            println!("{:?}", row);
        }
    }
    
    #[test]
    pub fn rref_simple() {
        
        let mut matrix = vec![
        vec![F2(1), F2(1), F2(0), F2(1)],
        vec![F2(1), F2(1), F2(0), F2(0)],
        vec![F2(1), F2(0), F2(0), F2(1)],
        ];
        
        println!("Original Matrix:");
        for row in &matrix {
            println!("{:?}", row);
        }

        matrix.rref();

        println!("\nReduced Row Echelon Form:");
        for row in &matrix {
            println!("{:?}", row);
        }
    }



    #[test]
    fn test_new() {
        let matrix: FieldMatrix<TestField> = FieldMatrix::zero(3, 3);
        assert_eq!(matrix.len(), 3);
        assert_eq!(matrix[0].len(), 3);
        assert!(matrix.iter().all(|row| row.iter().all(|&v| v == TestField::zero())));
    }

    #[test]
    fn test_identity() {
        let identity = Vec::<Vec<TestField>>::identity(3);
        for i in 0..3 {
            for j in 0..3 {
                if i == j {
                    assert_eq!(identity[i][j], TestField::one());
                } else {
                    assert_eq!(identity[i][j], TestField::zero());
                }
            }
        }
        let a = TestField{0: 2};
    }

    #[test]
    fn test_set_and_get() {
        let mut matrix = Vec::<Vec<TestField>>::zero(3, 3);
        matrix.set(1, 2, TestField{0:3});
        assert_eq!(matrix.get(1, 2), TestField{0:3});
    }

    #[test]
    fn test_transpose() {
        let mut matrix = Vec::<Vec<TestField>>::zero(2, 3);
        matrix.set(0, 1, TestField{0:2});
        matrix.set(1, 0, TestField{0:5});
        matrix.transpose();
        assert_eq!(matrix.get(1, 0), TestField{0:2});
        assert_eq!(matrix.get(0, 1), TestField{0:5});
    }

    #[test]
    fn test_rref_1() {
        let mut matrix = vec![
            vec![TestField{0:2}, TestField{0:4}, TestField{0:8}],
            vec![TestField{0:1}, TestField{0:2}, TestField{0:4}],
        ];
        matrix.rref();
        let expected = vec![
            vec![TestField{0:1}, TestField{0:2}, TestField{0:4}],
            vec![TestField{0:0}, TestField{0:0}, TestField{0:0}],
        ];
        assert_eq!(matrix, expected);
    }

    #[test]
    fn test_rref_2() {
        let mut matrix = vec![
            vec![TestField{0:0}, TestField{0:0}, TestField{0:8}],
            vec![TestField{0:0}, TestField{0:2}, TestField{0:4}],
        ];
        matrix.rref();
        let expected = vec![
            vec![TestField{0:0}, TestField{0:1}, TestField{0:0}],
            vec![TestField{0:0}, TestField{0:0}, TestField{0:1}],
        ];
        assert_eq!(matrix, expected);
    }


    #[test]
    fn test_rref_3() {
        let mut matrix = vec![
            vec![TestField{0:1}, TestField{0:1}, TestField{0:0}],
            vec![TestField{0:0}, TestField{0:1}, TestField{0:1}],
        ];
        matrix.rref();
        let expected = vec![
            vec![TestField{0:1}, TestField{0:0}, TestField{0:22}],
            vec![TestField{0:0}, TestField{0:1}, TestField{0:1}],
        ];
        assert_eq!(matrix, expected);
    }

    #[test]
    fn test_rref_kernel() {
        let matrix = vec![
            vec![TestField{0:1}, TestField{0:0}, TestField{0:22}],
            vec![TestField{0:0}, TestField{0:1}, TestField{0:1}],
        ];
        let kernel = matrix.rref_kernel();
        let expected = vec![vec![TestField{0: 22}, TestField{0:0}, TestField{0:1}]];
        assert_eq!(kernel, expected);
    }

    #[test]
    fn test_compose() {
        let mut matrix1 = vec![
            vec![TestField{0:1}, TestField{0:2}],
            vec![TestField{0:3}, TestField{0:4}],
        ];
        let mut matrix2 = vec![
            vec![TestField{0:2}, TestField{0:0}],
            vec![TestField{0:1}, TestField{0:3}],
        ];
        let result = matrix1.compose(&mut matrix2);
        let expected = vec![
            vec![TestField{0:4}, TestField{0:6}],
            vec![TestField{0:10}, TestField{0:12}],
        ];
        assert_eq!(result, expected);
    }

    #[test]
    fn test_pivots() {
        let matrix = vec![
            vec![TestField{0:1}, TestField{0:0}, TestField{0:3}],
            vec![TestField{0:0}, TestField{0:1}, TestField{0:4}],
        ];
        let pivots = matrix.pivots();
        assert_eq!(pivots, vec![(0, 0), (1, 1)]);
    }

    #[test]
    fn test_vstack() {
        let mut matrix1 = vec![vec![TestField{0:1}], vec![TestField{0:2}]];
        let mut matrix2 = vec![vec![TestField{0:3}], vec![TestField{0:4}]];
        matrix1.vstack(&mut matrix2);
        let expected = vec![
            vec![TestField{0:1}],
            vec![TestField{0:2}],
            vec![TestField{0:3}],
            vec![TestField{0:4}],
        ];
        assert_eq!(matrix1, expected);
    }

    #[test]
    fn test_block_sum_1() {
        let mut matrix1 = vec![vec![TestField{0:1}], vec![TestField{0:2}]];
        let mut matrix2 = vec![vec![TestField{0:3}], vec![TestField{0:4}]];
        matrix1.block_sum(&mut matrix2);
        let expected = vec![
            vec![TestField{0:1}, TestField{0:0}],
            vec![TestField{0:2}, TestField{0:0}],
            vec![TestField{0:0}, TestField{0:3}],
            vec![TestField{0:0}, TestField{0:4}],
        ];
        assert_eq!(matrix1, expected);
    }

    #[test]
    fn test_block_sum_2() {
        let mut matrix1 = vec![vec![TestField{0:1}, TestField{0:2}]];
        let mut matrix2 = vec![vec![TestField{0:3}], vec![TestField{0:4}]];
        matrix1.block_sum(&mut matrix2);
        let expected = vec![
            vec![TestField{0:1}, TestField{0:2}, TestField{0:0}],
            vec![TestField{0:0}, TestField{0:0}, TestField{0:3}],
            vec![TestField{0:0}, TestField{0:0}, TestField{0:4}],
        ];
        assert_eq!(matrix1, expected);
    }
}
