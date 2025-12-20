#[cfg(test)]
mod tests {
    use crate::{matrices::f2_matrix::F2Matrix, matrix::Matrix, ring::CRing, rings::finite_fields::F2};

    #[test]
    fn test_vstack_1() {
        let mut matrix1 = F2Matrix {
            data: vec![0b10], // F2(1), F2(0)
            domain: 2,
            codomain: 1,
        };
        matrix1.set_element(0, 0, F2::one());
        matrix1.set_element(1, 0, F2::zero());
        
        let mut matrix2 = F2Matrix::zero(2, 1);
        matrix2.set_element(0, 0, F2::zero());
        matrix2.set_element(1, 0, F2::one());
        
        matrix1.vstack(&mut matrix2);
        
        assert_eq!(matrix1.domain(), 2);
        assert_eq!(matrix1.codomain(), 2);
        assert_eq!(matrix1.get_element(0, 0), F2::one());
        assert_eq!(matrix1.get_element(1, 0), F2::zero());
        assert_eq!(matrix1.get_element(0, 1), F2::zero());
        assert_eq!(matrix1.get_element(1, 1), F2::one());
    }

    #[test]
    fn test_vstack_2() {
        let mut matrix1 = F2Matrix::zero(2, 2);
        matrix1.set_element(0, 0, F2::one());
        matrix1.set_element(1, 0, F2::one());
        matrix1.set_element(0, 1, F2::zero());
        matrix1.set_element(1, 1, F2::zero());
        
        let mut matrix2 = F2Matrix::zero(2, 1);
        matrix2.set_element(0, 0, F2::zero());
        matrix2.set_element(1, 0, F2::one());
        
        matrix1.vstack(&mut matrix2);
        
        assert_eq!(matrix1.domain(), 2);
        assert_eq!(matrix1.codomain(), 3);
        assert_eq!(matrix1.get_element(0, 0), F2::one());
        assert_eq!(matrix1.get_element(1, 0), F2::one());
        assert_eq!(matrix1.get_element(0, 1), F2::zero());
        assert_eq!(matrix1.get_element(1, 1), F2::zero());
        assert_eq!(matrix1.get_element(0, 2), F2::zero());
        assert_eq!(matrix1.get_element(1, 2), F2::one());
    }

    #[test]
    fn test_block_sum_1() {
        let mut matrix1 = F2Matrix::zero(1, 2);
        matrix1.set_element(0, 0, F2::one());
        matrix1.set_element(0, 1, F2::zero());
        
        let matrix2 = F2Matrix::zero(1, 2);
        let mut matrix2 = matrix2;
        matrix2.set_element(0, 0, F2::one());
        matrix2.set_element(0, 1, F2::zero());
        
        matrix1.block_sum(&matrix2);
        
        assert_eq!(matrix1.domain(), 2);
        assert_eq!(matrix1.codomain(), 4);
        assert_eq!(matrix1.get_element(0, 0), F2::one());
        assert_eq!(matrix1.get_element(1, 0), F2::zero());
        assert_eq!(matrix1.get_element(0, 1), F2::zero());
        assert_eq!(matrix1.get_element(1, 1), F2::zero());
        assert_eq!(matrix1.get_element(0, 2), F2::zero());
        assert_eq!(matrix1.get_element(1, 2), F2::one());
        assert_eq!(matrix1.get_element(0, 3), F2::zero());
        assert_eq!(matrix1.get_element(1, 3), F2::zero());
    }

    #[test]
    fn test_block_sum_2() {
        let mut matrix1 = F2Matrix::zero(2, 1);
        matrix1.set_element(0, 0, F2::one());
        matrix1.set_element(1, 0, F2::zero());
        
        let mut matrix2 = F2Matrix::zero(1, 2);
        matrix2.set_element(0, 0, F2::one());
        matrix2.set_element(0, 1, F2::zero());
        
        matrix1.block_sum(&matrix2);
        
        assert_eq!(matrix1.domain(), 3);
        assert_eq!(matrix1.codomain(), 3);
        assert_eq!(matrix1.get_element(0, 0), F2::one());
        assert_eq!(matrix1.get_element(1, 0), F2::zero());
        assert_eq!(matrix1.get_element(2, 0), F2::zero());
        assert_eq!(matrix1.get_element(0, 1), F2::zero());
        assert_eq!(matrix1.get_element(1, 1), F2::zero());
        assert_eq!(matrix1.get_element(2, 1), F2::one());
        assert_eq!(matrix1.get_element(0, 2), F2::zero());
        assert_eq!(matrix1.get_element(1, 2), F2::zero());
        assert_eq!(matrix1.get_element(2, 2), F2::zero());
    }

    #[test]
    fn test_get() {
        let mut matrix = F2Matrix::zero(2, 2);
        matrix.set_element(0, 0, F2::one());
        matrix.set_element(1, 0, F2::zero());
        matrix.set_element(0, 1, F2::one());
        matrix.set_element(1, 1, F2::zero());

        assert_eq!(matrix.get(0, 0), F2::one());
        assert_eq!(matrix.get(1, 0), F2::zero());
        assert_eq!(matrix.get(0, 1), F2::one());
        assert_eq!(matrix.get(1, 1), F2::zero());
    }

    #[test]
    fn test_set() {
        let mut matrix = F2Matrix::zero(2, 2);

        matrix.set(0, 0, F2::one());
        matrix.set(1, 1, F2::one());

        assert_eq!(matrix.get(0, 0), F2::one());
        assert_eq!(matrix.get(1, 1), F2::one());
        assert_eq!(matrix.get(0, 1), F2::zero());
        assert_eq!(matrix.get(1, 0), F2::zero());
    }

    #[test]
    fn test_add_at() {
        let mut matrix = F2Matrix::zero(2, 2);
        matrix.set_element(0, 0, F2::one());
        matrix.set_element(1, 0, F2::zero());
        matrix.set_element(0, 1, F2::one());
        matrix.set_element(1, 1, F2::zero());

        matrix.add_at(0, 0, F2::one()); // 1 + 1 = 0 in F2
        matrix.add_at(1, 1, F2::one()); // 0 + 1 = 1 in F2

        assert_eq!(matrix.get(0, 0), F2::zero());
        assert_eq!(matrix.get(1, 1), F2::one());
    }

    #[test]
    fn test_set_row() {
        let mut matrix = F2Matrix::zero(2, 2);

        let new_row = vec![F2::one(), F2::zero()];
        matrix.set_row(0, &new_row);

        assert_eq!(matrix.get(0, 0), F2::one());
        assert_eq!(matrix.get(1, 0), F2::zero());
    }

    #[test]
    fn test_domain_codomain() {
        let matrix = F2Matrix::zero(2, 3);

        assert_eq!(matrix.domain(), 2);
        assert_eq!(matrix.codomain(), 3);
    }

    #[test]
    fn test_compose() {
        // Create a 2x2 matrix A
        let mut matrix_a = F2Matrix::zero(2, 2);
        matrix_a.set_element(0, 0, F2::one());
        matrix_a.set_element(1, 0, F2::zero());
        matrix_a.set_element(0, 1, F2::one());
        matrix_a.set_element(1, 1, F2::one());

        // Create a 2x1 matrix B
        let mut matrix_b = F2Matrix::zero(2, 1);
        matrix_b.set_element(0, 0, F2::one());
        matrix_b.set_element(1, 0, F2::one());

        // Compute A * B
        let result = matrix_a.compose(&matrix_b);

        assert_eq!(result.domain(), 2);
        assert_eq!(result.codomain(), 2);
        // First row: [1, 0] * [1, 1]^T = 1*1 + 0*1 = 1
        assert_eq!(result.get_element(0, 0), F2::one());
        // Second row: [1, 1] * [1, 1]^T = 1*1 + 1*1 = 0 (in F2)
        assert_eq!(result.get_element(1, 0), F2::zero());
    }

    #[test]
    fn test_transpose() {
        let mut matrix = F2Matrix::zero(2, 3);
        matrix.set_element(0, 0, F2::one());
        matrix.set_element(1, 0, F2::zero());
        matrix.set_element(0, 1, F2::one());
        matrix.set_element(1, 1, F2::one());
        matrix.set_element(0, 2, F2::zero());
        matrix.set_element(1, 2, F2::one());

        let transposed = matrix.transpose();

        assert_eq!(transposed.domain(), 3);
        assert_eq!(transposed.codomain(), 2);
        assert_eq!(transposed.get_element(0, 0), F2::one());
        assert_eq!(transposed.get_element(0, 1), F2::zero());
        assert_eq!(transposed.get_element(1, 0), F2::one());
        assert_eq!(transposed.get_element(1, 1), F2::one());
        assert_eq!(transposed.get_element(2, 0), F2::zero());
        assert_eq!(transposed.get_element(2, 1), F2::one());
    }

    #[test]
    fn test_extend_one_row() {
        let mut matrix = F2Matrix::zero(2, 2);
        matrix.set_element(0, 0, F2::one());
        matrix.set_element(1, 1, F2::one());

        matrix.extend_one_row();

        assert_eq!(matrix.domain(), 2);
        assert_eq!(matrix.codomain(), 3);
        assert_eq!(matrix.get_element(0, 0), F2::one());
        assert_eq!(matrix.get_element(1, 1), F2::one());
        assert_eq!(matrix.get_element(0, 2), F2::zero());
        assert_eq!(matrix.get_element(1, 2), F2::zero());
    }

    #[test]
    fn test_large_matrix() {
        // Test with a matrix larger than 64 bits to verify bit packing
        let mut matrix = F2Matrix::zero(100, 50);
        
        // Set some elements
        matrix.set_element(0, 0, F2::one());
        matrix.set_element(63, 10, F2::one());
        matrix.set_element(64, 20, F2::one());
        matrix.set_element(99, 49, F2::one());

        // Check they were set correctly
        assert_eq!(matrix.get_element(0, 0), F2::one());
        assert_eq!(matrix.get_element(63, 10), F2::one());
        assert_eq!(matrix.get_element(64, 20), F2::one());
        assert_eq!(matrix.get_element(99, 49), F2::one());
        
        // Check other elements are zero
        assert_eq!(matrix.get_element(1, 0), F2::zero());
        assert_eq!(matrix.get_element(62, 10), F2::zero());
        assert_eq!(matrix.get_element(65, 20), F2::zero());
        assert_eq!(matrix.get_element(98, 49), F2::zero());
    }
}