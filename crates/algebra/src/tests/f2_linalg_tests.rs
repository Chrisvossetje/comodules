#[cfg(test)]
mod tests {
    use crate::{matrices::f2_matrix::F2Matrix, matrix::Matrix, ring::CRing, rings::finite_fields::F2};

    #[test]
    fn test_rref_simple() {
        let mut matrix = F2Matrix::zero(3, 2);
        // Create matrix:
        // [1 0 1]
        // [0 1 1] 
        matrix.set_element(0, 0, F2::one());  // (0,0) = 1
        matrix.set_element(2, 0, F2::one());  // (2,0) = 1
        matrix.set_element(1, 1, F2::one());  // (1,1) = 1  
        matrix.set_element(2, 1, F2::one());  // (2,1) = 1

        matrix.rref();

        // Should be in RREF form:
        // [1 0 1]
        // [0 1 1]
        assert_eq!(matrix.get_element(0, 0), F2::one());
        assert_eq!(matrix.get_element(1, 0), F2::zero());
        assert_eq!(matrix.get_element(2, 0), F2::one());
        assert_eq!(matrix.get_element(0, 1), F2::zero());
        assert_eq!(matrix.get_element(1, 1), F2::one());
        assert_eq!(matrix.get_element(2, 1), F2::one());

        assert!(matrix.is_rref());
    }

    #[test]
    fn test_rref_with_row_swaps() {
        let mut matrix = F2Matrix::zero(3, 3);
        // Create matrix:
        // [0 1 1]
        // [1 0 1]
        // [1 1 0]
        matrix.set_element(1, 0, F2::one());  // (1,0) = 1
        matrix.set_element(2, 0, F2::one());  // (2,0) = 1
        matrix.set_element(0, 1, F2::one());  // (0,1) = 1
        matrix.set_element(2, 1, F2::one());  // (2,1) = 1
        matrix.set_element(0, 2, F2::one());  // (0,2) = 1
        matrix.set_element(1, 2, F2::one());  // (1,2) = 1

        matrix.rref();

        assert!(matrix.is_rref());
        assert_eq!(matrix.rank(), 3); // Full rank
    }

    #[test]
    fn test_pivots() {
        let mut matrix = F2Matrix::zero(4, 3);
        // Create matrix:
        // [1 0 0 1]
        // [0 1 0 1]
        // [0 0 1 1]
        matrix.set_element(0, 0, F2::one());
        matrix.set_element(3, 0, F2::one());
        matrix.set_element(1, 1, F2::one());
        matrix.set_element(3, 1, F2::one());
        matrix.set_element(2, 2, F2::one());
        matrix.set_element(3, 2, F2::one());

        let pivots = matrix.pivots();
        assert_eq!(pivots, vec![(0, 0), (1, 1), (2, 2)]);
    }

    #[test]
    fn test_rank_and_nullity() {
        let mut matrix = F2Matrix::zero(3, 2);
        // Create matrix:
        // [1 0 1]
        // [0 1 1]
        matrix.set_element(0, 0, F2::one());
        matrix.set_element(2, 0, F2::one());
        matrix.set_element(1, 1, F2::one());
        matrix.set_element(2, 1, F2::one());

        assert_eq!(matrix.rank(), 2);
        assert_eq!(matrix.nullity(), 1);
    }

    #[test]
    fn test_kernel() {
        let mut matrix = F2Matrix::zero(3, 2);
        // Create matrix:
        // [1 0 1]
        // [0 1 1]
        matrix.set_element(0, 0, F2::one());
        matrix.set_element(2, 0, F2::one());
        matrix.set_element(1, 1, F2::one());
        matrix.set_element(2, 1, F2::one());
        
        matrix.rref();
        let kernel = matrix.rref_kernel();

        // Kernel should be 1-dimensional
        assert_eq!(kernel.codomain(), 1);
        assert_eq!(kernel.domain(), 3);

        // The kernel vector should be [1, 1, 1]^T (since [1 0 1; 0 1 1] * [1; 1; 1] = [0; 0])
        assert_eq!(kernel.get_element(0, 0), F2::one());
        assert_eq!(kernel.get_element(1, 0), F2::one());
        assert_eq!(kernel.get_element(2, 0), F2::one());

        // Verify it's actually in the kernel
        let result = matrix.compose(&kernel);
        for i in 0..result.domain() {
            for j in 0..result.codomain() {
                assert_eq!(result.get_element(i, j), F2::zero());
            }
        }
    }

    #[test]
    fn test_swap_rows() {
        let mut matrix = F2Matrix::zero(3, 3);
        matrix.set_element(0, 0, F2::one());
        matrix.set_element(1, 0, F2::zero());
        matrix.set_element(2, 0, F2::one());
        matrix.set_element(0, 1, F2::zero());
        matrix.set_element(1, 1, F2::one());
        matrix.set_element(2, 1, F2::zero());

        matrix.swap_rows(0, 1);

        // Row 0 and row 1 should be swapped
        assert_eq!(matrix.get_element(0, 0), F2::zero());
        assert_eq!(matrix.get_element(1, 0), F2::one());
        assert_eq!(matrix.get_element(2, 0), F2::zero());
        assert_eq!(matrix.get_element(0, 1), F2::one());
        assert_eq!(matrix.get_element(1, 1), F2::zero());
        assert_eq!(matrix.get_element(2, 1), F2::one());
    }

    #[test]
    fn test_add_row_to_row() {
        let mut matrix = F2Matrix::zero(3, 2);
        matrix.set_element(0, 0, F2::one());
        matrix.set_element(1, 0, F2::zero());
        matrix.set_element(2, 0, F2::one());
        matrix.set_element(0, 1, F2::zero());
        matrix.set_element(1, 1, F2::one());
        matrix.set_element(2, 1, F2::one());

        matrix.add_row_to_row(0, 1);

        // Row 1 = Row 1 XOR Row 0
        assert_eq!(matrix.get_element(0, 0), F2::one());   // unchanged
        assert_eq!(matrix.get_element(1, 0), F2::zero());  // unchanged
        assert_eq!(matrix.get_element(2, 0), F2::one());   // unchanged
        assert_eq!(matrix.get_element(0, 1), F2::one());   // 0 XOR 1 = 1
        assert_eq!(matrix.get_element(1, 1), F2::one());   // 1 XOR 0 = 1
        assert_eq!(matrix.get_element(2, 1), F2::zero());  // 1 XOR 1 = 0
    }

    #[test]
    fn test_solve_system_unique_solution() {
        let mut a = F2Matrix::zero(3, 3);
        // Create identity matrix
        a.set_element(0, 0, F2::one());
        a.set_element(1, 1, F2::one());
        a.set_element(2, 2, F2::one());

        let mut b = F2Matrix::zero(1, 3);
        b.set_element(0, 0, F2::one());
        b.set_element(0, 1, F2::zero());
        b.set_element(0, 2, F2::one());

        let solution = a.solve(&b).expect("Should have a solution");

        assert_eq!(solution.get_element(0, 0), F2::one());
        assert_eq!(solution.get_element(1, 0), F2::zero());
        assert_eq!(solution.get_element(2, 0), F2::one());

        // Verify solution: A * x = b
        let result = a.compose(&solution);
        for i in 0..b.domain() {
            for j in 0..b.codomain() {
                assert_eq!(result.get_element(i, j), b.get_element(i, j));
            }
        }
    }

    #[test]
    fn test_solve_system_no_solution() {
        let mut a = F2Matrix::zero(2, 2);
        // Create matrix [1 1; 1 1] (rank 1)
        a.set_element(0, 0, F2::one());
        a.set_element(1, 0, F2::one());
        a.set_element(0, 1, F2::one());
        a.set_element(1, 1, F2::one());

        let mut b = F2Matrix::zero(1, 2);
        // b = [1; 0] - inconsistent with rank-1 matrix
        b.set_element(0, 0, F2::one());
        b.set_element(0, 1, F2::zero());

        let solution = a.solve(&b);
        assert!(solution.is_none());
    }

    #[test]
    fn test_is_rref() {
        // Test a matrix that is in RREF
        let mut rref_matrix = F2Matrix::zero(3, 2);
        rref_matrix.set_element(0, 0, F2::one());  // pivot at (0,0)
        rref_matrix.set_element(2, 0, F2::one());
        rref_matrix.set_element(1, 1, F2::one());  // pivot at (1,1)
        rref_matrix.set_element(2, 1, F2::zero());

        assert!(rref_matrix.is_rref());

        // Test a matrix that is NOT in RREF
        let mut non_rref_matrix = F2Matrix::zero(3, 2);
        non_rref_matrix.set_element(0, 0, F2::one());
        non_rref_matrix.set_element(2, 0, F2::one());
        non_rref_matrix.set_element(0, 1, F2::one()); // Should be 0 for RREF
        non_rref_matrix.set_element(1, 1, F2::one());
        non_rref_matrix.set_element(2, 1, F2::one());

        assert!(!non_rref_matrix.is_rref());
    }

    #[test]
    fn test_first_non_zero_entry() {
        let mut matrix = F2Matrix::zero(3, 3);
        
        // All zeros
        assert_eq!(matrix.first_non_zero_entry(), None);
        
        // Add a non-zero entry
        matrix.set_element(1, 2, F2::one());
        assert_eq!(matrix.first_non_zero_entry(), Some((2, 1)));
        
        // Add another non-zero entry that comes earlier
        matrix.set_element(0, 1, F2::one());
        assert_eq!(matrix.first_non_zero_entry(), Some((1, 0)));
    }

    #[test]
    fn test_large_matrix_rref() {
        // Test with a larger matrix to verify bit-packing works correctly
        let mut matrix = F2Matrix::zero(100, 50);
        
        // Create an identity-like matrix
        for i in 0..50 {
            matrix.set_element(i, i, F2::one());
        }
        
        // Add some additional entries
        for i in 50..100 {
            matrix.set_element(i, i % 50, F2::one());
        }

        matrix.rref();
        
        assert!(matrix.is_rref());
        assert_eq!(matrix.rank(), 50);
        assert_eq!(matrix.nullity(), 50);
    }
}