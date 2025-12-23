#[cfg(test)]
mod tests {
    use crate::{matrices::{f2_matrix::F2Matrix, flat_matrix::FlatMatrix}, matrix::Matrix, ring::CRing, rings::finite_fields::F2};

    #[test]
    fn test_rref_one_vector() {
        let mut matrix = F2Matrix::zero(1, 3);
        // Create matrix:
        // [1]
        // [1]
        // [1]
        matrix.set_element(0, 0, F2::one()); 
        matrix.set_element(0, 1, F2::one()); 
        matrix.set_element(0, 2, F2::one()); 

        matrix.rref();


        // Should be in rref form:
        // [1 0 1]
        // [0 1 1]
        assert_eq!(matrix.get_element(0, 0), F2::one());
        assert_eq!(matrix.get_element(0, 1), F2::zero());
        assert_eq!(matrix.get_element(0, 2), F2::zero());

        assert!(matrix.is_rref());
    }

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

        // Should be in rref form:
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

        println!("{:?}", matrix);


        assert_eq!(matrix.rank(), 2); // Full rank
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
        assert_eq!(kernel.domain(), 3);
        assert_eq!(kernel.codomain(), 1);

        // The kernel vector should be [1, 1, 1]^T (since [1 0 1; 0 1 1] * [1; 1; 1] = [0; 0])
        assert_eq!(kernel.get_element(0, 0), F2::one());
        assert_eq!(kernel.get_element(1, 0), F2::one());
        assert_eq!(kernel.get_element(2, 0), F2::one());

        // Verify it's actually in the kernel
        let result = matrix.compose(&kernel.transpose());
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
    fn test_is_rref() {
        // Test a matrix that is in rref
        let mut rref_matrix = F2Matrix::zero(3, 2);
        rref_matrix.set_element(0, 0, F2::one());  // pivot at (0,0)
        rref_matrix.set_element(2, 0, F2::one());
        rref_matrix.set_element(1, 1, F2::one());  // pivot at (1,1)
        rref_matrix.set_element(2, 1, F2::zero());

        assert!(rref_matrix.is_rref());

        // Test a matrix that is NOT in rref
        let mut non_rref_matrix = F2Matrix::zero(3, 2);
        non_rref_matrix.set_element(0, 0, F2::one());
        non_rref_matrix.set_element(2, 0, F2::one());
        non_rref_matrix.set_element(0, 1, F2::one()); // Should be 0 for rref
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
        assert_eq!(matrix.first_non_zero_entry(), Some((1, 2)));
        
        // Add another non-zero entry that comes earlier
        matrix.set_element(0, 1, F2::one());
        assert_eq!(matrix.first_non_zero_entry(), Some((0, 1)));
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

    #[test]
    fn test_rref_4() {
        let mut matrix = FlatMatrix {
            data: vec![F2 { 0: 1 }, F2 { 0: 0 }],
            domain: 1,
            codomain: 2,
        };
        let expected = FlatMatrix {
            data: vec![F2 { 0: 1 }, F2 { 0: 0 }],
            domain: 1,
            codomain: 2,
        };
        let mut f2_matrix = F2Matrix::from_flat(matrix);
        f2_matrix.rref();
        let f2_expected = F2Matrix::from_flat(expected);
        assert_eq!(f2_matrix, f2_expected);
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
        let mut f2_matrix = F2Matrix::from_flat(matrix);
        f2_matrix.rref();
        let f2_expected = F2Matrix::from_flat(expected);
        assert_eq!(f2_matrix, f2_expected);
    }
}