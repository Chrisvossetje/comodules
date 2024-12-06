#[cfg(test)]
mod tests {
    use crate::{linalg::field::{Field, F2}, linalg::matrix::{kernel_from_rref, row_reduce_to_rref}};
    
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
        row_reduce_to_rref(&mut matrix);
        println!("Rref Matrix:");
        for row in matrix.clone() {
            println!("{:?}", row);
        }
        let kern = kernel_from_rref(&matrix);

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

        row_reduce_to_rref(&mut matrix);

        println!("\nReduced Row Echelon Form:");
        for row in &matrix {
            println!("{:?}", row);
        }
    }
}