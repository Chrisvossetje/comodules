use crate::fp::Field;


pub trait Matrix<F: Field> {
    fn kernel(&self);

    // Faster but requires self to be in rref ?
    fn rref_kernel(&self);
    
    fn transpose(&mut self);
    fn rref(&mut self);

    // OR JUST A USIZE, NOT KNOWING WHICH COLUMN IT ORIGINATED FROM ?
    fn pivots(&mut self) -> Vec<(usize,usize)>; 

    fn vstack(&mut self, other: &Self);
    fn block_sum(&mut self, other: &Self);

    fn identity(d: usize) -> Self;
}


pub struct FieldMatrix<F: Field> {

}


fn row_reduce_to_rref<F: Field>(matrix: &mut Vec<Vec<F>>) {
    let rows = matrix.len();
    let cols = matrix[0].len();

    for i in 0..rows {
        // Step 0: Ensure the pivot is nonzero by swapping rows if needed
        if matrix[i][i].is_zero() {
            for k in (i + 1)..rows {
                if !matrix[k][i].is_zero() {
                    matrix.swap(i, k); // Swap rows
                    break;
                }
            }
        }

        // Step 1: Normalize the pivot row
        let pivot = matrix[i][i];
        if !pivot.is_zero() {
            let inv = pivot.inv();
            for j in 0..cols {
                matrix[i][j] *= inv;
            }
        }

        // Step 2: Eliminate all entries below the pivot
        for k in (i + 1)..rows {
            let factor = matrix[k][i];
            for j in 0..cols {
                let res = factor * matrix[i][j];
                matrix[k][j] -= res;
            }
        }
    }

    let mut id = 0;
    let mut pivots = vec![];
    for i in (0..cols) {
        if !matrix[id][i].is_zero() {
            pivots.push((id, i));
            id += 1;
        }
        if id == rows {
            break;
        }
    }

    // Step 3: Eliminate entries above each pivot
    for (id, i) in pivots {
        for k in (0..id) {
            let factor = matrix[k][i];
            for j in (i..cols) {
                let res = factor * matrix[id][j];
                matrix[k][j] -= res;

            }
        }
    }
}

fn kernel_from_rref<F: Field>(matrix: &Vec<Vec<F>>) -> Vec<Vec<F>> {
    let rows = matrix.len();
    let cols = matrix[0].len();
    let mut pivots = vec![None; cols]; // Track pivot columns (exclude last column for augmented matrix)
    let mut free_vars = Vec::new();
    
    // Identify pivot columns
    for i in 0..rows {
        for j in 0..cols - 1 {
            if matrix[i][j] == F::one() && pivots[j].is_none() {
                pivots[j] = Some(i);
                break;
            }
        }
    }

    // Identify free variables (non-pivot columns)
    for j in 0..cols {
        if pivots[j].is_none() {
            free_vars.push(j);
        }
    }

    println!("{:?}", pivots);
    println!("{:?}", free_vars);

    // Generate kernel basis vectors
    let mut kernel = Vec::new();
    for &free_var in &free_vars {
        let mut kernel_vector = vec![F::zero(); cols];
        kernel_vector[free_var] = F::one(); // Free variable set to 1

        // Back-substitute to calculate the dependent variables
        for i in 0..rows {
            if !matrix[i][free_var].is_zero() {
                for j in 0..free_var {
                    match pivots[j] {
                        None => {},
                        Some(k) => {
                            if i == k {
                                let res =  matrix[i][free_var] * kernel_vector[free_var];
                                kernel_vector[j] -= res;
                                break;
                            }
                        }
                    }
                }
                break;
            }
        }


        kernel.push(kernel_vector);
    }

    kernel
}


#[cfg(test)]
mod tests {
    use crate::{fp::{Field, F2}, matrix::{kernel_from_rref, row_reduce_to_rref}};
    
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
        vec![F2::one(), F2::one(), F2::zero(), F2::one()],
        vec![F2::one(), F2::one(), F2::zero(), F2::zero()],
        vec![F2::one(), F2::zero(), F2::zero(), F2::one()],
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