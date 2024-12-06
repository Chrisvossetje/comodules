use crate::linalg::field::{Field};

pub type FieldMatrix<F: Field> = Vec<Vec<F>>;

pub trait Matrix<F: Field> : Clone {
    fn zero(domain: usize, codomain: usize) -> Self;
    fn identity(d: usize) -> Self;
    
    fn get(&self, x: usize, y: usize) -> F;
    fn set(&mut self, x: usize, y: usize, f: F);
    
    // domain l == codomain r, l \circ r
    fn compose(self, rhs: &mut Self) -> Self;

    fn transpose(&mut self);
    fn get_transpose(&self) -> Self;

    fn domain(&self) -> usize;
    fn codomain(&self) -> usize;

    // OR JUST A USIZE, NOT KNOWING WHICH COLUMN IT ORIGINATED FROM ?
    // nah, both
    fn pivots(&self) -> Vec<(usize,usize)>; 

    fn vstack(&mut self, other: &mut Self);
    fn block_sum(&mut self, other: &mut Self);

    
    fn rref(&mut self);
    
    // Faster but requires self to be in rref ?
    fn rref_kernel(&self) -> Self;
    
    fn cokernel(&self) -> Self;
    fn kernel(&self) -> Self;
}

impl<F: Field> Matrix<F> for Vec<Vec<F>> {
    // TODO: THIS DOES NOT WORK
    // See comment at kernel
    fn cokernel(&self) -> Self {
        let mut trans = self.get_transpose();
        trans.rref();
        trans.rref_kernel()
    }

    // TODO: THIS DOES NOT WORK
    // This gives back the kernel of rref of self
    // rref should include the translation to make this work ?
    fn kernel(&self) -> Self {
        let mut clone = self.clone();
        clone.rref();
        clone.rref_kernel()
    }

    fn rref_kernel(&self) -> Self {
        let rows = self.len();
        let cols = self[0].len();
        let mut pivots = vec![None; cols]; // Track pivot columns (exclude last column for augmented matrix)
        let mut free_vars = Vec::new();
        
        // Identify pivot columns
        for i in 0..rows {
            for j in 0..cols - 1 {
                if self[i][j] == F::one() && pivots[j].is_none() {
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
                match self[i][free_var].inv() {
                    None => {},
                    Some(inv) => {
                        for j in 0..free_var {
                            match pivots[j] {
                                None => {},
                                Some(k) => {
                                    if i == k {
                                        kernel_vector[j] += inv;
                                        break;
                                    }
                                }
                            }
                        }
                        break;
                    },
                }
            }


            kernel.push(kernel_vector);
        }

        kernel
    }

    fn transpose(&mut self) {
        let rows = self.len();
        let cols = self[0].len();
        let mut new_matrix = vec![vec![F::zero(); rows]; cols];

        for i in 0..rows {
            for j in 0..cols {
                new_matrix[j][i] = self[i][j];
            }
        }

        *self = new_matrix;
    }

    fn get_transpose(&self) -> Self {
        let rows = self.len();
        let cols = self[0].len();
        let mut new_matrix = vec![vec![F::zero(); rows]; cols];

        for i in 0..rows {
            for j in 0..cols {
                new_matrix[j][i] = self[i][j];
            }
        }

        new_matrix
    }
    
    fn rref(&mut self) {
        let rows = self.len();
        let cols = if rows > 0 { self[0].len() } else { 0 };
    
        let mut lead = 0; // Index of the leading column
    
        for r in 0..rows {
            if lead >= cols {
                break;
            }
    
            let mut i = r;
            while self[i][lead].is_zero() {
                i += 1;
                if i == rows {
                    i = r;
                    lead += 1;
                    if lead == cols {
                        return;
                    }
                }
            }
    
            // Swap rows to move the pivot row to the current row
            self.swap(i, r);
    
            // Normalize the pivot row (make the leading coefficient 1)
            let pivot = self[r][lead];
            if !pivot.is_zero() {
                let pivot_inv = pivot.inv().expect("Pivot should be invertible");
                for j in 0..cols {
                    self[r][j] *= pivot_inv;
                }
            }
    
            // Eliminate all other entries in the leading column
            for i in 0..rows {
                if i != r {
                    let factor = self[i][lead];
                    for j in 0..cols {
                        let temp = factor * self[r][j];
                        self[i][j] -= temp;
                    }
                }
            }
    
            lead += 1;
        }
    }

    fn pivots(&self) -> Vec<(usize,usize)> {
        let mut id = 0;
        let mut pivots = vec![];
        for i in (0..self[0].len()) {
            if !self[id][i].is_zero() {
                pivots.push((id, i));
                id += 1;
            }
            if id == self.len() {
                break;
            }
        }
        pivots
    }

    fn vstack(&mut self, other: &mut Self) {
        assert_eq!(self.domain(), other.domain(), "Domains of the two matrices do not have the same dimension");

        self.append(other);
    }
    
    fn block_sum(&mut self, other: &mut Self) {
        if other[0].len() == 0 {
            return;
        }

        let self_domain = self[0].len();
        let other_domain = other[0].len();
        for el in self.iter_mut() {
            let mut zeros = vec![F::zero(); other_domain]; 
            el.append(&mut zeros);
        };

        for oth in other.iter_mut() {
            let mut zeros = vec![F::zero(); self_domain];
            zeros.append(oth);
            
            self.push(zeros);
        }
    }

    fn identity(d: usize) -> Self {
        let mut matrix = vec![vec![F::zero(); d]; d];
        for i in 0..d {
            matrix[i][i] = F::one();
        }
        matrix
    }
    
    fn get(&self, x: usize, y: usize) -> F {
        self[y][x]
    }
    
    fn set(&mut self, x: usize, y: usize, f: F) {
        self[y][x] = f;
    }
    
    fn domain(&self) -> usize {
        self.len()
    }
    fn codomain(&self) -> usize {
        self.len()
    }

    // domain l == codomain r, l \circ r
    fn compose(self, rhs: &mut Self) -> Self {
        assert_eq!(self.domain(), rhs.codomain(), "Matrix domain not equal to codomain");
        
        rhs.transpose();
        let mut compose = vec![vec![F::zero(); self.codomain()]; rhs.domain()];

        for x in (0..self.codomain()) {
            for y in (0..rhs.domain()) {
                compose[x][y] = dot_product(&self[x] , &rhs[y])
            }
        }

        compose
    }
    
    fn zero(domain: usize, codomain: usize) -> Self {
        vec![vec![F::zero(); domain]; codomain]
    }
    
    
}

fn dot_product<F: Field>(l: &Vec<F>, r: &Vec<F>) -> F {
    l.iter().zip(r.iter()).map(|(x,y)| *x * *y).sum()
}