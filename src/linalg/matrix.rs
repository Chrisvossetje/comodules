use crate::linalg::field::Field;

#[derive(Debug, Clone, PartialEq)]
pub struct FieldMatrix<F: Field> {
    pub data: Vec<Vec<F>>,
    pub domain: usize,
    pub codomain: usize,
}

pub trait Matrix<F: Field>: Clone {
    fn zero(domain: usize, codomain: usize) -> Self;
    fn identity(d: usize) -> Self;

    fn get(&self, x: usize, y: usize) -> F;
    fn set(&mut self, x: usize, y: usize, f: F);

    // domain l == codomain r, l \circ r
    fn compose(&self, rhs: &mut Self) -> Self;

    fn transpose(&self) -> Self;

    fn domain(&self) -> usize;
    fn codomain(&self) -> usize;

    // OR JUST A USIZE, NOT KNOWING WHICH COLUMN IT ORIGINATED FROM ?
    // nah, both
    fn pivots(&self) -> Vec<(usize, usize)>;

    fn vstack(&mut self, other: &mut Self);
    fn block_sum(&mut self, other: &mut Self);

    fn rref(&mut self);

    // Faster but requires self to be in rref ?
    fn rref_kernel(&self) -> Self;

    fn cokernel(&self) -> Self;
    fn kernel(&self) -> Self;

    fn first_non_zero_entry(&self) -> Option<(usize, usize)>;
}

impl<F: Field> Matrix<F> for FieldMatrix<F> {
    // TODO: THIS DOES NOT WORK
    // See comment at kernel
    fn cokernel(&self) -> Self {
        let mut trans = self.transpose();
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
        if self.domain == 0 {
            return Self {
                data: vec![],
                domain: 0,
                codomain: 0,
            }
        }
        // let rows = self.codomain;
        // let cols = self.domain;
        let mut pivots = vec![None; self.domain];
        let mut free_vars = Vec::new();

        let pivs = self.pivots();

        for (i, j) in pivs {
            pivots[i] = Some(j);
        }

        // Identify free variables (non-pivot columns)
        for j in 0..self.domain {
            if pivots[j].is_none() {
                free_vars.push(j);
            }
        }

        // Generate kernel basis vectors
        let mut kernel = Vec::new();
        for &free_var in &free_vars {
            let mut kernel_vector = vec![F::zero(); self.domain];
            kernel_vector[free_var] = F::one(); // Free variable set to 1

            // Back-substitute to calculate the dependent variables
            for i in 0..self.codomain {
                match self.data[i][free_var].inv() {
                    None => {}
                    Some(inv) => {
                        for j in 0..free_var {
                            match pivots[j] {
                                None => {}
                                Some(k) => {
                                    if i == k {
                                        kernel_vector[j] += inv;
                                        break;
                                    }
                                }
                            }
                        }
                        break;
                    }
                }
            }

            kernel.push(kernel_vector);
        }
        let codomain = kernel.len();
        Self {
            data: kernel,
            domain: self.domain,
            codomain,
        }
    }

    // fn transpose(&mut self) {
    //     let rows = self.codomain;
    //     let cols = self.domain;
    //     let mut new_matrix = vec![vec![F::zero(); rows]; cols];

    //     for i in 0..rows {
    //         for j in 0..cols {
    //             new_matrix[j][i] = self.data[i][j];
    //         }
    //     }

    //     self.data = new_matrix;
    //     self.codomain = cols;
    //     self.domain = rows;
    // }

    fn transpose(&self) -> Self {
        let rows = self.codomain;
        let cols = self.domain;
        let mut new_matrix = vec![vec![F::zero(); rows]; cols];

        for i in 0..rows {
            for j in 0..cols {
                new_matrix[j][i] = self.data[i][j];
            }
        }

        Self {
            data: new_matrix,
            domain: rows,
            codomain: cols,
        }
    }

    fn rref(&mut self) {
        let rows = self.codomain;
        let cols = self.domain;

        let mut lead = 0; // Index of the leading column

        for r in 0..rows {
            if lead >= cols {
                break;
            }

            let mut i = r;
            while self.data[i][lead].is_zero() {
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
            self.data.swap(i, r);

            // Normalize the pivot row (make the leading coefficient 1)
            let pivot = self.data[r][lead];
            if !pivot.is_zero() {
                let pivot_inv = pivot.inv().expect("Pivot should be invertible");
                for j in 0..cols {
                    self.data[r][j] *= pivot_inv;
                }
            }

            // Eliminate all other entries in the leading column
            for i in 0..rows {
                if i != r {
                    let factor = self.data[i][lead];
                    for j in 0..cols {
                        let temp = factor * self.data[r][j];
                        self.data[i][j] -= temp;
                    }
                }
            }

            lead += 1;
        }
    }

    // (Row, Column)
    fn pivots(&self) -> Vec<(usize, usize)> {
        let mut id = 0;
        let mut pivots = vec![];
        for i in 0..self.codomain {
            if !self.data[id][i].is_zero() {
                pivots.push((id, i));
                id += 1;
            }
            if id == self.data.len() {
                break;
            }
        }
        pivots
    }

    fn vstack(&mut self, other: &mut Self) {
        assert_eq!(
            self.domain(),
            other.domain(),
            "Domains of the two matrices do not have the same dimension"
        );

        self.data.append(&mut other.data);
        self.codomain += other.codomain;
    }

    fn block_sum(&mut self, other: &mut Self) {
        let self_domain = self.domain;
        let other_domain = other.domain;
        for el in self.data.iter_mut() {
            let mut zeros = vec![F::zero(); other_domain];
            el.append(&mut zeros);
        }

        for oth in other.data.iter_mut() {
            let mut zeros = vec![F::zero(); self_domain];
            zeros.append(oth);

            self.data.push(zeros);
        }

        self.domain += other_domain;
        self.codomain += other.codomain;
    }

    fn identity(d: usize) -> Self {
        let mut matrix = vec![vec![F::zero(); d]; d];
        for i in 0..d {
            matrix[i][i] = F::one();
        }
        Self {
            data: matrix,
            domain: d,
            codomain: d,
        }
    }

    fn get(&self, x: usize, y: usize) -> F {
        self.data[y][x]
    }

    fn set(&mut self, x: usize, y: usize, f: F) {
        self.data[y][x] = f;
    }

    fn domain(&self) -> usize {
        self.domain
    }
    fn codomain(&self) -> usize {
        self.codomain
    }

    // domain l == codomain r, l \circ r
    fn compose(&self, rhs: &mut Self) -> Self {
        assert_eq!(
            self.domain,
            rhs.codomain,
            "Matrix domain not equal to codomain"
        );

        let trans = rhs.transpose();
        let mut compose = vec![vec![F::zero(); trans.codomain]; self.codomain];

        for x in 0..self.codomain {
            for y in 0..trans.codomain {
                compose[x][y] = dot_product(&self.data[x], &trans.data[y])
            }
        }
        Self {
            data: compose,
            domain: trans.codomain,
            codomain: self.codomain,
        }
    }

    fn zero(domain: usize, codomain: usize) -> Self {
        Self {
            data: vec![vec![F::zero(); domain]; codomain],
            domain,
            codomain,
        }
    }

    // (Row, Column) == (codomain, domain)
    fn first_non_zero_entry(&self) -> Option<(usize, usize)> {
        for codom_id in 0..self.codomain {
            for dom_id in 0..self.domain {
                if !self.data[codom_id][dom_id].is_zero() {
                    return Some((codom_id, dom_id));
                }
            }
        }
        None
    }
}

fn dot_product<F: Field>(l: &Vec<F>, r: &Vec<F>) -> F {
    l.iter().zip(r.iter()).map(|(x, y)| *x * *y).sum()
}
