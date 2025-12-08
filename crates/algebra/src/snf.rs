use crate::{matrices::flat_matrix::FlatMatrix, matrix::Matrix, ring::{CRing, ValuationRing}};


pub(crate) trait SmithNormalForm<R: CRing>: Matrix<R> + Clone {
    /// snf(A) -> (U,S,V) 
    /// U A V = S
    fn snf(&self) -> (Self, Self, Self) {
        let (u,s,v, _, _) = self.full_snf();
        (u,s,v)
    }

    /// full_snf(A) -> (U,S,V,U^-1, V^-1) 
    /// U A V = S
    fn full_snf(&self) -> (Self, Self, Self, Self, Self);
}

impl<R: ValuationRing> SmithNormalForm<R> for FlatMatrix<R> {
    fn full_snf(&self) -> (Self, Self, Self, Self, Self) {
        enum Action<R: ValuationRing> {
            Swap(usize, usize),
            Add(usize, usize, R)
        }

        
        let mut s = self.clone();
        let mut u = Self::identity(self.codomain());
        let mut v = Self::identity(self.domain());

        let mut u_inv_actions: Vec<Action<R>> = vec![];
        let mut v_inv_actions = vec![];
        
        let m = s.codomain();
        let n = s.domain();
        
        let min_r = m.min(n);

        for r in 0..min_r {

            // Find element in (sub)matrix which divides all others
            let mut candidate = (r,r);
            let mut candidate_val = s.get(r, r); 
            for y in r..m {
                for x in r..n {
                    let el = s.get(x, y);
                    if !candidate_val.divides(&el) {
                        candidate_val = el;
                        candidate = (x,y);
                    } 
                }
            }

            if candidate_val.is_zero() {
                break;
            }

            // swap rows and colums to get candidate in correct position
            u.swap_rows(r, candidate.1);
            u_inv_actions.push(Action::Swap(r, candidate.1));
            s.swap_rows(r, candidate.1);

            s.swap_cols(r, candidate.0);
            v.swap_cols(r, candidate.0);
            v_inv_actions.push(Action::Swap(r, candidate.0));


            
            let pivot = candidate_val;

            // reduce rows below (r,r)
            for k in (r + 1)..m {
                let entry = s.get(r, k);
                if !entry.is_zero() {
                    let factor = -entry.unsafe_divide(pivot);

                    u.add_row_multiple(k, r, factor);
                    u_inv_actions.push(Action::Add(k,r,factor));
                    s.add_row_multiple(k, r, factor);
                }
            }

            // reduce columns next to (r,r)
            for l in (r+1)..s.domain {
                let entry = s.get(l, r);
                if !entry.is_zero() {
                    let factor = -entry.unsafe_divide(pivot);
                    s.add_col_multiple(l, r, factor);
                    v.add_col_multiple(l, r, factor);
                    v_inv_actions.push(Action::Add(l,r,factor));
                }
            }
        }
        
        // construct u and v inverse
        let mut u_inv = Self::identity(self.codomain());
        let mut v_inv = Self::identity(self.domain());
        
        for a in u_inv_actions.iter().rev() {
            match a {
                Action::Swap(b, c) => {
                    u_inv.swap_rows(*b, *c);
                },
                Action::Add(b, c, factor) => {
                    u_inv.add_row_multiple(*b, *c, factor.neg());
                },
            }
        }
        for a in v_inv_actions.iter().rev() {
            match a {
                Action::Swap(b, c) => {
                    v_inv.swap_cols(*b, *c);
                },
                Action::Add(b, c, factor) => {
                    v_inv.add_col_multiple(*b, *c, factor.neg());
                },
            }
        }

        (u, s, v, u_inv, v_inv)
    }
}

// TODO: THIS CAN BE REMOVED / IMPROVED TO WORK WITH PID's
// PROBABLY ONLY RELEVANT FOR THE PID Z (the integers) (or maybe a generic k[x])
// k[x,y] will require a different method (groebner :(, maybe not tho !)

// impl<F: Field> SmithNormalForm<UniPolRing<F>> for FlatMatrix<UniPolRing<F>> {    
//     fn snf(&self) -> (Self, Self, Self) 
//     where 
//         Self: Clone,
//     {
//         let mut s = self.clone();
//         let mut u = Self::identity(self.codomain());
//         let mut v = Self::identity(self.domain());
        
//         let m = s.codomain();
//         let n = s.domain();
//         // let mut column_indices = Vec::new();
//         let mut last_j = 0;
        
//         // Main algorithm: iterate t from 1 to m
//         for t in 0..m {
//             // println!("{}:::!!!\n",t);
//             // print_matrix(&s);

//             // Step I: Choosing a pivot
//             // Find smallest column index with non-zero entry starting from last_j            
//             let mut search = None;
//             'outer: for j in last_j..n { // column is choice of domain
//                 for a in t..m { // row is choice of codomain
//                     if !s.get(j, a).is_zero() {
//                         search = Some((a, j));
//                         break 'outer;
//                     }
//                 }
//             }

//             let (a, j_t) = match search {
//                 Some((a, j)) => (a, j),
//                 None => break, // No more non-zero entries
//             };
//             last_j = j_t + 1;
            
//             u.swap_rows(t, a);
//             s.swap_rows(t, a);

//             // As all "prime" factors are finite, this loop is finite
//             loop {
//                 // println!("Loop");
//                 // print_matrix(&s);

//                 let mut improved = false;
                
                
//                 // FIRST DO ROWS
                
//                 // Step II: Improving the pivot
//                 // Check if pivot divides all entries in column j_t
//                 let mut smallest_index = t;
//                 for k in (t + 1)..m { // rows
//                     if !s.get(j_t, smallest_index).divides(&s.get(j_t, k)) {
//                         smallest_index = k;
//                     }
//                 }

//                 u.swap_rows(t, smallest_index);
//                 s.swap_rows(t, smallest_index);
//                 if t != smallest_index {
//                     improved = true;
//                 }

//                 // println!("Row Step II");
//                 // print_matrix(&s);


//                 // Step III: Eliminating entries
//                 // Clear column j_t below position (t, j_t)
//                 let pivot = s.get(j_t, t);
//                 for k in (t + 1)..m {
//                     let entry = s.get(j_t, k);
//                     if !entry.is_zero() {
//                         let factor = -entry.unsafe_divide(pivot);
//                         // TODO
//                         u.add_row_multiple(k, t, factor);
//                         s.add_row_multiple(k, t, factor);
//                     }
//                 }
//                 // println!("Row Step III");
//                 // print_matrix(&s);


//                 // SECOND DO COLUMNS

//                 // Step II: Improving the pivot
//                 // Check if pivot divides all entries in column j_t
//                 let mut smallest_index = j_t;
//                 for k in (j_t + 1)..n {
//                     if !s.get(smallest_index, t).divides(&s.get(k, t)) {
//                         smallest_index = k;
//                     }
//                 }

                
//                 s.swap_cols(j_t, smallest_index);
//                 v.swap_cols(j_t, smallest_index);
//                 if j_t != smallest_index {
//                     improved = true;
//                 }
                
//                 // println!("Column Step II");
//                 // print_matrix(&s);


//                 // Step III: Eliminating entries
//                 // Clear column j_t below position (t, j_t)

//                 let pivot = s.get(j_t, t);
//                 for k in (j_t + 1)..n {
//                     let entry = s.get(k, t);
//                     if !entry.is_zero() {
//                         let factor = -entry.unsafe_divide(pivot);
//                         s.add_col_multiple(k, j_t, factor);
//                         v.add_col_multiple(k, j_t, factor);
//                     }
//                 }

//                 // println!("Column Step III");
//                 // print_matrix(&s);

//                 if !improved {
//                     break;
//                 }
//             }
            
            
//             // // Clear row t to the right of position (t, j_t) using column operations
//             // loop {
//             //     let mut cleared_row = true;
//             //     for j in (j_t + 1)..n {
//             //         let entry = s.get(j, t);
//             //         if !entry.is_zero() {
//             //             cleared_row = false;
                        
//             //             // Simple column elimination: subtract multiple of pivot column
//             //             let pivot = s.get(j_t, t);
//             //             if !pivot.is_zero() {
//             //                 // Use simple approach: subtract pivot column from current column
//             //                 s.add_col_multiple(j, j_t, -R::one());
//             //                 v.add_col_multiple(j, j_t, -R::one());
//             //             }
//             //             break;
//             //         }
//             //     }
                
//             //     if cleared_row {
//             //         break;
//             //     }
                
//             //     // Re-clear column if needed (entries may have become non-zero again)
//             //     let pivot = s.get(j_t, t);
//             //     for k in (t + 1)..m {
//             //         let entry = s.get(j_t, k);
//             //         if !entry.is_zero() && !pivot.is_zero() {
//             //             // Simple row elimination
//             //             s.add_row_multiple(k, t, -R::one());
//             //             u.add_row_multiple(k, t, -R::one());
//             //         }
//             //     }
//             // }
//         }
        
//         // Final step: ensure diagonal entries satisfy divisibility condition
//         // self.ensure_divisibility(&mut s, &mut u, &mut v);
        
//         (u, s, v)
//     }
// }
