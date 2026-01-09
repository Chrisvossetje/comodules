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
