use crate::{matrices::f2_matrix::F2Matrix, matrix::Matrix};


#[derive(Debug)]
pub(crate) struct M4riTableU64 {
    rows: Vec<usize>,     // pivot row indices (in current block)
    cols: Vec<usize>,     // pivot columns (same order as rows)
    data: Vec<u64>,       // flattened table (2^k - 1) * wpr words
    min_word: usize,      // smallest word index among pivot columns
    wpr: usize,
}

impl M4riTableU64 {
    fn new(k: usize, wpr: usize) -> Self {
        Self {
            rows: Vec::with_capacity(k),
            cols: Vec::with_capacity(k),
            data: Vec::with_capacity(((1usize << k) - 1) * wpr),
            min_word: usize::MAX,
            wpr,
        }
    }

    fn len(&self) -> usize { self.cols.len() }
    fn is_empty(&self) -> bool { self.cols.is_empty() }
    fn rows(&self) -> &[usize] { &self.rows }

    fn clear(&mut self) {
        self.rows.clear();
        self.cols.clear();
        self.data.clear();
        self.min_word = usize::MAX;
    }

    fn add(&mut self, pivot_col: usize, row: usize) {
        self.cols.push(pivot_col);
        self.rows.push(row);
    }

    fn xor_from_word(dst: &mut [u64], src: &[u64], start_word: usize) {
        for i in start_word..dst.len() {
            dst[i] ^= src[i];
        }
    }

    /// Build table of all nonzero XOR combinations of pivot rows in current block.
    /// Matches `m4ri.rs::generate` doubling pattern.
    fn generate(&mut self, matrix: &F2Matrix) {
        debug_assert_eq!(self.wpr, matrix.words_per_row);
        let wpr = self.wpr;

        for (n, (&c, &r)) in self.cols.iter().zip(self.rows.iter()).enumerate() {
            let row_words = matrix.get_row(r);
            let old_len = self.data.len();

            // append base row
            self.data.extend_from_slice(row_words);
            // duplicate previous combinations
            self.data.extend_from_within(0..old_len);

            // xor base row into the duplicated block
            let start = 1usize << n;
            let end = (1usize << (n + 1)) - 1;
            let start_word = c >> 6;

            for idx in start..end {
                let dst = &mut self.data[idx * wpr .. (idx + 1) * wpr];
                Self::xor_from_word(dst, row_words, start_word);
            }

            self.min_word = self.min_word.min(start_word);
        }
    }

    /// Reduce a row using the precomputed combination table.
    /// Matches `m4ri.rs::reduce`: index from pivot bits, then XOR one table row.
    fn reduce(&self, row: &mut [u64]) {
        let mut index: usize = 0;

        // IMPORTANT: match m4ri.rs: iterate cols in reverse
        for &c in self.cols.iter().rev() {
            let word = c >> 6;
            let bit  = c & 63;
            index <<= 1;
            index |= ((row[word] >> bit) & 1) as usize;
        }

        if index != 0 {
            let base = (index - 1) * self.wpr;
            let src = &self.data[base .. base + self.wpr];
            Self::xor_from_word(row, src, self.min_word);
        }
    }

    /// Reduce matrix[target] by current pivot rows *naively* (like `m4ri.rs::reduce_naive`)
    fn reduce_naive(&self, matrix: &mut F2Matrix, target: usize) {
        for (&r, &c) in self.rows.iter().zip(self.cols.iter()) {
            debug_assert_ne!(r, target);
            let word = c >> 6;
            let bit  = c & 63;
            let mask = 1u64 << bit;
            if (matrix.data[target  * self.wpr + word] & mask) != 0 {
                matrix.xor_row_from_word(target, r, word);
            }
        }
    }
}

impl F2Matrix {
    fn first_one_in_row(&self, row: usize) -> Option<usize> {
        for (col, word) in self.get_row(row).iter().enumerate() {
            if *word != 0 {
                let tz = word.trailing_zeros() as usize;
                return Some((col << 6) + tz);
            }
        }
        None
    }

    pub(crate) fn echelonize_m4ri(&mut self) {
        let rows = self.codomain;
        let cols = self.domain;
        let wpr  = self.words_per_row;

        // In matrix_inner.rs, pivots is column->pivot_row mapping.
        let mut empty_rows: Vec<usize> = Vec::with_capacity(rows);
        let mut col_to_pivot_row: Vec<isize> = vec![-1; cols];

        if rows == 0 || cols == 0 {
            self.pivots = vec![None; cols];
            return;
        }

        let k = 6; // This k seems to be the best
        let mut table = M4riTableU64::new(k, wpr);

        // === Main loop: exactly like matrix_inner.rs F2 path ===
        for i in 0..rows {
            // 1) reduce row i by current table pivots (naive)
            table.reduce_naive(self, i);

            // 2) find first nonzero in row i => pivot column
            if let Some(c) = self.first_one_in_row(i) {

                // 3) record pivot column -> pivot row
                col_to_pivot_row[c] = i as isize;

                // 4) clear pivot column from existing pivot rows in table (Gauss-Jordan step)
                let word = c >> 6;
                let bit  = c & 63;
                let mask = 1u64 << bit;

                // FOR REDUCED ??
                for &r in table.rows() {
                    if r == i { panic!("this row cannot be in the table yet"); }
                    if (self.get_row(r)[word] & mask) != 0 {
                        self.xor_row_from_word(r, i, word);
                    }
                }

                // 5) add pivot row to table
                table.add(c, i);

                // 6) once table is full: generate and reduce
                if table.len() == k {
                    table.generate(self);

                    // reduce rows above first pivot row in block
                    let first = table.rows()[0];
                    for j in 0..first {
                        let rowj = self.get_row_mut(j);
                        table.reduce(rowj);
                    }

                    // reduce rows below current i
                    for j in (i + 1)..rows {
                        let rowj = self.get_row_mut(j);
                        table.reduce(rowj);
                    }

                    table.clear();
                }
            } else {
                empty_rows.push(i);
            }
        }
        
        // === Flush leftover table: generate then reduce rows above first pivot row ===
        if !table.is_empty() {
            table.generate(self);
            let first = table.rows()[0];
            for j in 0..first {
                let rowj = self.get_row_mut(j);
                table.reduce(rowj);
            }
            table.clear();
        }

        self.pivots = vec![None; cols];
        for (col, &row) in col_to_pivot_row.iter().enumerate() {
            if row >= 0 {
                self.pivots[col] = Some(row as u32);
            }
        }
    }
}
