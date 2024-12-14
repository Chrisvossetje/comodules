#[derive(Debug, Clone, PartialEq)]
pub struct Page {
    pub name: String,
    pub id: usize,
    pub degrees: Vec<char>,
    pub x_formula: String,
    pub y_formula: String,
    pub generators: Vec<(usize, usize, Vec<i32>, Option<String>)>,
    pub structure_lines: Vec<((usize, usize), (usize, usize), usize, String)>,
    pub differentials: Vec<((usize, usize), (usize, usize))>,
}
