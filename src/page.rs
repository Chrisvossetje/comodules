use std::{
    fs::File,
    io::{self, Write},
};

use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
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

impl Page {
    pub fn save_to_json(&self, file_path: String) -> io::Result<()> {
        let content = serde_json::to_string(self)?;

        // Create or open the file at the specified path
        let mut file = File::create(file_path)?;

        // Write the string content to the file
        file.write_all(content.as_bytes())
    }
}
