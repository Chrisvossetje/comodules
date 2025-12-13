use std::{
    fs::File,
    io::{self, Write},
};

use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Page {
    pub id: usize,
    pub generators: Vec<(usize, usize, Vec<i64>, Option<String>)>,
    pub structure_lines: Vec<((usize, usize), (usize, usize), String, String)>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SSeq {
    pub name: String,
    pub degrees: Vec<char>,
    pub x_formula: String,
    pub y_formula: String,
    pub pages: Vec<Page>,
    pub differentials: Vec<((usize, usize), (usize, usize), usize)>,
}

impl SSeq {
    pub fn to_string(&self) -> String {
        serde_json::to_string(self).unwrap()
    }

    pub fn save_to_json(&self, file_path: &str) -> io::Result<()> {
        let content = serde_json::to_string(self)?;

        // Create or open the file at the specified path
        let mut file = File::create(file_path)?;

        // Write the string content to the file
        file.write_all(content.as_bytes())
    }
}
