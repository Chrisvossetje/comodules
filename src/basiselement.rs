use std::fmt::Debug;
use serde::{Deserialize, Serialize};

pub trait BasisElement: 'static + Debug + Clone + Default {}

#[derive(Debug, Clone, PartialEq, Default, Deserialize, Serialize)]
#[allow(non_camel_case_types)]
pub struct kBasisElement {
    pub name: String,
    pub generator: bool,
    pub primitive: Option<usize>,
    pub generated_index: usize,
}

impl BasisElement for kBasisElement {}