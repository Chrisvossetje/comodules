use deepsize::DeepSizeOf;
use serde::{Deserialize, Serialize};
use std::fmt::Debug;

pub trait BasisElement: 'static + Debug + Clone + Default + DeepSizeOf {}

#[derive(Debug, Clone, PartialEq, Default, Deserialize, Serialize, DeepSizeOf)]
#[allow(non_camel_case_types)]
pub struct kBasisElement {
    pub name: String,
    pub generator: bool,
    pub primitive: Option<usize>,
    pub generated_index: usize,
}

impl BasisElement for kBasisElement {}
