use crate::ring::CRing;

pub trait Field: CRing {
    fn inv(self) -> Option<Self>;
    fn get_characteristic() -> usize;
}