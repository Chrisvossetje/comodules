use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

pub trait PrimeField : Sized + Eq + Copy + Add<Output = Self> + Sub<Output = Self> + Neg<Output = Self> + Mul<Output = Self> + AddAssign + SubAssign + MulAssign + From<u64>  {
    fn get_characteristic(&self) -> usize;
    fn get_degree(&self) -> usize;
    fn get_order(&self) -> usize { self.get_characteristic().pow(self.get_degree() as u32) }

    
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct F2 {
    value: u8
}

impl Add for F2 {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        F2 { value: self.value ^ rhs.value }
    }
}

impl Sub for F2 {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        F2 { value: self.value ^ rhs.value }
    }
    
}

impl Neg for F2 {
    type Output = Self;

    fn neg(self) -> Self {
        self
    }
    
}

impl Mul for F2 {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        F2 { value: self.value & rhs.value }
    }
}

impl AddAssign for F2 {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl SubAssign for F2 {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl MulAssign for F2 {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl From<u64> for F2 {
    fn from(value: u64) -> Self {
        F2 { value: (value % 2) as u8 }
    }
}

impl PrimeField for F2 {
    fn get_characteristic(&self) -> usize {
        2
    }

    fn get_degree(&self) -> usize {
        1
    }
}