use crate::linalg::field::Field;

use std::cmp::Ordering;
use std::fmt::Debug;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

#[derive(Debug, Clone)]
pub struct Monomial<const N: usize>(pub [u16; N]);

impl<const N: usize> PartialEq for Monomial<N> {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl<const N: usize> PartialOrd for Monomial<N> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.0.partial_cmp(&other.0)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct MultivariatePolynomial<F: Field, const N: usize>(pub Vec<(Monomial<N>, F)>);

impl<F: Field, const N: usize> MultivariatePolynomial<F, N> {
    pub fn default() {
        Self::zero();
    }

    pub fn new(init: Vec<(Monomial<N>, F)>) -> Self {
        let mut n = Self(init);
        n.minimize();
        n
    }

    fn minimize(&mut self) {
        self.0
            .sort_by(|(m1, _), (m2, _)| m1.partial_cmp(m2).unwrap());
        let mut minimized: Vec<(Monomial<N>, F)> = Vec::new();
        for (monomial, coeff) in self.0.drain(..) {
            if let Some((last_m, last_c)) = minimized.pop() {
                if last_m == monomial {
                    let new_c = last_c + coeff;
                    minimized.push((monomial, new_c));
                    continue;
                }
            }
            minimized.push((monomial, coeff));
        }
        self.0 = minimized;
    }

    pub fn is_minimal(&self) -> bool {
        for i in 1..self.0.len() {
            if self.0[i - 1].0 > self.0[i].0 {
                return false;
            }
            if self.0[i - 1].0 == self.0[i].0 {
                return false;
            }
            if self.0[i].1 == F::zero() {
                return false;
            }
        }
        true
    }
}

impl<F: Field, const N: usize> Add for MultivariatePolynomial<F, N> {
    type Output = Self;
    fn add(mut self, rhs: Self) -> Self::Output {
        self.0.extend(rhs.0);
        self.minimize();
        self
    }
}

impl<F: Field, const N: usize> Sub for MultivariatePolynomial<F, N> {
    type Output = Self;

    fn sub(mut self, rhs: Self) -> Self::Output {
        self.0.extend(rhs.0.into_iter().map(|(m, c)| (m, -c)));
        self.minimize();
        self
    }
}

impl<F: Field, const N: usize> Neg for MultivariatePolynomial<F, N> {
    type Output = Self;
    fn neg(mut self) -> Self::Output {
        self.0.iter_mut().for_each(|(_, c)| *c = -c.clone());
        self
    }
}

impl<F: Field, const N: usize> Mul for MultivariatePolynomial<F, N> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        let mut result = Vec::new();
        for (m1, c1) in &self.0 {
            for (m2, c2) in &rhs.0 {
                let mut new_monomial = [0; N];
                for i in 0..N {
                    new_monomial[i] = m1.0[i] + m2.0[i];
                }
                result.push((Monomial(new_monomial), c1.clone() * c2.clone()));
            }
        }
        let mut product = MultivariatePolynomial(result);
        product.minimize();
        product
    }
}

impl<F: Field, const N: usize> AddAssign for MultivariatePolynomial<F, N> {
    fn add_assign(&mut self, rhs: Self) {
        self.0.extend(rhs.0);
        self.minimize();
    }
}

impl<F: Field, const N: usize> SubAssign for MultivariatePolynomial<F, N> {
    fn sub_assign(&mut self, rhs: Self) {
        self.0.extend(rhs.0.into_iter().map(|(m, c)| (m, -c)));
        self.minimize();
    }
}

impl<F: Field, const N: usize> MulAssign for MultivariatePolynomial<F, N> {
    fn mul_assign(&mut self, rhs: Self) {
        *self = self.clone() * rhs;
    }
}

impl<F: Field, const N: usize> MultivariatePolynomial<F, N> {
    pub fn eval(&self, _values: &[F; N]) -> F {
        todo!();
        // let mut result = F::zero();
        // for (monomial, coeff) in &self.0 {
        //     let term_value = monomial.0.iter().zip(values.iter()).map(|(&exp, &val)| val.clone().pow(exp as u32)).product();
        //     result = result + coeff.clone() * term_value;
        // }
        // result
    }

    pub fn scalar_mult(&mut self, scalar: F) {
        for (_, coeff) in &mut self.0 {
            *coeff = coeff.clone() * scalar.clone();
        }
        self.minimize();
    }

    pub fn is_zero(&self) -> bool {
        self.0.is_empty()
    }

    pub fn one() -> Self {
        MultivariatePolynomial(vec![(Monomial([0; N]), F::one())])
    }

    pub fn zero() -> Self {
        MultivariatePolynomial(vec![])
    }

    pub fn parse(_input: &str) -> Result<Self, String> {
        Err("Parsing not implemented yet".to_string())
    }
}
