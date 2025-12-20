// TODO : Fix this ?

// use std::{
//     fmt::{Debug, Display}, hash::Hash, iter::Sum, ops::{Add, AddAssign, Sub, SubAssign}, str::FromStr
// };

// use ahash::HashMap;

// use crate::{
//     basiselement::kBasisElement,
//     comodule::kcoalgebra::kCoalgebra,
//     grading::{Grading, OrderedGrading, Parse},
//     tensor::Tensor
// };

// use algebra::{self, field::Field, matrices::flat_matrix::FlatMatrix};

// pub trait Group:
//     'static
//     + Clone
//     + Hash
//     + Copy
//     + Debug
//     + Display
//     + Sized
//     + PartialEq
//     + Eq
//     + Parse
//     + FromStr
// {
//     fn mult(self, rhs: Self) -> Self;
//     fn inv(self) -> Self;

//     fn neutral() -> Self;

//     fn elements() -> Vec<Self>;

//     fn names() -> Vec<(usize, String)>;
//     fn generate_coalgebra<F: Field>() -> Result<kCoalgebra<Self, F, FlatMatrix<F>>, String>
//     where
//         Self: Grading;

//     fn toi32(self) -> i32;
// }

// macro_rules! define_group {
//     ($name:ident, $table:expr, $names:expr) => {
//         #[derive(Debug, Clone, Copy, Hash, PartialEq, Eq, PartialOrd, Ord, Default)]
//         pub struct $name(pub usize);

//         impl Display for $name {
//             fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
//                 write!(f, "{}", $names[self.0])
//             }
//         }

//         impl Parse for $name {
//             fn parse(s: &str) -> Result<Self, String> {
//                 for (i, name) in $names.iter().enumerate() {
//                     if s == *name {
//                         return Ok($name(i));
//                     }
//                 }
//                 Err(format!("Unknown element: {}", s))
//             }
//         }

//         impl FromStr for $name {
//             type Err = String;

//             fn from_str(s: &str) -> Result<Self, Self::Err> {
//                 Self::parse(s)
//             }
//         }

//         impl Group for $name {
//             fn mult(self, rhs: Self) -> Self {
//                 $name($table[self.0][rhs.0])
//             }

//             fn inv(self) -> Self {
//                 for i in 0..$names.len() {
//                     if $table[self.0][i] == 0 {
//                         return $name(i);
//                     }
//                 }
//                 panic!("No inverse found for element {}", self.0);
//             }

//             fn neutral() -> Self {
//                 $name(0)
//             }

//             fn elements() -> Vec<Self> {
//                 (0..$names.len()).map(|i| $name(i)).collect()
//             }

//             fn generate_coalgebra<F: Field>() -> Result<kCoalgebra<Self, F, FlatMatrix<F>>, String>
//             where Self: Grading
//             {

//                 let mut pre_space = HashMap::default();

//                 for (index, name) in $name::names() {
//                     pre_space.insert($name(index), vec![kBasisElement {
//                             name: name.clone(),
//                             generator: (index == 0),
//                             primitive: None,
//                             generated_index: 0,
//                     }]);
//                 }
//                 let space = GradedVectorSpace::<Self, kBasisElement>(pre_space);

//                 let tensor = Tensor::generate(&space, &space);

//                 let mut pre_coaction = HashMap::default();
//                 for g in Self::elements() {
//                     // All els (h,i) such that hi = g
//                     let mut els_to_g: Vec<(Self, Self)> = Vec::new();
//                     for h in Self::elements() {
//                         for i in Self::elements() {
//                             if Self::mult(h, i) == g {
//                                 els_to_g.push((h,i));
//                             }
//                         }
//                     }

//                     debug_assert!(els_to_g.len() >= 2, "At least the identity should go to g twice. (In group multiplication)");
//                     let tensor_dim = *tensor.dimensions.get(&g).ok_or("Tensor has no elements in this dimension")?;
//                     debug_assert_eq!(tensor_dim, els_to_g.len());
//                     let mut m = FlatMatrix::zero(1, tensor_dim);
//                     for c in 0..tensor_dim {
//                         let ((h, _), (i, _)) = tensor.deconstruct.get(&(g,c)).ok_or("Expected stuff in this element of g")?;
//                         if els_to_g.contains(&(*h,*i)) {
//                             m.set(0, c, F::one());
//                         }
//                     }
//                     pre_coaction.insert(g, m);
//                 }

//                 let coaction = GradedLinearMap::from(pre_coaction);

//                 let mut coalg = kCoalgebra {
//                     space,
//                     coaction,
//                     tensor,
//                 };
//                 coalg.reduce();
//                 Ok(coalg)
//             }

//             fn names() -> Vec<(usize, String)> {
//                 $names.iter().enumerate().map(|(index,s)| (index, s.to_string())).collect()
//             }

//             fn toi32(self) -> i32 {
//                 self.0 as i32
//             }
//         }

//         impl_group_arithmetic!($name);
//         impl_group_grading!($name);
//     };
// }

// macro_rules! impl_group_arithmetic {
//     ($type:ty) => {
//         impl Add for $type {
//             type Output = Self;

//             fn add(self, other: Self) -> Self::Output {
//                 self.mult(other)
//             }
//         }

//         impl AddAssign for $type {
//             fn add_assign(&mut self, other: Self) {
//                 *self = self.mult(other);
//             }
//         }

//         impl Sub for $type {
//             type Output = Self;

//             fn sub(self, other: Self) -> Self::Output {
//                 self.mult(other.inv())
//             }
//         }

//         impl SubAssign for $type {
//             fn sub_assign(&mut self, other: Self) {
//                 *self = self.mult(other.inv());
//             }
//         }

//         impl Sum for $type {
//             fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
//                 iter.fold(Self::neutral(), |a, b| a.mult(b))
//             }
//         }
//     };
// }

// macro_rules! impl_group_grading {
//     ($type:ty) => {
//         impl Grading for $type {
//             fn degree_names() -> Vec<char> {
//                 vec!['t']
//             }

//             fn default_formulas() -> (String, String) {
//                 ("t".to_string(), "s".to_string())
//             }

//             fn export_grade(self) -> Vec<i32> {
//                 vec![self.toi32()]
//             }

//             fn zero() -> Self {
//                 Self::neutral()
//             }

//             fn integer_multiplication(self, other: i32) -> Self {
//                 let mut e = Self::neutral();
//                 for _ in 0..other {
//                     e = Self::mult(self, e);
//                 }
//                 e
//             }
//             fn infty() -> Self {
//                 Self(usize::MAX - 100)
//             }
//         }
//     };
// }

// define_group!(Z2,
//     [
//         [0, 1],  // 0 * 0 = 0, 0 * 1 = 1
//         [1, 0],  // 1 * 0 = 1, 1 * 1 = 0
//     ],
//     ["0", "1"]
// );

// impl OrderedGrading for Z2 {
//     fn incr(self) -> Self {
//         Z2(self.0 + 1)
//     }

//     fn compare(self, rhs: &Self) -> std::cmp::Ordering {
//         self.0.cmp(&rhs.0)
//     }
// }

// define_group!(Z3,
//     [
//         [0, 1, 2],  // 0 * 0 = 0, 0 * 1 = 1
//         [1, 2, 0],  // 1 * 0 = 1, 1 * 1 = 0
//         [2, 0, 1],  // 1 * 0 = 1, 1 * 1 = 0
//     ],
//     ["0", "1", "2"]
// );
