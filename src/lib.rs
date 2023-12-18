//extern crate f128;
//extern crate itertools;
extern crate pest;
#[macro_use]
extern crate pest_derive;
//
//use num_traits::{Float, FloatConst, FromPrimitive, Num, Signed, ToPrimitive, Zero};
//
//use itertools::Itertools;
//
//use num::traits::{NumAssign, NumOps, NumRef};
//use num::NumCast;
//use num::{BigInt, BigRational};
//use std::fmt::{Debug, Display, LowerExp};
//use std::iter::Sum;
//use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
////use utils::Signum;
////pub mod test;
//
// TODO: use it
pub const MAX_VARIABLE: usize = 10;

pub mod qcoeff;
pub mod qmpoly;
pub mod qmrat;
pub mod parser;
//use utils::{multinomial, next_combination_with_replacement};
