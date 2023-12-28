//! Polynomial manipulation using the powerful [Flint] (Fast Library for Number Theory)
//! library to deal with the arithmetic over a variaty of rings, such as,
//! - __Integers__:
//!     - [Coefficients](ZCoeff)
//!     - [Polynomials](ZMPoly)
//!     - [Rational Functions](ZMRat)
//! - __Rationals__:
//!     - [Coefficients](QCoeff)
//!     - [Polynomials](QMPoly)
//!     - [Rational Functions](QMRat)
//! - __Finite Fields__:
//!     - [Coefficients](ModCoeff)
//!     - ~~Polynomials~~
//!     - ~~Rational Functions~~
//!
//! [Flint]: https://flintlib.org
//!
extern crate pest;
#[macro_use]
extern crate pest_derive;

//use itertools::Itertools;
// TODO: use it
pub const MAX_VARIABLE: usize = 10;

mod parser;
/// Rational
mod qcoeff;
mod qmpoly;
mod qmrat;
/// Integer
mod zcoeff;
mod zmpoly;
mod zmrat;
/// Modulus
mod modcoeff;

pub use crate::parser::{parse_mpoly, parse_mrat};
/// Rational
pub use crate::qcoeff::QCoeff;
pub use crate::qmpoly::QMPoly;
pub use crate::qmrat::QMRat;
/// Integer
pub use crate::zcoeff::ZCoeff;
pub use crate::zmpoly::ZMPoly;
pub use crate::zmrat::ZMRat;
/// Modulus
pub use crate::modcoeff::ModCoeff;
