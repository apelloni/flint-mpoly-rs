//extern crate f128;
//extern crate itertools;
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

pub use crate::parser::{parse_mpoly, parse_mrat};
/// Rational
pub use crate::qcoeff::QCoeff;
pub use crate::qmpoly::QMPoly;
pub use crate::qmrat::QMRat;
/// Integer
pub use crate::zcoeff::ZCoeff;
pub use crate::zmpoly::ZMPoly;
pub use crate::zmrat::ZMRat;
