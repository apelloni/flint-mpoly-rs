//extern crate f128;
//extern crate itertools;
extern crate pest;
#[macro_use]
extern crate pest_derive;

//use itertools::Itertools;
// TODO: use it
pub const MAX_VARIABLE: usize = 10;

mod parser;
mod qcoeff;
mod qmpoly;
mod qmrat;

pub use crate::parser::{parse_mpoly, parse_mrat};
pub use crate::qcoeff::QCoeff;
pub use crate::qmpoly::QMPoly;
pub use crate::qmrat::QMRat;
