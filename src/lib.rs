//extern crate f128;
//extern crate itertools;
extern crate pest;
#[macro_use]
extern crate pest_derive;

//use itertools::Itertools;
// TODO: use it
pub const MAX_VARIABLE: usize = 10;

pub mod qcoeff;
pub mod qmpoly;
pub mod qmrat;
pub mod parser;
//use utils::{multinomial, next_combination_with_replacement};
