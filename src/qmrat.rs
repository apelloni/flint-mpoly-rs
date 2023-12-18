use regex::Regex;
use crate::qcoeff::QCoeff;
use crate::qmpoly::QMPoly;
use flint_sys::fmpq_mpoly::*;
use std::fmt;
use std::ops::{Add, Div, Mul, Sub};

/// Container for rational polynomials
pub struct QMRat {
    pub num: QMPoly,       // numerator
    pub den: QMPoly,       // Polynomial context
    pub vars: Vec<String>, // Name of variables
    _reduced: bool,
    _coeff: QCoeff, // coefficient container
    _gcd: QMPoly,   // gcd container
    _res: QMPoly,   // partial operations container
}

impl QMRat {
    /// Create a rational function equal to 0
    pub fn new(vars: &[String]) -> Self {
        let num = QMPoly::new(vars);
        let mut den = QMPoly::new(vars);
        den.set_to_one();
        QMRat {
            num,
            den,
            vars: vars.to_vec(),
            _reduced: false,
            _coeff: QCoeff::default(),
            _gcd: QMPoly::new(vars),
            _res: QMPoly::new(vars),
        }
    }
    /// Create a rational function from a polynomial
    pub fn from_mpoly(mpoly: &QMPoly) -> Self {
        let mut den = QMPoly::new(&mpoly.vars);
        den.set_to_one();
        QMRat {
            num: QMPoly::clone_from(mpoly),
            den,
            vars: mpoly.vars.clone(),
            _reduced: false,
            _coeff: QCoeff::default(),
            _gcd: QMPoly::new(&mpoly.vars),
            _res: QMPoly::new(&mpoly.vars),
        }
    }
    /// Create a rational function from two polynomials
    pub fn from_mpolys(num: &QMPoly, den: &QMPoly) -> Self {
        QMRat {
            num: QMPoly::clone_from(num),
            den: QMPoly::clone_from(den),
            vars: num.vars.clone(),
            _reduced: false,
            _coeff: QCoeff::default(),
            _gcd: QMPoly::new(&num.vars),
            _res: QMPoly::new(&num.vars),
        }
    }
    /// Clone will copy the pointers and we will have two objects pointing at the same memory
    /// address. With this function we can create a clone in a safe way
    pub fn clone_from(other: &Self) -> Self {
        let mut this = QMRat::new(&other.vars);
        unsafe {
            // copy numerator
            fmpq_mpoly_set(
                &mut this.num.raw as *mut _,
                &other.num.raw as *const _ as *mut _,
                &other.num.ctx as *const _ as *mut _,
            );
            // copy denominator
            fmpq_mpoly_set(
                &mut this.den.raw as *mut _,
                &other.den.raw as *const _ as *mut _,
                &other.den.ctx as *const _ as *mut _,
            );
        }
        this
    }
    /// Chec if the function is zero
    pub fn is_zero(&self) -> bool{
        self.num.is_zero()
    }
    /// Check if the function is const
    pub fn is_const(&self) -> bool {
        self.num.is_const() && self.den.is_const()
    }
    /// Compute gcd of num end den
    fn gcd(&self) {
        unsafe {
            fmpq_mpoly_gcd(
                &self._gcd.raw as *const _ as *mut _,
                &self.num.raw as *const _ as *mut _,
                &self.den.raw as *const _ as *mut _,
                &self.num.ctx as *const _ as *mut _,
            );
        }
    }
    /// Compute gcd of num end den
    pub fn reduce(&mut self) {
        self.gcd();
        let mut exact = true;
        unsafe {
            // Simplify the numerator
            exact &= fmpq_mpoly_divides(
                &self._res.raw as *const _ as *mut _,
                &self.num.raw as *const _ as *mut _,
                &self._gcd.raw as *const _ as *mut _,
                &self.num.ctx as *const _ as *mut _,
            ) == 1;
            fmpq_mpoly_swap(
                &self._res.raw as *const _ as *mut _,
                &self.num.raw as *const _ as *mut _,
                &self.num.ctx as *const _ as *mut _,
            );
            // Simplify the denominator
            exact &= fmpq_mpoly_divides(
                &self._res.raw as *const _ as *mut _,
                &self.den.raw as *const _ as *mut _,
                &self._gcd.raw as *const _ as *mut _,
                &self.den.ctx as *const _ as *mut _,
            ) == 1;
            fmpq_mpoly_swap(
                &self._res.raw as *const _ as *mut _,
                &self.den.raw as *const _ as *mut _,
                &self.den.ctx as *const _ as *mut _,
            );
        }
        assert!(exact);
        self._reduced = true;
    }
    /// Normalize to have a monic polynomial at the denominator
    pub fn normalize(&mut self) {
        unsafe {
            let length =
                fmpq_mpoly_length(&mut self.den.raw as *mut _, &mut self.den.ctx as *mut _);
            // Get leading coefficient
            if length > 0 {
                fmpq_mpoly_get_term_coeff_fmpq(
                    &mut self._coeff.raw as *mut _,
                    &mut self.den.raw as *mut _,
                    0,
                    &mut self.den.ctx as *const _ as *mut _,
                );
                if !self._coeff.is_zero() && !self._coeff.is_one() {
                    // scale numerator
                    fmpq_mpoly_scalar_div_fmpq(
                        &mut self._res.raw as *mut _,
                        &mut self.num.raw as *mut _,
                        &mut self._coeff.raw as *mut _,
                        &mut self.num.ctx as *mut _,
                    );
                    fmpq_mpoly_swap(
                        &mut self._res.raw as *mut _,
                        &mut self.num.raw as *mut _,
                        &mut self.num.ctx as *mut _,
                    );
                    // scale denominator
                    fmpq_mpoly_scalar_div_fmpq(
                        &mut self._res.raw as *mut _,
                        &mut self.den.raw as *mut _,
                        &mut self._coeff.raw as *mut _,
                        &mut self.den.ctx as *mut _,
                    );
                    fmpq_mpoly_swap(
                        &mut self._res.raw as *mut _,
                        &mut self.den.raw as *mut _,
                        &mut self.den.ctx as *mut _,
                    );
                }
            }
        }
    }
    /// Fully reduce by finding the gcd and normalizing to the denominator leading coefficient
    pub fn full_reduce(&mut self) {
        self.reduce();
        self.normalize();
    }
    /// Rise function to integer power
    pub fn pown(&mut self, n: u64) {
        // Bring franction to the simplest form before applying the power
        if !self._reduced {
            self.full_reduce();
        }
        unsafe {
            // Rise numerator
            fmpq_mpoly_pow_ui(
                &mut self._res.raw as *mut _,
                &mut self.num.raw as *mut _,
                n,
                &mut self.num.ctx as *mut _,
            );
            fmpq_mpoly_swap(
                &mut self._res.raw as *mut _,
                &mut self.num.raw as *mut _,
                &mut self.num.ctx as *mut _,
            );
            // Rise denominator
            fmpq_mpoly_pow_ui(
                &mut self._res.raw as *mut _,
                &mut self.den.raw as *mut _,
                n,
                &mut self.den.ctx as *mut _,
            );
            fmpq_mpoly_swap(
                &mut self._res.raw as *mut _,
                &mut self.den.raw as *mut _,
                &mut self.den.ctx as *mut _,
            );
        }
    }
    /// Clear the function and set it to zero
    pub fn clear(&mut self) {
        self.num.set_to_zero();
        self.den.set_to_one();
        // Free the coefficient of the helper functions
        self._gcd.set_to_zero();
        self._res.set_to_zero();
    }
    /// Convert the function to a human readable string
    pub fn to_str(&self) -> String {
        let rg = Regex::new(r"([+-])1\*").unwrap();
        rg.replace_all(format!("{}", self).as_str(), "$1").to_string()
    }
}

///Implement addition with operator `+`
impl<'a, 'b> Add<&'b QMRat> for &'a QMRat {
    type Output = QMRat;

    fn add(self, other: &'b QMRat) -> QMRat {
        let mut this = QMRat::new(&self.vars);
        this.num = &(&self.num * &other.den) + &(&other.num * &self.den);
        this.den = &self.den * &other.den;
        this.full_reduce();
        this
    }
}

///Implement addition with operator `-`
impl<'a, 'b> Sub<&'b QMRat> for &'a QMRat {
    type Output = QMRat;

    fn sub(self, other: &'b QMRat) -> QMRat {
        let mut this = QMRat::new(&self.vars);
        this.num = &(&self.num * &other.den) - &(&other.num * &self.den);
        this.den = &self.den * &other.den;
        this.full_reduce();
        this
    }
}
///Implement addition with operator `*`
impl<'a, 'b> Mul<&'b QMRat> for &'a QMRat {
    type Output = QMRat;

    fn mul(self, other: &'b QMRat) -> QMRat {
        let mut this = QMRat::new(&self.vars);
        this.num = &self.num * &other.num;
        this.den = &self.den * &other.den;
        this.full_reduce();
        this
    }
}
///Implement addition with operator `/`
impl<'a, 'b> Div<&'b QMRat> for &'a QMRat {
    type Output = QMRat;

    fn div(self, other: &'b QMRat) -> QMRat {
        let mut this = QMRat::new(&self.vars);
        this.num = &self.num * &other.den;
        this.den = &self.den * &other.num;
        this.full_reduce();
        this
    }
}

// To use the `{}` for the structure QMRat
impl fmt::Display for QMRat {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "({})/({})", self.num, self.den)
    }
}
