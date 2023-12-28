use crate::parse_mrat;
use crate::ZCoeff;
use crate::ZMPoly;
use flint_sys::fmpz_mpoly::*;
use regex::Regex;
use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};

/// Container for rational polynomials over rational numbers
#[derive(Debug)]
pub struct ZMRat {
    pub num: ZMPoly,       // numerator
    pub den: ZMPoly,       // Polynomial context
    pub vars: Vec<String>, // Name of variables
    _reduced: bool,
    _coeff: ZCoeff, // coefficient container
    _gcd: ZMPoly,   // gcd container
    _res: ZMPoly,   // partial operations container
}

impl ZMRat {
    /// Create a rational function equal to 0
    /// ```
    /// use flint_mpoly::ZMRat;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mrat = ZMRat::new(&vars);
    /// assert_eq!("0",mrat.to_str());
    /// ```
    pub fn new(vars: &[String]) -> Self {
        let num = ZMPoly::new(vars);
        let mut den = ZMPoly::new(vars);
        den.set_to_one();
        ZMRat {
            num,
            den,
            vars: vars.to_vec(),
            _reduced: false,
            _coeff: ZCoeff::default(),
            _gcd: ZMPoly::new(vars),
            _res: ZMPoly::new(vars),
        }
    }

    /// Create a rational function from a polynomial
    /// ```
    /// use flint_mpoly::ZMPoly;
    /// use flint_mpoly::ZMRat;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mpoly= ZMPoly::from_str("x1+3 x2",&vars).unwrap();
    /// let mrat = ZMRat::from_mpoly(&mpoly);
    /// assert_eq!("+x1+3*x2",mrat.to_str());
    /// ```
    pub fn from_mpoly(mpoly: &ZMPoly) -> Self {
        let mut den = ZMPoly::new(&mpoly.vars);
        den.set_to_one();
        ZMRat {
            num: ZMPoly::clone_from(mpoly),
            den,
            vars: mpoly.vars.clone(),
            _reduced: false,
            _coeff: ZCoeff::default(),
            _gcd: ZMPoly::new(&mpoly.vars),
            _res: ZMPoly::new(&mpoly.vars),
        }
    }

    /// Create a rational function from two polynomials defining the numerator and denominator,
    /// respectively
    /// ```
    /// use flint_mpoly::ZMPoly;
    /// use flint_mpoly::ZMRat;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let num= ZMPoly::from_str("x1+3 x2",&vars).unwrap();
    /// let den= ZMPoly::from_str("x1-x2",&vars).unwrap();
    /// let mrat = ZMRat::from_mpolys(&num,&den);
    /// assert_eq!("(+x1+3*x2)/(+x1-x2)",mrat.to_str());
    /// ```
    pub fn from_mpolys(num: &ZMPoly, den: &ZMPoly) -> Self {
        ZMRat {
            num: ZMPoly::clone_from(num),
            den: ZMPoly::clone_from(den),
            vars: num.vars.clone(),
            _reduced: false,
            _coeff: ZCoeff::default(),
            _gcd: ZMPoly::new(&num.vars),
            _res: ZMPoly::new(&num.vars),
        }
    }

    /// Initialize rational from string assuming the given variables
    /// ```
    /// use flint_mpoly::ZMRat;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mrat = ZMRat::from_str("(x1+x2)/x1", &vars).unwrap();
    /// assert_eq!("(+x1+x2)/(+x1)",mrat.to_str());
    /// ```
    pub fn from_str(expr_str: &str, vars: &[String]) -> Result<Self, &'static str> {
        parse_mrat(expr_str, vars).expect("parse rational").to_zmrat()
    }

    /// Set rational from string assuming the variables already stored in ZMRat
    /// ```
    /// use flint_mpoly::ZMRat;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mut mrat = ZMRat::new(&vars);
    /// mrat.set_from_str("(x1+x2)/x1").unwrap();
    /// assert_eq!("(+x1+x2)/(+x1)",mrat.to_str());
    /// ```
    pub fn set_from_str(&mut self, expr_str: &str) -> Result<(), &'static str> {
        *self = ZMRat::from_str(expr_str, &self.vars)?;
        Ok(())
    }

    /// Clone will copy the pointers and we will have two objects pointing at the same memory
    /// address. With this function we can create a clone in a safe way
    /// ```
    /// use flint_mpoly::ZMRat;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mrat_1 = ZMRat::from_str("x1/(x2+1)",&vars).unwrap();
    /// let mrat_2 = ZMRat::clone_from(&mrat_1);
    /// assert_eq!(mrat_1.to_str(),mrat_2.to_str());
    /// ```
    pub fn clone_from(other: &Self) -> Self {
        let mut this = ZMRat::new(&other.vars);
        unsafe {
            // copy numerator
            fmpz_mpoly_set(
                &mut this.num.raw as *mut _,
                &other.num.raw as *const _ as *mut _,
                &other.num.ctx as *const _ as *mut _,
            );
            // copy denominator
            fmpz_mpoly_set(
                &mut this.den.raw as *mut _,
                &other.den.raw as *const _ as *mut _,
                &other.den.ctx as *const _ as *mut _,
            );
        }
        this
    }
    /// Swap the structure content with another ZMRat provided they have the same context
    /// ```
    /// use flint_mpoly::ZMRat;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mut mrat_a = ZMRat::from_str("(x1+x2)/x1",&vars).unwrap();
    /// let mut mrat_b = ZMRat::from_str("x2/(x1-x2)",&vars).unwrap();
    /// assert_eq!("(+x1+x2)/(+x1)",mrat_a.to_str());
    /// assert_eq!("(+x2)/(+x1-x2)",mrat_b.to_str());
    /// mrat_a.swap(&mut mrat_b);
    /// assert_eq!("(+x1+x2)/(+x1)",mrat_b.to_str());
    /// assert_eq!("(+x2)/(+x1-x2)",mrat_a.to_str());
    /// ```
    pub fn swap(&mut self, other: &mut Self) {
        assert!(
            self.vars == other.vars,
            "Unmatch variables while swapping expressions"
        );
        unsafe {
            fmpz_mpoly_swap(
                &mut self.num.raw as *mut _,
                &mut other.num.raw as *mut _,
                &mut self.num.ctx as *mut _,
            );
            fmpz_mpoly_swap(
                &mut self.den.raw as *mut _,
                &mut other.den.raw as *mut _,
                &mut self.den.ctx as *mut _,
            );
        }
    }
    /// Set function to zero and return the previous value
    /// ```
    /// use flint_mpoly::ZMRat;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mut mrat_a = ZMRat::from_str("(x1+x2)/x1",&vars).unwrap();
    /// assert_eq!("(+x1+x2)/(+x1)",mrat_a.to_str());
    /// let mut mrat_b = mrat_a.zero_move();
    /// assert_eq!("0",mrat_a.to_str());
    /// assert_eq!("(+x1+x2)/(+x1)",mrat_b.to_str());
    /// ```
    pub fn zero_move(&mut self) -> Self {
        let mut out = Self::new(&self.vars);
        out.swap(self);
        out
    }
    /// Chec if the function is zero
    pub fn is_zero(&self) -> bool {
        self.num.is_zero()
    }
    /// Check if the function is const
    pub fn is_const(&self) -> bool {
        self.num.is_const() && self.den.is_const()
    }
    /// Compute gcd of num end den
    fn gcd(&self) {
        unsafe {
            fmpz_mpoly_gcd(
                &self._gcd.raw as *const _ as *mut _,
                &self.num.raw as *const _ as *mut _,
                &self.den.raw as *const _ as *mut _,
                &self.num.ctx as *const _ as *mut _,
            );
        }
    }
    /// Get GCD bewtween numerator and denominator
    /// ```
    /// use flint_mpoly::ZMRat;
    /// use flint_mpoly::parse_mpoly;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mut mrat = ZMRat::new(&vars);
    /// mrat.num.set_from_str("(1+x1)*(x2-x1)").unwrap();
    /// mrat.den.set_from_str("(x1+x2)*(x2-x1)").unwrap();
    /// let gcd = mrat.get_gcd();
    /// assert_eq!("+x1-x2", gcd.to_str());
    /// ```
    pub fn get_gcd(&self) -> ZMPoly {
        self.gcd();
        ZMPoly::clone_from(&self._gcd)
    }
    /// Reduce by dividing both numerator and denominators by the gcd.
    /// This expression is not guaranteed to return a denominaor or numerator in monic form.
    /// ```
    /// use flint_mpoly::ZMRat;
    /// use flint_mpoly::parse_mpoly;
    /// let vars = [String::from("x1"),String::from("x2")];
    /// let mut mrat = ZMRat::new(&vars);
    /// mrat.num.set_from_str("(1+x1)*(x2-x1)").unwrap();
    /// mrat.den.set_from_str("(x1+x2)*(x2-x1)").unwrap();
    /// mrat.reduce();
    /// assert_eq!("(+x1+1)/(+x1+x2)", mrat.to_str());
    /// ```
    pub fn reduce(&mut self) {
        self.gcd();
        if self._gcd.is_one() {
            return;
        }
        let mut exact = true;
        unsafe {
            // Simplify the numerator
            exact &= fmpz_mpoly_divides(
                &self._res.raw as *const _ as *mut _,
                &self.num.raw as *const _ as *mut _,
                &self._gcd.raw as *const _ as *mut _,
                &self.num.ctx as *const _ as *mut _,
            ) == 1;
            fmpz_mpoly_swap(
                &self._res.raw as *const _ as *mut _,
                &self.num.raw as *const _ as *mut _,
                &self.num.ctx as *const _ as *mut _,
            );
            // Simplify the denominator
            exact &= fmpz_mpoly_divides(
                &self._res.raw as *const _ as *mut _,
                &self.den.raw as *const _ as *mut _,
                &self._gcd.raw as *const _ as *mut _,
                &self.den.ctx as *const _ as *mut _,
            ) == 1;
            fmpz_mpoly_swap(
                &self._res.raw as *const _ as *mut _,
                &self.den.raw as *const _ as *mut _,
                &self.den.ctx as *const _ as *mut _,
            );
        }
        assert!(exact, "Failed to reduce!");
        self._reduced = true;
    }

    /// Rise function to integer power
    pub fn pown(&mut self, n: u64) {
        // Bring franction to the simplest form before applying the power
        if !self._reduced {
            self.reduce();
        }
        unsafe {
            // Rise numerator
            fmpz_mpoly_pow_ui(
                &mut self._res.raw as *mut _,
                &mut self.num.raw as *mut _,
                n,
                &mut self.num.ctx as *mut _,
            );
            fmpz_mpoly_swap(
                &mut self._res.raw as *mut _,
                &mut self.num.raw as *mut _,
                &mut self.num.ctx as *mut _,
            );
            // Rise denominator
            fmpz_mpoly_pow_ui(
                &mut self._res.raw as *mut _,
                &mut self.den.raw as *mut _,
                n,
                &mut self.den.ctx as *mut _,
            );
            fmpz_mpoly_swap(
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
    /// Clear all memory allocated to ZMRat
    pub fn free(&mut self) {
        self.num.free();
        self.den.free();
        // Free the coefficient of the helper functions
        self._gcd.free();
        self._res.free();
    }
    /// Convert the function to a human readable string
    pub fn to_str(&self) -> String {
        let rg = Regex::new(r"([+-])1\*").unwrap();
        rg.replace_all(format!("{}", self).as_str(), "$1")
            .to_string()
    }
}

///Implement addition with operator `+`
/// ```
/// use flint_mpoly::ZMRat;
/// let vars = [String::from("x1"), String::from("x2")];
/// let mrat_a = ZMRat::from_str("1/x1",&vars).unwrap();
/// let mrat_b = ZMRat::from_str("1/x2",&vars).unwrap();
/// let mrat_ab = &mrat_a + &mrat_b;
/// assert_eq!("(+x1+x2)/(+x1*x2)", mrat_ab.to_str());
/// ```
impl<'a, 'b> Add<&'b ZMRat> for &'a ZMRat {
    type Output = ZMRat;

    fn add(self, other: &'b ZMRat) -> ZMRat {
        let mut this = ZMRat::new(&self.vars);
        this.num = &(&self.num * &other.den) + &(&other.num * &self.den);
        this.den = &self.den * &other.den;
        this.reduce();
        this
    }
}

///Implement addition with operator `-`
/// ```
/// use flint_mpoly::ZMRat;
/// let vars = [String::from("x1"), String::from("x2")];
/// let mrat_a = ZMRat::from_str("1/x1",&vars).unwrap();
/// let mrat_b = ZMRat::from_str("1/x2",&vars).unwrap();
/// let mrat_ab = &mrat_a - &mrat_b;
/// assert_eq!("(-x1+x2)/(+x1*x2)", mrat_ab.to_str());
/// ```
impl<'a, 'b> Sub<&'b ZMRat> for &'a ZMRat {
    type Output = ZMRat;

    fn sub(self, other: &'b ZMRat) -> ZMRat {
        let mut this = ZMRat::new(&self.vars);
        this.num = &(&self.num * &other.den) - &(&other.num * &self.den);
        this.den = &self.den * &other.den;
        this.reduce();
        this
    }
}
///Implement addition with operator `*`
/// ```
/// use flint_mpoly::ZMRat;
/// let vars = [String::from("x1"), String::from("x2")];
/// let mrat_a = ZMRat::from_str("1/x1",&vars).unwrap();
/// let mrat_b = ZMRat::from_str("1/x2",&vars).unwrap();
/// let mrat_ab = &mrat_a * &mrat_b;
/// assert_eq!("(+1)/(+x1*x2)", mrat_ab.to_str());
/// ```
impl<'a, 'b> Mul<&'b ZMRat> for &'a ZMRat {
    type Output = ZMRat;

    fn mul(self, other: &'b ZMRat) -> ZMRat {
        let mut this = ZMRat::new(&self.vars);
        this.num = &self.num * &other.num;
        this.den = &self.den * &other.den;
        this.reduce();
        this
    }
}
///Implement addition with operator `/`
/// ```
/// use flint_mpoly::ZMRat;
/// let vars = [String::from("x1"), String::from("x2")];
/// let mrat_a = ZMRat::from_str("1/x1",&vars).unwrap();
/// let mrat_b = ZMRat::from_str("1/x2",&vars).unwrap();
/// let mrat_ab = &mrat_a / &mrat_b;
/// assert_eq!("(+x2)/(+x1)", mrat_ab.to_str());
/// ```
impl<'a, 'b> Div<&'b ZMRat> for &'a ZMRat {
    type Output = ZMRat;

    fn div(self, other: &'b ZMRat) -> ZMRat {
        let mut this = ZMRat::new(&self.vars);
        this.num = &self.num * &other.den;
        this.den = &self.den * &other.num;
        this.reduce();
        this
    }
}

///Implement negative sign
/// ```
/// use flint_mpoly::ZMRat;
/// // Define new polynomial
/// let vars = [String::from("x1"), String::from("x2")];
/// let mrat = ZMRat::from_str("x1/x2",&vars).unwrap();
/// let mrat_neg = -&mrat;
/// assert_eq!("(-x1)/(+x2)", mrat_neg.to_str());
/// ```
impl<'a> Neg for &'a ZMRat {
    type Output = ZMRat;
    fn neg(self) -> ZMRat {
        // TODO: Check if the two contexts are the same
        let mut this = ZMRat::new(&self.vars);
        this.den = self.den.clone();
        unsafe {
            fmpz_mpoly_neg(
                &mut this.num.raw as *mut _,
                &self.num.raw as *const _ as *mut _,
                &self.num.ctx as *const _ as *mut _,
            );
        }
        this
    }
}

/// Clone for ZMRat
/// ```
/// use flint_mpoly::ZMRat;
/// use std::str::FromStr;
/// let vars = [String::from("x1"),String::from("x2"),String::from("x3")];
/// let mrat_1 = ZMRat::from_str("x1+x3/x2",&vars).unwrap();
/// let mrat_2 = mrat_1.clone();
/// assert_eq!(mrat_1.to_str(),mrat_2.to_str());
/// ```
impl Clone for ZMRat {
    fn clone(&self) -> Self {
        ZMRat::clone_from(self)
    }
}

// To use the `{}` for the structure ZMRat
impl fmt::Display for ZMRat {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.is_zero() {
            write!(f, "0")
        } else if self.den.is_one() {
            write!(f, "{}", self.num)
        } else {
            write!(f, "({})/({})", self.num, self.den)
        }
    }
}

// Send and Sync
unsafe impl Send for ZMRat {}
unsafe impl Sync for ZMRat {}
