use flint_sys::fmpq::*;
use flint_sys::fmpq_mpoly::*;
use flint_sys::mpoly::ordering_t_ORD_DEGLEX;
use std::borrow::Borrow;
use std::fmt;
use std::mem::MaybeUninit;
use std::ops::{Add, Div, Mul, Sub};

struct QCoeff {
    raw: fmpq, // coefficient container
}

impl QCoeff {
    // Initialize coefficient to 0.
    fn new() -> Self {
        unsafe {
            let mut c = MaybeUninit::uninit();
            fmpq_init(c.as_mut_ptr());
            fmpq_zero(c.as_mut_ptr());
            QCoeff {
                raw: c.assume_init(),
            }
        }
    }
    // Initialize coefficient to the canonical form of the fraction p / q.
    fn from_int(p: i64, q: u64) -> Self {
        unsafe {
            let mut c = MaybeUninit::uninit();
            fmpq_init(c.as_mut_ptr());
            fmpq_set_si(c.as_mut_ptr(), p, q);
            QCoeff {
                raw: c.assume_init(),
            }
        }
    }
    //Set a coefficient to the canonical form of the fraction p / q.
    fn set_from_int(&mut self, p: i64, q: u64) {
        unsafe {
            fmpq_set_si(&self.raw as *const _ as *mut _, p, q);
        }
    }
    /// Check if the coefficient is zero
    fn is_zero(&self) -> bool {
        unsafe { fmpq_is_zero(&self.raw as *const _ as *mut _) == 1 }
    }
    /// Check if the coefficient is one
    fn is_one(&self) -> bool {
        unsafe { fmpq_is_one(&self.raw as *const _ as *mut _) == 1 }
    }
}

struct QMPoly {
    //pub raw: MaybeUninit<fmpq_mpoly_struct>,
    pub raw: fmpq_mpoly_struct,     // Polynomial
    pub ctx: fmpq_mpoly_ctx_struct, // Polynomial context
    pub vars: Vec<String>,          // Name of variables
    _coeff1: QCoeff,                // coefficient container
    _coeff2: QCoeff,                // coefficient container
    _res: QCoeff,                   // coefficient container
    _exp: [u64; 10],                // exponent container
}

/// Container for polynomial that uses Flint as backend
impl QMPoly {
    /// New 0 polynomial with selected variables
    fn new(vars: &[String]) -> Self {
        unsafe {
            let mut raw = MaybeUninit::uninit();
            let mut ctx = MaybeUninit::uninit();
            let nvars = vars.len() as i64;
            fmpq_mpoly_ctx_init(ctx.as_mut_ptr(), nvars, ordering_t_ORD_DEGLEX);
            fmpq_mpoly_init(raw.as_mut_ptr(), ctx.as_mut_ptr());
            //fmpq_mpoly_zero(raw.as_mut_ptr(), ctx.as_mut_ptr());
            QMPoly {
                //raw: raw,
                raw: raw.assume_init(),
                ctx: ctx.assume_init(),
                vars: vars.to_vec(),
                _coeff1: QCoeff::new(),
                _coeff2: QCoeff::new(),
                _res: QCoeff::new(),
                _exp: [0; 10],
            }
        }
    }
    /// Clone will copy the pointers and we will have two objects pointing at the same memory
    /// address. With this function we can create a clone in a safe way
    fn clone_from(other: &Self) -> Self {
        let mut this = QMPoly::new(&other.vars);
        unsafe {
            fmpq_mpoly_set(
                &mut this.raw as *mut _,
                &other.raw as *const _ as *mut _,
                &other.ctx as *const _ as *mut _,
            );
        }
        this
    }
    /// Set polynomial to zero
    fn set_to_zero(&mut self) {
        unsafe {
            fmpq_mpoly_zero(&mut self.raw as *mut _, &mut self.ctx as *mut _);
        }
    }
    /// Set polynomial to one
    fn set_to_one(&mut self) {
        unsafe {
            fmpq_mpoly_one(&mut self.raw as *mut _, &mut self.ctx as *mut _);
        }
    }
    /// Check if the polynomial is zero
    fn is_zero(&self) -> bool {
        unsafe {
            fmpq_mpoly_is_zero(
                &self.raw as *const _ as *mut _,
                &self.ctx as *const _ as *mut _,
            ) == 1
        }
    }
    /// Check if the polynomial is one
    fn is_one(&self) -> bool {
        unsafe {
            fmpq_mpoly_is_one(
                &self.raw as *const _ as *mut _,
                &self.ctx as *const _ as *mut _,
            ) == 1
        }
    }
    /// Add to the polynomial the coefficient p/q with the specified powers
    fn add_monomial(&mut self, monimial_powers: &[u64], p: i64, q: u64) {
        unsafe {
            fmpq_mpoly_get_coeff_fmpq_ui(
                &mut self._coeff1.raw as *mut _,
                &mut self.raw as *mut _,
                monimial_powers.as_ptr(),
                &mut self.ctx as *mut _,
            );
            self._coeff2.set_from_int(p, q);
            if self._coeff1.is_zero() {
                fmpq_swap(
                    &mut self._res.raw as *mut _,
                    &mut self._coeff2.raw as *mut _,
                )
            } else {
                fmpq_add(
                    &mut self._res.raw as *mut _,
                    &mut self._coeff1.raw as *mut _,
                    &mut self._coeff2.raw as *mut _,
                );
            }
            fmpq_mpoly_set_coeff_fmpq_ui(
                &mut self.raw as *mut _,
                &mut self._res.raw as *mut _,
                monimial_powers.as_ptr(),
                &mut self.ctx as *mut _,
            );
        }
    }
    /// Convert the polynomial to a human readable string
    fn to_str(&self) -> String {
        unsafe {
            let length = fmpq_mpoly_length(
                &self.raw as *const _ as *mut _,
                &self.ctx as *const _ as *mut _,
            );

            let mut spoly = String::new();
            let nvars = self.vars.len();
            for ci in 0..length {
                // get coeff
                fmpq_mpoly_get_term_coeff_fmpq(
                    &self._coeff1.raw as *const _ as *mut _,
                    &self.raw as *const _ as *mut _,
                    ci,
                    &self.ctx as *const _ as *mut _,
                );
                // get exponent
                fmpq_mpoly_get_term_exp_ui(
                    &self._exp as *const _ as *mut _,
                    &self.raw as *const _ as *mut _,
                    ci,
                    &self.ctx as *const _ as *mut _,
                );

                if fmpq_is_zero(&self._coeff1.raw as *const _ as *mut _) == 0 {
                    std::fmt::write(&mut spoly, format_args!("{}", self._coeff1)).unwrap();

                    for (n, p) in self._exp[..nvars].iter().enumerate() {
                        match *p {
                            0 => continue,
                            1 => {
                                std::fmt::write(&mut spoly, format_args!("*{}", self.vars[n]))
                                    .unwrap();
                            }
                            _ => {
                                std::fmt::write(
                                    &mut spoly,
                                    format_args!("*{}^{}", self.vars[n], p),
                                )
                                .unwrap();
                            }
                        }
                    }
                }
            }
            if spoly.is_empty() {
                spoly = String::from("0");
            }
            spoly
        }
    }
}

/// Container for rational polynomials
struct QMRat {
    pub num: QMPoly,       // numerator
    pub den: QMPoly,       // Polynomial context
    pub vars: Vec<String>, // Name of variables
    _coeff: QCoeff,        // coefficient container
    _gcd: QMPoly,          // gcd container
    _res: QMPoly,          // partial operations container
}

impl QMRat {
    /// Create a rational function equal to 0
    fn new(vars: &[String]) -> Self {
        let num = QMPoly::new(vars);
        let mut den = QMPoly::new(vars);
        den.set_to_one();
        QMRat {
            num,
            den,
            vars: vars.to_vec(),
            _coeff: QCoeff::new(),
            _gcd: QMPoly::new(vars),
            _res: QMPoly::new(vars),
        }
    }
    /// Create a rational function from a polynomial
    fn from_mpoly(mpoly: &QMPoly) -> Self {
        let mut den = QMPoly::new(&mpoly.vars);
        den.set_to_one();
        QMRat {
            num: QMPoly::clone_from(mpoly),
            den,
            vars: mpoly.vars.clone(),
            _coeff: QCoeff::new(),
            _gcd: QMPoly::new(&mpoly.vars),
            _res: QMPoly::new(&mpoly.vars),
        }
    }
    /// Create a rational function from two polynomials
    fn from_mpolys(num: &QMPoly, den: &QMPoly) -> Self {
        QMRat {
            num: QMPoly::clone_from(num),
            den: QMPoly::clone_from(den),
            vars: num.vars.clone(),
            _coeff: QCoeff::new(),
            _gcd: QMPoly::new(&num.vars),
            _res: QMPoly::new(&num.vars),
        }
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
    fn reduce(&mut self) {
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
    }
    /// Normalize to have a monic polynomial at the denominator
    fn normalize(&mut self) {
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
    fn full_reduce(&mut self) {
        self.reduce();
        self.normalize();
    }
}

///Implement addition with operator `+`
impl<'a, 'b> Add<&'b QCoeff> for &'a QCoeff {
    type Output = QCoeff;

    fn add(self, other: &'b QCoeff) -> QCoeff {
        let mut this = QCoeff::new();
        unsafe {
            fmpq_add(
                &mut this.raw as *mut _,
                &self.raw as *const _ as *mut _,
                &other.raw as *const _ as *mut _,
            );
        }
        this
    }
}

///Implement addition with operator `-`
impl<'a, 'b> Sub<&'b QCoeff> for &'a QCoeff {
    type Output = QCoeff;

    fn sub(self, other: &'b QCoeff) -> QCoeff {
        let mut this = QCoeff::new();
        unsafe {
            fmpq_sub(
                &mut this.raw as *mut _,
                &self.raw as *const _ as *mut _,
                &other.raw as *const _ as *mut _,
            );
        }
        this
    }
}
///Implement addition with operator `*`
impl<'a, 'b> Mul<&'b QCoeff> for &'a QCoeff {
    type Output = QCoeff;

    fn mul(self, other: &'b QCoeff) -> QCoeff {
        let mut this = QCoeff::new();
        unsafe {
            fmpq_mul(
                &mut this.raw as *mut _,
                &self.raw as *const _ as *mut _,
                &other.raw as *const _ as *mut _,
            );
        }
        this
    }
}
///Implement addition with operator `/`
impl<'a, 'b> Div<&'b QCoeff> for &'a QCoeff {
    type Output = QCoeff;

    fn div(self, other: &'b QCoeff) -> QCoeff {
        let mut this = QCoeff::new();
        unsafe {
            fmpq_div(
                &mut this.raw as *mut _,
                &self.raw as *const _ as *mut _,
                &other.raw as *const _ as *mut _,
            );
        }
        this
    }
}
///Implement addition with operator `+`
impl<'a, 'b> Add<&'b QMPoly> for &'a QMPoly {
    type Output = QMPoly;

    fn add(self, other: &'b QMPoly) -> QMPoly {
        // TODO: Check if the two contexts are the same
        let mut this = QMPoly::new(&self.vars);
        unsafe {
            fmpq_mpoly_add(
                &mut this.raw as *mut _,
                &self.raw as *const _ as *mut _,
                &other.raw as *const _ as *mut _,
                &other.ctx as *const _ as *mut _,
            );
        }
        this
    }
}
///Implement addition with operator `-`
impl<'a, 'b> Sub<&'b QMPoly> for &'a QMPoly {
    type Output = QMPoly;

    fn sub(self, other: &'b QMPoly) -> QMPoly {
        // TODO: Check if the two contexts are the same
        let mut this = QMPoly::new(&self.vars);
        unsafe {
            fmpq_mpoly_sub(
                &mut this.raw as *mut _,
                &self.raw as *const _ as *mut _,
                &other.raw as *const _ as *mut _,
                &other.ctx as *const _ as *mut _,
            );
        }
        this
    }
}
///Implement addition with operator `*`
impl<'a, 'b> Mul<&'b QMPoly> for &'a QMPoly {
    type Output = QMPoly;

    fn mul(self, other: &'b QMPoly) -> QMPoly {
        // TODO: Check if the two contexts are the same
        let mut this = QMPoly::new(&self.vars);
        unsafe {
            fmpq_mpoly_mul(
                &mut this.raw as *mut _,
                &self.raw as *const _ as *mut _,
                &other.raw as *const _ as *mut _,
                &other.ctx as *const _ as *mut _,
            );
        }
        this
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

/// Clear the content of all raw pointers before droping QMPoly
impl Drop for QMPoly {
    fn drop(&mut self) {
        unsafe {
            fmpq_mpoly_clear(&mut self.raw as *mut _, &mut self.ctx as *mut _);
            fmpq_mpoly_ctx_clear(&mut self.ctx as *mut _);
        }
    }
}
/// Clear the content of all raw pointers before dropping QCoeff
impl Drop for QCoeff {
    fn drop(&mut self) {
        unsafe {
            fmpq_clear(&mut self.raw as *mut _);
        }
    }
}

// To use the `{}` for the structure QCoeff
impl fmt::Display for QCoeff {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.raw.den.0 == 1 {
            write!(f, "{:+}", self.raw.num.0)
        } else {
            write!(f, "{:+}/{}", self.raw.num.0, self.raw.den.0)
        }
    }
}
// To use the `{}` for the structure QMPoly
impl fmt::Display for QMPoly {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.is_zero() {
            return write!(f, "0");
        }
        unsafe {
            let length = fmpq_mpoly_length(
                &self.raw as *const _ as *mut _,
                &self.ctx as *const _ as *mut _,
            );

            let nvars = self.vars.len();
            for ci in 0..length {
                // get coeff
                fmpq_mpoly_get_term_coeff_fmpq(
                    &self._coeff1.raw as *const _ as *mut _,
                    &self.raw as *const _ as *mut _,
                    ci,
                    &self.ctx as *const _ as *mut _,
                );
                // get exponent
                fmpq_mpoly_get_term_exp_ui(
                    &self._exp as *const _ as *mut _,
                    &self.raw as *const _ as *mut _,
                    ci,
                    &self.ctx as *const _ as *mut _,
                );

                if fmpq_is_zero(&self._coeff1.raw as *const _ as *mut _) == 0 {
                    write!(f, "{}", self._coeff1)?;

                    for (n, p) in self._exp[..nvars].iter().enumerate() {
                        match *p {
                            0 => continue,
                            1 => {
                                write!(f, "*{}", self.vars[n])?;
                            }
                            _ => {
                                write!(f, "*{}^{}", self.vars[n], p)?;
                            }
                        }
                    }
                }
            }
        }
        Ok(())
    }
}
// To use the `{}` for the structure QMRat
impl fmt::Display for QMRat {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "({})/({})", self.num, self.den)
    }
}

fn main() {
    // Define variables
    let vars = [String::from("x"), String::from("y")];

    //Set a coefficient to the canonical form of the fraction p / q.
    let mut x = QCoeff::from_int(99, 223);
    println!("x = {}", x);
    x.set_from_int(1, 2);
    println!("x = {}", x);

    // Create Polynomial
    let mut f = QMPoly::new(&vars);
    f.add_monomial(&[1, 1], 3, 4);
    f.add_monomial(&[2, 1], 5, 4);
    f.add_monomial(&[0, 0], 7, 3);
    f.add_monomial(&[0, 0], 7, 3);
    println!("f   = {}", f.to_str());

    // Add Polynomials
    let g = QMPoly::clone_from(&f);
    let h = &f + &g;
    println!("g   = {}", g);
    println!("g+f = {}", h);

    // Subtract Polynomials
    let h = &f - &g;
    println!("g-f = {}", h);

    // Multiply Polynomials
    let h = &f * &g;
    println!("g*f = {}", h);

    // Add Coefficients
    let c1 = QCoeff::from_int(2, 3);
    let c2 = QCoeff::from_int(3, 5);
    let c3 = &c1 + &c2;
    assert_eq!(c3.raw.num.0, 19);
    assert_eq!(c3.raw.den.0, 15);

    // Subtract Coefficients
    let c3 = &c1 - &c2;
    assert_eq!(c3.raw.num.0, 1);
    assert_eq!(c3.raw.den.0, 15);

    // Multiply Coefficients
    let c3 = &c1 * &c2;
    assert_eq!(c3.raw.num.0, 2);
    assert_eq!(c3.raw.den.0, 5);

    // Divide Coefficients
    let c3 = &c1 * &c2;
    assert_eq!(c3.raw.num.0, 2);
    assert_eq!(c3.raw.den.0, 5);

    // ==================================
    // Rational functions
    let mut r1 = QMRat::from_mpolys(&f, &h);
    let mut r2 = QMRat::from_mpolys(&f, &(&g + &h));
    println!("r1 = {}", r1);
    println!("r2 = {}", r2);
//    r1.full_reduce();
//    println!("r1 -> gcd = {}", r1._gcd.to_str());
//    println!("r1 -> reduce = {}", r1);
//    r2.full_reduce();
//    println!("r2 -> gcd = {}", r2._gcd.to_str());
//    println!("r2 -> reduce = {}", r2);

    let mut r3 = &r1/&r2;
    println!("r1/r2 = {}", r3);
    r3.full_reduce();
    println!("r1/r2 = {}", r3);

    println!("Hello, world!");
}
