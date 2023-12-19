use flint_mpoly::*;
//use flint_mpoly::parser::*;

use std::str::FromStr;

use flint_sys::fmpq::*;
use flint_sys::fmpz::*;
use num::BigInt;

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
    f.add_coeff_pq(&[1, 1], 3, 4);
    f.add_coeff_pq(&[2, 1], 5, 4);
    f.add_coeff_pq(&[0, 0], 7, 3);
    f.add_coeff(&[3, 1], &x);
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
    println!();

    // When constructing a new rational expression there is no automatic simplifycation
    let mut r1 = QMRat::from_mpolys(&f, &h);
    let mut r2 = QMRat::from_mpolys(&f, &(&g + &h));
    println!("r1 = {}", r1);
    println!("r2 = {}", r2);

    let mut r3 = &r1 / &r2;
    println!("r1/r2 = {}", r3);
    r3.full_reduce();
    println!("r1/r2 = {}", r3);

    // To simplify the original expression we need to call for full_reduce()
    r1.full_reduce();
    r2.full_reduce();
    println!("r1 -> reduced = {}", r1);
    println!("r2 -> reduced = {}", r2);

    // ==================================
    // Read From string
    println!();

    let coeff_str = "1391941";
    let c = QCoeff::from_str(coeff_str).expect("bad string");
    println!("read from str: {}", coeff_str);
    println!(" - c.num = {}", c.raw.num.0);
    println!(" - c.den = {}", c.raw.den.0);
    println!(" => c = {}", c);
    let coeff_str = "1391941/218653";
    let c = QCoeff::from_str(coeff_str).expect("bad string");
    println!("read from str: {}", coeff_str);
    println!(" - c.num = {}", c.raw.num.0);
    println!(" - c.den = {}", c.raw.den.0);
    println!(" => c = {}", c);
    let coeff_str = "-1234567899876543219991823/218653";
    let c = QCoeff::from_str(coeff_str).expect("bad string");
    println!("read from str: {}", coeff_str);
    println!(" - c.num = {}", c.raw.num.0);
    println!(" - c.den = {}", c.raw.den.0);
    println!(" => c = {}", c);

    // Create Polynomial
    let mut f = QMPoly::new(&vars);
    f.add_coeff_pq(&[0, 0], 1, 3);
    f.add_coeff_str(&[10, 7], coeff_str);
    println!("f={}", f);
    //    unsafe {
    //        println!("{}", fmpz_size(&mut c.raw.num as *mut _));
    //        println!("{}", fmpz_bits(&mut c.raw.num as *mut _));
    //        println!("{}", fmpz_sizeinbase(&mut c.raw.num as *mut _, 10));
    //        let mut n = BigInt::default();
    //        let sgn = (1 + fmpz_sgn(&mut c.raw.num as *mut _)) / 2;
    //        println!("sgn = {}", sgn);
    //        for i in 0..fmpz_bits(&mut c.raw.num as *mut _) {
    //            n.set_bit(i, fmpz_tstbit(&mut c.raw.num as *mut _, i) == sgn);
    //        }
    //        if sgn == 0 {
    //            n.set_bit(0, true);
    //            n = -n;
    //        }
    //        println!("{}", n);
    //    }

    println!("Hello, world!");
}
