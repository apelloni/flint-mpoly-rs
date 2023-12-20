# flint-mpoly-rs
multivariate polynomial for rust with FLINT backend


## Coefficients
 - [x] `QCoeff`: Coefficients that take values in $\mathbb{Q}$ (from: `fmpq.h`)
 - [ ] `ZCoeff`: Coefficients that take values in $\mathbb{Z}$ (from: `fmpz.h`)
 - [ ] `ModCoeff`: Coefficients that take values in $\mathbb{Z}/p\mathbb{Z}$ (from: `fmpz_mod.h`)

## Polynomials
 - [x] `QMPoly`: Multivariate Polynomials in $\mathbb{Q}[x_1,x_2,..]$ (from: `fmpq_mpoly.h`)
 - [ ] `ZMPoly`: Multivariate Polynomials in $\mathbb{Z}[x_1,x_2,..]$ (from: `fmpz_mpoly.h`)
 - [ ] `ModMPoly`: Multivariate Polynomials in $\mathbb{Z}/p\mathbb{Z}[x_1,x_2,..]$ (from: `fmpz_mod_mpoly.h`)
 - [ ] `QPoly`: Univariate Polynomials in $\mathbb{Q}[x]$ (from: `fmpq_poly.h`)
 - [ ] `ZPoly`: Univariate Polynomials in $\mathbb{Z}[x]$ (from: `fmpz_poly.h`)
 - [ ] `ModPoly`: Univariate Polynomials in $\mathbb{Z}/p\mathbb{Z}[x]$ (from: `fmpz_mod_poly.h`)

## Rational Polynomials
 - [x] `QMRat`: Multivariate Rational in $\mathbb{Q}(x_1,x_2,..)$
 - [ ] `ZMRat`: Multivariate Rational in $\mathbb{Z}(x_1,x_2,..)$
 - [ ] `ModMRat`: Multivariate Rational in $\mathbb{Z}/p\mathbb{Z}(x_1,x_2,..)$
 - [ ] `QRat`: Univariate Rational in $\mathbb{Q}(x)$
 - [ ] `ZRat`: Univariate Rational in $\mathbb{Z}(x)$
 - [ ] `ModRat`: Univariate Rational in $\mathbb{Z}/p\mathbb{Z}(x)$
