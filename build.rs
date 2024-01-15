use std::env;
use std::ffi::OsString;
use std::path::PathBuf;

const FLINT_VER: &str = ">3.0.1";
const FLINT_LIB: &str = "libflint.a";
const FLINT_DIR: &str = "lib/";

fn main() {
    // Pass `-fopenmp` to the linker.
    //println!("cargo:rustc-link-arg=-fkkkkk");
    // Get library from gmp-mpfr-sys crate
    let gmp_mpfr_dir = PathBuf::from(cargo_env("DEP_GMP_OUT_DIR"));

    println!("cargo:rustc-link-search={}", gmp_mpfr_dir.display());
    println!("cargo:rustc-link-search={}", FLINT_DIR);
    println!("cargo:rustc-link-lib=static=flint");
    println!("cargo:rustc-link-lib=static=gmp");
    println!("cargo:rustc-link-lib=static=mpfr");
}

fn cargo_env(name: &str) -> OsString {
    env::var_os(name)
        .unwrap_or_else(|| panic!("environment variable not found: {}, please use cargo", name))
}
