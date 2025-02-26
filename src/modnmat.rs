use flint_sys::nmod_mat::*;
use std::mem::MaybeUninit;

/// Container for polynomials over rational numbers
#[derive(Debug)]
pub struct NModMat {
    pub raw: nmod_mat_struct, // Polynomial
    pub nmod: u64,            // mod for the integer
}

/// Container for polynomial that uses Flint as backend
impl NModMat {
    /// New 0 matrix of given size
    /// ```
    /// use flint_mpoly::NModMat;
    /// let cols = 10;
    /// let rows = 10;
    /// let mut mat = NModMat::new(cols,cols,3);
    /// for i in 0..rows{
    ///     for j in 0..cols{
    ///         assert_eq!(0,mat.get_entry(i,j));
    ///     }
    /// }
    /// ```
    pub fn new(rows: i64, cols: i64, nmod: u64) -> Self {
        unsafe {
            let mut raw = MaybeUninit::uninit();
            nmod_mat_init(raw.as_mut_ptr(), rows, cols, nmod);
            Self {
                //raw: raw,
                raw: raw.assume_init(),
                nmod,
            }
        }
    }

    /// Set enty (i,j) to given value
    /// ```
    /// use flint_mpoly::NModMat;
    /// let cols = 10;
    /// let rows = 10;
    /// let mut mat = NModMat::new(cols,cols,3);
    /// for i in 0..rows{
    ///     for j in 0..cols{
    ///         mat.set_entry(i, j, (i*j) as u64);
    ///     }
    /// }
    /// for i in 0..rows{
    ///     for j in 0..cols{
    ///         assert_eq!(i*j, mat.get_entry(i,j) as i64);
    ///     }
    /// }
    /// ```
    pub fn set_entry(&mut self, row: i64, col: i64, val: u64) {
        unsafe {
            nmod_mat_set_entry(&mut self.raw as *mut _, row, col, val);
        }
    }
    pub fn get_entry(&mut self, row: i64, col: i64) -> u64 {
        unsafe { nmod_mat_get_entry(&mut self.raw as *mut _, row, col) }
    }

    /// Set all entry of the matrix to 0
    /// ```
    /// use flint_mpoly::NModMat;
    /// let cols = 10;
    /// let rows = 10;
    /// let mut mat = NModMat::new(cols,cols,3);
    /// for i in 0..rows{
    ///     for j in 0..cols{
    ///         mat.set_entry(i, j, (i*j) as u64);
    ///     }
    /// }
    /// mat.to_zero();
    /// for i in 0..rows{
    ///     for j in 0..cols{
    ///         assert_eq!(0, mat.get_entry(i,j) as i64);
    ///     }
    /// }
    /// ```
    pub fn to_zero(&mut self) {
        unsafe {
            nmod_mat_zero(&mut self.raw as *mut _);
        }
    }

    /// Tranform matrix in Reduced Row Echelon Form (RREF)
    /// ```
    /// use flint_mpoly::NModMat;
    /// let cols = 10;
    /// let rows = 10;
    /// let module = 3;
    /// let mut mat = NModMat::new(cols,cols,module);
    /// for i in 0..rows{
    ///     for j in 0..cols{
    ///         mat.set_entry(i, j, (i+j) as u64 % module);
    ///     }
    /// }
    /// for i in 0..rows{
    ///     for j in 0..cols{
    ///         print!("{: >width$}", mat.get_entry(i,j), width=4);
    ///         //assert_eq!(0, mat.get_entry(i,j) as i64);
    ///     }
    ///     println!("\n");
    /// }
    ///
    /// let rank = mat.rref();
    /// assert_eq!(rank, 2);
    /// //for i in 0..rank{
    /// //    for j in 0..mat.ncols(){
    /// //        print!("{: >width$}", mat.get_entry(i,j), width=4);
    /// //        //assert_eq!(0, mat.get_entry(i,j) as i64);
    /// //    }
    /// //    println!("\n");
    /// //}
    /// ```
    pub fn rref(&mut self) -> i64 {
        unsafe { nmod_mat_rref(&mut self.raw as *mut _) }
    }

    /// Compute rank of matrix
    /// ```
    /// use flint_mpoly::NModMat;
    /// let cols = 10;
    /// let rows = 10;
    /// let module = 3;
    /// let mut mat = NModMat::new(cols,cols,module);
    /// for i in 0..rows{
    ///     for j in 0..cols{
    ///         mat.set_entry(i, j, (i+j) as u64 % module);
    ///     }
    /// }
    /// for i in 0..rows{
    ///     for j in 0..cols{
    ///         print!("{: >width$}", mat.get_entry(i,j), width=4);
    ///         //assert_eq!(0, mat.get_entry(i,j) as i64);
    ///     }
    ///     println!("\n");
    /// }
    ///
    /// let rank = mat.rank();
    /// assert_eq!(rank, 2);
    /// //for i in 0..rank{
    /// //    for j in 0..mat.ncols(){
    /// //        print!("{: >width$}", mat.get_entry(i,j), width=4);
    /// //        //assert_eq!(0, mat.get_entry(i,j) as i64);
    /// //    }
    /// //    println!("\n");
    /// //}
    /// ```
    pub fn rank(&mut self) -> i64 {
        unsafe { nmod_mat_rank(&mut self.raw as *mut _) }
    }

    /// Get Number of cols
    pub fn ncols(&mut self) -> i64 {
        unsafe { nmod_mat_ncols(&mut self.raw as *mut _) }
    }
    /// Get Number of rows
    pub fn nrows(&mut self) -> i64 {
        unsafe { nmod_mat_nrows(&mut self.raw as *mut _) }
    }
    /// Check if matrix is zero
    pub fn is_zero(&mut self) -> bool {
        unsafe { nmod_mat_is_zero(&mut self.raw as *mut _) > 0 }
    }
    /// Check if row is zero
    pub fn is_zero_row(&mut self, row: i64) -> bool {
        unsafe { nmod_mat_is_zero_row(&mut self.raw as *mut _, row) > 0 }
    }
}

/// Clear the content of all raw pointers before dropping NModMat
impl Drop for NModMat {
    fn drop(&mut self) {
        unsafe {
            nmod_mat_clear(&mut self.raw as *mut _);
            //flint_cleanup();
        }
    }
}

// Send and Sync
unsafe impl Send for NModMat {}
unsafe impl Sync for NModMat {}
