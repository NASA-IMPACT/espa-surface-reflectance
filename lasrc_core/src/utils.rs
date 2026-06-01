//! Utility functions ported from quick_select.c and poly_coeff.c in C LaSRC.

use crate::constants::NCOEF;

/// Find the median of a floating point array using the quickselect algorithm.
///
/// Based on "Numerical Recipes in C", Second Edition, Section 8.5.
/// This code by Nicolas Devillard - 1998. Public domain.
///
/// Modifies the input array order as a side effect.
/// Returns the element at index n/2 (lower-median for even-length arrays).
pub fn quick_select(arr: &mut [f32]) -> f32 {
    let n = arr.len();
    if n == 0 {
        panic!("quick_select called on empty slice");
    }

    let mut low: usize = 0;
    let mut high: usize = n - 1;
    let median: usize = n / 2;

    loop {
        // One element
        if high <= low {
            return arr[median];
        }

        // Two elements
        if high == low + 1 {
            if arr[low] > arr[high] {
                arr.swap(low, high);
            }
            return arr[median];
        }

        // Median-of-three: sort arr[middle], arr[high], arr[low] pairwise
        let middle = (low + high) / 2;
        if arr[middle] > arr[high] {
            arr.swap(middle, high);
        }
        if arr[low] > arr[high] {
            arr.swap(low, high);
        }
        if arr[middle] > arr[low] {
            arr.swap(middle, low);
        }

        // Swap low item (now in position middle) into position low+1
        arr.swap(middle, low + 1);

        // Partition: nibble from each end towards middle
        let mut ll = low + 1;
        let mut hh = high;
        loop {
            loop {
                ll += 1;
                if arr[low] <= arr[ll] {
                    break;
                }
            }
            loop {
                hh -= 1;
                if arr[hh] <= arr[low] {
                    break;
                }
            }

            if hh < ll {
                break;
            }

            arr.swap(ll, hh);
        }

        // Swap middle item (in position low) back into correct position
        arr.swap(low, hh);

        // Re-set active partition
        if hh <= median {
            low = ll;
        }
        if hh >= median {
            high = hh - 1;
        }
    }
}

/// Fit a cubic polynomial to the given data points via LU decomposition.
///
/// Builds the normal equations Z'Z * coeff = Z'y using the design matrix
/// where x[i][j] = aot[i]^(3-j) for j=0..3, then solves with LU decomposition
/// and partial pivoting.
///
/// Uses f32 arithmetic internally to match the C implementation (which uses
/// `float` throughout `poly_coeff.c`), then casts back to f64 for the return.
///
/// Returns `[a3, a2, a1, a0]` for the polynomial `a3*x^3 + a2*x^2 + a1*x + a0`.
pub fn get_3rd_order_poly_coeff(aot: &[f64], atm: &[f64]) -> [f64; NCOEF] {
    let n_atm = aot.len().min(atm.len());

    // Cast inputs to f32 to match C's float arithmetic
    let aot_f32: Vec<f32> = aot[..n_atm].iter().map(|&v| v as f32).collect();
    let atm_f32: Vec<f32> = atm[..n_atm].iter().map(|&v| v as f32).collect();

    // Build design matrix x[i][j] = aot[i]^(3-j)
    let mut x = vec![[0.0f32; NCOEF]; n_atm];
    for i in 0..n_atm {
        x[i][0] = aot_f32[i] * aot_f32[i] * aot_f32[i];
        x[i][1] = aot_f32[i] * aot_f32[i];
        x[i][2] = aot_f32[i];
        x[i][3] = 1.0;
    }

    // Build normal equations: z = X'X (NCOEF x NCOEF), y1 = X'y (NCOEF)
    let mut z = [[0.0f32; NCOEF]; NCOEF];
    let mut y1 = [0.0f32; NCOEF];

    for i in 0..NCOEF {
        for j in 0..NCOEF {
            for k in 0..n_atm {
                z[i][j] += x[k][i] * x[k][j];
            }
        }
        for j in 0..n_atm {
            y1[i] += x[j][i] * atm_f32[j];
        }
    }

    // LU decompose with partial pivoting
    let mut perm = [0usize; NCOEF];
    lup_decompose_f32(&mut z, &mut perm, 1e-10);

    // Solve for coefficients
    let mut coeff_f32 = [0.0f32; NCOEF];
    lup_solve_f32(&z, &perm, &y1, &mut coeff_f32);

    // Cast back to f64 for the return type
    [
        coeff_f32[0] as f64,
        coeff_f32[1] as f64,
        coeff_f32[2] as f64,
        coeff_f32[3] as f64,
    ]
}

/// LU decomposition with partial pivoting (f32, matching C's float).
fn lup_decompose_f32(a: &mut [[f32; NCOEF]; NCOEF], perm: &mut [usize; NCOEF], tol: f32) {
    let n = NCOEF;

    for i in 0..n {
        perm[i] = i;
    }

    for i in 0..n {
        let mut max_a = 0.0f32;
        let mut imax = i;
        for k in i..n {
            let abs_a = a[k][i].abs();
            if abs_a > max_a {
                max_a = abs_a;
                imax = k;
            }
        }

        if max_a < tol {
            continue;
        }

        if imax != i {
            perm.swap(i, imax);
            a.swap(i, imax);
        }

        for j in i + 1..n {
            a[j][i] /= a[i][i];
            for k in i + 1..n {
                a[j][k] -= a[j][i] * a[i][k];
            }
        }
    }
}

/// Forward/back substitution (f32 version, matching C's float).
fn lup_solve_f32(
    a: &[[f32; NCOEF]; NCOEF],
    perm: &[usize; NCOEF],
    b: &[f32; NCOEF],
    x: &mut [f32; NCOEF],
) {
    let n = NCOEF;

    for i in 0..n {
        x[i] = b[perm[i]];
        for j in 0..i {
            x[i] -= a[i][j] * x[j];
        }
    }

    for i in (0..n).rev() {
        for j in i + 1..n {
            x[i] -= a[i][j] * x[j];
        }
        x[i] /= a[i][i];
    }
}

/// Evaluate a cubic polynomial at x.
///
/// Returns `coeff[0]*x^3 + coeff[1]*x^2 + coeff[2]*x + coeff[3]`.
pub fn eval_polynomial(coeff: &[f64; NCOEF], x: f64) -> f64 {
    coeff[0] * x * x * x + coeff[1] * x * x + coeff[2] * x + coeff[3]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_quick_select_odd() {
        let mut arr = vec![9.0, 1.0, 5.0, 3.0, 7.0];
        assert_eq!(quick_select(&mut arr), 5.0);
    }

    #[test]
    fn test_quick_select_even() {
        let mut arr = vec![4.0, 2.0, 6.0, 8.0];
        let med = quick_select(&mut arr);
        assert_eq!(med, 6.0); // index n/2 = 2, which is 6 when sorted
    }

    #[test]
    fn test_quick_select_single() {
        let mut arr = vec![42.0];
        assert_eq!(quick_select(&mut arr), 42.0);
    }

    #[test]
    fn test_quick_select_two() {
        let mut arr = vec![10.0, 5.0];
        let med = quick_select(&mut arr);
        assert_eq!(med, 10.0); // index 1 of sorted [5, 10]
    }

    #[test]
    fn test_poly_coeff_linear() {
        let aot: Vec<f64> = (0..10).map(|i| i as f64 * 0.5).collect();
        let atm: Vec<f64> = aot.iter().map(|&x| 2.0 * x + 1.0).collect();
        let coeff = get_3rd_order_poly_coeff(&aot, &atm);
        // Tolerances relaxed to 1e-4 because polynomial fitting now uses f32
        // arithmetic (matching C's float precision).
        assert!(coeff[0].abs() < 1e-4, "a3 should be ~0: {}", coeff[0]);
        assert!(coeff[1].abs() < 1e-4, "a2 should be ~0: {}", coeff[1]);
        assert!((coeff[2] - 2.0).abs() < 1e-4, "a1 should be ~2: {}", coeff[2]);
        assert!((coeff[3] - 1.0).abs() < 1e-4, "a0 should be ~1: {}", coeff[3]);
    }

    #[test]
    fn test_eval_polynomial() {
        let coeff = [1.0, -2.0, 3.0, 4.0];
        assert!((eval_polynomial(&coeff, 0.0) - 4.0).abs() < 1e-10);
        assert!((eval_polynomial(&coeff, 1.0) - 6.0).abs() < 1e-10);
        assert!((eval_polynomial(&coeff, 2.0) - 10.0).abs() < 1e-10);
    }
}
