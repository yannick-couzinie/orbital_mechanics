//! Bisection algorithm for finding the root

use super::errors::RootError;
use super::structs::RootResult;

fn evaluate(func: &impl Fn(f64) -> f64, input: f64) -> Result<f64, RootError> {
    let output = func(input);

    if output.is_finite() {
        Ok(output)
    } else {
        Err(RootError::NonFiniteEvaluation { input, output })
    }
}

pub fn bisect(
    func: impl Fn(f64) -> f64,
    initial_lower: f64,
    initial_upper: f64,
    tolerance: f64,
) -> Result<RootResult, RootError> {
    if tolerance <= 0.0 || !tolerance.is_finite() {
        return Err(RootError::InvalidTolerance { tolerance });
    }

    if !initial_lower.is_finite() || !initial_upper.is_finite() {
        return Err(RootError::InvalidSearchInterval {
            lower: initial_lower,
            upper: initial_upper,
        });
    }

    let iterations: usize =
        (((initial_upper - initial_lower).abs() / tolerance).log2()).ceil() as usize;
    let mut lower: f64 = initial_lower;
    let mut upper: f64 = initial_upper;
    let mut middle: f64 = (lower + upper) / 2.0;

    let mut flower = evaluate(&func, initial_lower)?;
    let fupper = evaluate(&func, initial_upper)?;
    let mut fmiddle = evaluate(&func, middle)?;

    if (flower < 0.0 && fupper < 0.0) || (flower > 0.0 && fupper > 0.0) {
        return Err(RootError::RootNotBracketed { lower, upper });
    }

    if flower.abs() == 0.0 {
        return Ok(RootResult {
            root: initial_lower,
            residual: flower.abs(),
            iterations: 0,
        });
    }

    if fupper == 0.0 {
        return Ok(RootResult {
            root: initial_upper,
            residual: fupper.abs(),
            iterations: 0,
        });
    }

    for _ in 0..iterations {
        if fmiddle == 0.0 {
            return Ok(RootResult {
                root: middle,
                residual: fmiddle.abs(),
                iterations,
            });
        }
        if (flower >= 0.0 && fmiddle >= 0.0) || (flower <= 0.0 && fmiddle <= 0.0) {
            lower = middle;
            flower = fmiddle;
            middle = (lower + upper) / 2.0;
            fmiddle = evaluate(&func, middle)?;
        } else {
            upper = middle;
            middle = (lower + upper) / 2.0;
            fmiddle = evaluate(&func, middle)?;
        }
    }

    Ok(RootResult {
        root: middle,
        residual: fmiddle.abs(),
        iterations,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_utils::assert_approx_eq;

    #[test]
    fn test_root_in_interval() {
        let result = bisect(|x| x.powi(2) - 1.0, 0.0, 2.0, 1e-10).unwrap();
        assert_approx_eq(result.root, 1.0, 1e-9, "test_root_in_interval");
    }

    #[test]
    fn test_root_at_lower_border() {
        let result = bisect(|x| x.ln(), 1.0, 2.0, 1e-10).unwrap();
        assert_approx_eq(result.root, 1.0, 1e-9, "test_root_at_lower_border");
    }

    #[test]
    fn test_root_at_upper_border() {
        let result = bisect(|x| x.ln(), 0.5, 1.0, 1e-10).unwrap();
        assert_approx_eq(result.root, 1.0, 1e-9, "test_root_at_upper_border");
    }

    #[test]
    fn test_root_at_midpoint() {
        let result = bisect(|x| x.powi(3), -1.0, 1.0, 1e-10).unwrap();
        assert_approx_eq(result.root, 0.0, 1e-9, "test_root_at_midpoint");
    }

    #[test]
    fn test_not_bracketing_root() {
        assert!(matches!(
            bisect(|x| x.powi(2) + 1.0, -1.0, 1.0, 1.0e-6),
            Err(RootError::RootNotBracketed { .. })
        ));
    }

    #[test]
    fn test_infinite_function() {
        assert!(matches!(
            bisect(|x| x.signum() * f64::INFINITY, -1.0, 1.0, 1.0e-6),
            Err(RootError::NonFiniteEvaluation { .. })
        ));
    }

    #[test]
    fn test_invalid_tolerance() {
        assert!(matches!(
            bisect(|x| x.powi(2) - 1.0, 0.0, 1.0, -1.0e-6),
            Err(RootError::InvalidTolerance { .. })
        ));
    }
}
