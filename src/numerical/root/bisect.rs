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

    let mut flower = evaluate(&func, initial_lower)?;
    let fupper = evaluate(&func, initial_upper)?;

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

    let mut middle: f64 = (lower + upper) / 2.0;
    let mut fmiddle = evaluate(&func, middle)?;

    for cur_iter in 0..iterations {
        if fmiddle == 0.0 {
            return Ok(RootResult {
                root: middle,
                residual: fmiddle.abs(),
                iterations: cur_iter,
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
        let result = bisect(|x| x.powi(2) - 2.0, 1.0, 2.0, 1e-10).unwrap();
        assert_approx_eq(result.root, 2.0_f64.sqrt(), 1e-9, "test_root_in_interval");
    }

    #[test]
    fn test_root_of_transcendental() {
        // This is Problem 3.2: Find the three positive roots of the equation
        // 10 exp(sin(x)) = x^2 - 5x + 4
        // So I am looking for the roots of f(x) = x^2 - 5x + 4 - 10 exp(sin(x))

        // In my calculator I made a table from 0 to 20 with step width 1 and find
        // f(4) = -4.691
        // f(5) = 0.1669
        // f(6) = 2.4377
        // f(7) = -1.289
        // f(8) = 1.1049
        fn func(x: f64) -> f64 {
            x.powi(2) - 5.0 * x + 4.0 - 10.0 * x.sin().exp()
        }
        let result = bisect(func, 4.0, 5.0, 1e-10).unwrap();
        assert_approx_eq(result.root, 4.95777677, 1e-8, "test_transcendental_1");

        let result = bisect(func, 6.0, 7.0, 1e-10).unwrap();
        assert_approx_eq(result.root, 6.77321280, 1e-8, "test_transcendental_2");

        let result = bisect(func, 7.0, 8.0, 1e-10).unwrap();
        assert_approx_eq(result.root, 7.91978355, 1e-8, "test_transcendental_3");
    }

    #[test]
    fn test_root_of_transcendental_tanh_tan() {
        // This is Problem 3.3: Find the first four nonnegative roots of the equation
        // tan(x) = tanh(x) -> 0 = tanh(x) - tan(x)

        // Besides 0 which is an obvious 0 let's find the others
        // In my calculator I made a table from 0 to pi/2 but there was no sign change
        // from pi/2 to 3pi/2:
        // f(3.7707) = 0.271
        // f(3.9707) = -0.092
        //
        // from 3pi/2 to 5pi/2:
        // f(6.9123) = 0.2721
        // f(7.1123) = -0.091

        // from 5pi/2 to 7pi/2:
        // f(10.053) = 0.2736
        // f(10.253) = -0.0895
        fn func(x: f64) -> f64 {
            x.tanh() - x.tan()
        }
        let result = bisect(func, 0.0, 1.0, 1e-10).unwrap();
        assert_approx_eq(result.root, 0.0, 1e-8, "test_transcendental_1");

        let result = bisect(func, 3.7707, 3.9707, 1e-10).unwrap();
        assert_approx_eq(result.root, 3.92660231, 1e-8, "test_transcendental_1");

        let result = bisect(func, 6.9123, 7.1123, 1e-10).unwrap();
        assert_approx_eq(result.root, 7.06858275, 1e-8, "test_transcendental_2");

        let result = bisect(func, 10.053, 10.253, 1e-10).unwrap();
        assert_approx_eq(result.root, 10.21017612, 1e-8, "test_transcendental_3");
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
