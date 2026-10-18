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

    let flower = evaluate(&func, initial_lower)?;
    let fupper = evaluate(&func, initial_upper)?;
    if (flower <= 0.0 && fupper <= 0.0) || (flower >= 0.0 && fupper >= 0.0) {
        return Err(RootError::RootNotBracketed { lower, upper });
    }

    if flower.abs() <= tolerance {
        return Ok(RootResult {
            root: initial_lower,
            residual: flower.abs(),
            iterations: 0,
        });
    }

    if fupper.abs() <= tolerance {
        return Ok(RootResult {
            root: initial_upper,
            residual: fupper.abs(),
            iterations: 0,
        });
    }

    for _ in 0..iterations {
        let flower = evaluate(&func, lower)?;

        // at first iteration this is calculated twice but we need middle in the scope outside to
        // return it
        middle = (lower + upper) / 2.0;

        if !flower.is_finite() {
            return Err(RootError::NonFiniteEvaluation {
                input: lower,
                output: flower,
            });
        }

        let fmiddle = evaluate(&func, middle)?;
        if !fmiddle.is_finite() {
            return Err(RootError::NonFiniteEvaluation {
                input: middle,
                output: fmiddle,
            });
        }

        if (flower >= 0.0 && fmiddle >= 0.0) || (flower <= 0.0 && fmiddle <= 0.0) {
            lower = middle;
        } else {
            upper = middle;
        }
    }

    Ok(RootResult {
        root: middle,
        residual: evaluate(&func, middle)?.abs(),
        iterations,
    })
}
