//! Newton algorithm for finding the root

use super::errors::RootError;
use super::structs::RootResult;

const MAX_ITERATIONS: usize = 1000;

fn evaluate(func: &impl Fn(f64) -> f64, input: f64, derivative: bool) -> Result<f64, RootError> {
    let output = func(input);

    if output.is_finite() {
        Ok(output)
    } else {
        if derivative {
            Err(RootError::NonFiniteDerivative { input, output })
        } else {
            Err(RootError::NonFiniteEvaluation { input, output })
        }
    }
}

pub fn newton(
    func: impl Fn(f64) -> f64,
    deriv: impl Fn(f64) -> f64,
    initial_guess: f64,
    tolerance: f64,
) -> Result<RootResult, RootError> {
    if tolerance <= 0.0 || !tolerance.is_finite() {
        return Err(RootError::InvalidTolerance { tolerance });
    }

    if !initial_guess.is_finite() {
        return Err(RootError::InvalidInitialValue { initial_guess });
    }

    let mut ratio = f64::INFINITY;
    let mut current_guess = initial_guess;
    let mut iterations = 0;
    let mut fval = evaluate(&func, current_guess, false)?;

    while ratio.abs() > tolerance {
        if fval == 0.0 {
            break;
        }

        let deriv_val = evaluate(&deriv, current_guess, true)?;

        if deriv_val == 0.0 {
            return Err(RootError::ZeroDerivative {
                input: current_guess,
                fval,
                deriv_val,
            });
        }
        ratio = fval / deriv_val;

        if ratio.abs() < tolerance {
            break;
        }
        current_guess -= ratio;
        fval = evaluate(&func, current_guess, false)?;

        iterations += 1;
        if iterations >= MAX_ITERATIONS {
            return Err(RootError::DidNotConverge {
                current_guess,
                fval,
            });
        }
    }

    Ok(RootResult {
        root: current_guess,
        residual: fval.abs(),
        iterations,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_utils::assert_approx_eq;

    #[test]
    fn test_find_root() {
        let result = newton(|x| x.powi(2) - 2.0, |x| 2.0 * x, 1.4, 1e-10).unwrap();
        assert_approx_eq(result.root, 2.0_f64.sqrt(), 1e-9, "test_root_in_interval");
    }

    #[test]
    fn test_root_of_transcendental() {
        // Same as for bisection
        fn func(x: f64) -> f64 {
            x.powi(2) - 5.0 * x + 4.0 - 10.0 * x.sin().exp()
        }

        fn deriv(x: f64) -> f64 {
            2.0 * x - 5.0 - 10.0 * x.sin().exp() * x.cos()
        }
        let result = newton(func, deriv, 4.9, 1e-10).unwrap();
        assert_approx_eq(result.root, 4.95777677, 1e-8, "test_transcendental_1");

        let result = newton(func, deriv, 6.7, 1e-10).unwrap();
        assert_approx_eq(result.root, 6.77321280, 1e-8, "test_transcendental_2");

        let result = newton(func, deriv, 7.9, 1e-10).unwrap();
        assert_approx_eq(result.root, 7.91978355, 1e-8, "test_transcendental_3");
    }

    #[test]
    fn test_root_of_transcendental_tanh_tan() {
        // Same as for bisection
        fn func(x: f64) -> f64 {
            x.tanh() - x.tan()
        }
        fn deriv(x: f64) -> f64 {
            1.0 - x.tanh().powi(2) - 1.0 / x.cos().powi(2)
        }
        let result = newton(func, deriv, 3.9707, 1e-10).unwrap();
        assert_approx_eq(result.root, 3.92660231, 1e-8, "test_transcendental_1");

        let result = newton(func, deriv, 7.1123, 1e-10).unwrap();
        assert_approx_eq(result.root, 7.06858275, 1e-8, "test_transcendental_2");

        let result = newton(func, deriv, 10.253, 1e-10).unwrap();
        assert_approx_eq(result.root, 10.21017612, 1e-8, "test_transcendental_3");
    }

    #[test]
    fn test_init_guess_is_root() {
        let result = newton(|x| x.ln(), |x| 1.0 / x, 1.0, 1e-10).unwrap();
        assert_approx_eq(result.root, 1.0, 1e-9, "test_init_guess_is_root");
    }

    #[test]
    fn test_infinite_function() {
        assert!(matches!(
            newton(
                |x| x.signum() * f64::INFINITY,
                |x| x.signum() * f64::INFINITY,
                1.0,
                1.0e-6
            ),
            Err(RootError::NonFiniteEvaluation { .. })
        ));
    }

    #[test]
    fn test_infinite_derivative() {
        assert!(matches!(
            newton(|x| x.signum(), |x| x.signum() * f64::INFINITY, 1.0, 1.0e-6),
            Err(RootError::NonFiniteDerivative { .. })
        ));
    }

    #[test]
    fn test_invalid_tolerance() {
        assert!(matches!(
            newton(|x| x.powi(2) - 1.0, |x| 2.0 * x, 1.0, -1.0e-6),
            Err(RootError::InvalidTolerance { .. })
        ));
    }

    #[test]
    fn test_did_not_converge() {
        assert!(matches!(
            newton(
                |x| x.powi(3) - 2.0 * x + 2.0,
                |x| 3.0 * x.powi(2) - 2.0,
                0.0,
                1.0e-6
            ),
            Err(RootError::DidNotConverge { .. })
        ));
    }

    #[test]
    fn test_zero_derivative() {
        assert!(matches!(
            newton(|x| x.powi(3) - 2.0 * x + 2.0, |_| 0.0, 0.0, 1.0e-6),
            Err(RootError::ZeroDerivative { .. })
        ));
    }

    #[test]
    fn test_problem_39a() {
        //! An earth-orbiting satellite has a perigee radius of 7000km and an apogee radius of
        //! 10_000km. What true anomaly delta_theta is swept out between t = 0.5h and t = 1.5h after
        //! perigee passage?
        //! Answer should be 128.7 degrees.

        let eccentricity: f64 = (10_000.0 - 7_000.0) / (10_000.0 + 7_000.0);

        // roughly 2.1 hours
        let period: f64 =
            2.0 * std::f64::consts::PI / 398_600.0_f64.sqrt() * (17_000.0_f64 / 2.0).powf(1.5);

        let mean_anomaly_1 = 2.0 * std::f64::consts::PI / period * 0.5 * 60.0 * 60.0;

        let eccentric_anomaly_1 = newton(
            |x| x - eccentricity * x.sin() - mean_anomaly_1,
            |x| 1.0 - eccentricity * x.cos(),
            mean_anomaly_1 + eccentricity / 2.0,
            1.0e-6,
        )
        .unwrap()
        .root;

        let mean_anomaly_2 = 2.0 * std::f64::consts::PI / period * 1.5 * 60.0 * 60.0;

        // second point is beyond apogee, i.e. mean_anomaly_2 > pi
        let eccentric_anomaly_2 = newton(
            |x| x - eccentricity * x.sin() - mean_anomaly_2,
            |x| 1.0 - eccentricity * x.cos(),
            mean_anomaly_2 - eccentricity / 2.0,
            1.0e-6,
        )
        .unwrap()
        .root;

        let theta_1 = 2.0
            * (((1.0 + eccentricity) / (1.0 - eccentricity)).sqrt()
                * (eccentric_anomaly_1 / 2.0).tan())
            .atan();
        println!("{}", theta_1);

        // theta_2 is negative so subtract it from 2pi
        let theta_2 = 2.0 * std::f64::consts::PI
            + 2.0
                * (((1.0 + eccentricity) / (1.0 - eccentricity)).sqrt()
                    * (eccentric_anomaly_2 / 2.0).tan())
                .atan();
        println!("{}", theta_2);

        assert_approx_eq(
            theta_2 - theta_1,
            128.7 * 2.0 * std::f64::consts::PI / 360.0,
            1e-1,
            "Problem 3.9 solution",
        );
    }
}
