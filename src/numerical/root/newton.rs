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

        // theta_2 is negative so subtract it from 2pi
        let theta_2 = 2.0 * std::f64::consts::PI
            + 2.0
                * (((1.0 + eccentricity) / (1.0 - eccentricity)).sqrt()
                    * (eccentric_anomaly_2 / 2.0).tan())
                .atan();

        assert_approx_eq(
            theta_2 - theta_1,
            128.7 * 2.0 * std::f64::consts::PI / 360.0,
            1e-3,
            "Problem 3.9 solution",
        );
    }

    #[test]
    fn test_problem_310a() {
        //! An earth-orbiting satellite has a period of 14 h and a perigee radius of 10.000km. At
        //! time t=10h after perigee passage, determine
        //! (a) The radial position. (42356km)
        //! (b) The speed. (2.303 km/s)
        //! (c) The radial component of the velocity. (-1.271 km/s)

        let sma: f64 = ((14.0 * 60.0 * 60.0 * 398600.0_f64.sqrt()) / (2.0 * std::f64::consts::PI))
            .powf(2.0 / 3.0);
        let eccentricity: f64 = 1.0 - 10_000.0 / sma;

        let mean_anomaly = 2.0 * std::f64::consts::PI * 10.0 / 14.0;

        let eccentric_anomaly = newton(
            |x| x - eccentricity * x.sin() - mean_anomaly,
            |x| 1.0 - eccentricity * x.cos(),
            mean_anomaly - eccentricity / 2.0,
            1.0e-6,
        )
        .unwrap()
        .root;

        assert_approx_eq(
            sma * (1.0 - eccentricity * eccentric_anomaly.cos()),
            42354.92107805533,
            1e-3,
            "Problem 3.10 solution",
        );
    }

    #[test]
    fn test_problem_311() {
        //! A satellite in earth orbit has perigee and apogee radii of rp = 7500km and ra = 16000km,
        //! respectively. Find its true anomaly 40 minute after passing the true anomaly of 80
        //! degrees.
        //! { Ans.: 174.7deg }

        // we have already found that the time of passing of 80 degrees is 1473.5788s and so the
        // time of interest is
        let time_s: f64 = 3873.578838;
        let sma_km: f64 = 11750.0;
        let eccentricity: f64 = (16000.0 - 7500.0) / (16000.0 + 7500.0);

        // roughly 2.1 hours
        let period: f64 = 2.0 * std::f64::consts::PI / 398_600.0_f64.sqrt() * sma_km.powf(1.5);

        let mean_anomaly = 2.0 * std::f64::consts::PI / period * time_s;

        let eccentric_anomaly = newton(
            |x| x - eccentricity * x.sin() - mean_anomaly,
            |x| 1.0 - eccentricity * x.cos(),
            mean_anomaly + eccentricity / 2.0,
            1.0e-6,
        )
        .unwrap()
        .root;

        let theta_1 = 2.0
            * (((1.0 + eccentricity) / (1.0 - eccentricity)).sqrt()
                * (eccentric_anomaly / 2.0).tan())
            .atan();

        assert_approx_eq(
            theta_1,
            142.03 * 2.0 * std::f64::consts::PI / 360.0,
            1e-4,
            "Problem 3.11 solution",
        );
    }

    #[test]
    fn test_problem_316() {
        //! A spacecraft on a hyperbolic trajectory around the earth has a perigee ardius of 6600km
        //! and a perigee speed of 1.2v_esc. How far is the spacecraft from the center of the earth
        //! 24 h after passing through perigee.
        //! { Ans.: 656_610km }

        let rp_km: f64 = 6_600.0;
        // Velocity vp at perigee is sqrt(2mu/r)*1.2 and angular momentum is rp*vp
        let angular_momentum_km2s: f64 = (2.0 * 398600.0 * rp_km).sqrt() * 1.2;
        let time_s: f64 = 24.0 * 60.0 * 60.0;
        let eccentricity: f64 = angular_momentum_km2s.powi(2) / (398600.0 * rp_km) - 1.0;

        let hyperbolic_mean_anomaly = 398600.0_f64.powi(2) / angular_momentum_km2s.powi(3)
            * (eccentricity.powi(2) - 1.0).powf(1.5)
            * time_s;

        let eccentric_anomaly = newton(
            |x| eccentricity * x.sinh() - x - hyperbolic_mean_anomaly,
            |x| eccentricity * x.cosh() - 1.0,
            hyperbolic_mean_anomaly,
            1.0e-6,
        )
        .unwrap()
        .root;

        let true_anomaly = 2.0
            * (((eccentricity + 1.0) / (eccentricity - 1.0)).sqrt()
                * (eccentric_anomaly / 2.0).tanh())
            .atan();

        let distance_at_24h: f64 = angular_momentum_km2s.powi(2) / 398600.0 * 1.0
            / (1.0 + eccentricity * true_anomaly.cos());

        assert_approx_eq(distance_at_24h, 656610.0, 1e-0, "Problem 3.16 solution");
    }

    #[test]
    fn test_problem_317() {
        //! A trajectory has a perigee velocity 1.1vesc and a perigee altitude of 200km. If at 10
        //! a.m. the satellite is traveling toward the earth with a speed of 8km/s, how far will it
        //! be from the earth's surface at 5 p.m. the same day?
        //! { Ans.: 136_250km }

        let rp_km: f64 = 200.0 + 6378.0;
        let vp_kms: f64 = 1.1 * (2.0 * 398600.0 / rp_km).sqrt();

        let angular_momentum_km2s: f64 = rp_km * vp_kms;
        let eccentricity: f64 = angular_momentum_km2s.powi(2) / (398600.0 * rp_km) - 1.0;

        // we need to find the true anomaly when the trajectory had a velocity of 8km/s from that we
        // can calculate the time to perigee and then we can use 7 hours - time to perigee to get
        // the time after perigee and then we calculate the true anomaly and distance

        // The energy is constant v^2/2 - mu/r so use the perigee data and velocity to get the
        // distance
        let distance_at_10am_km =
            -(398600.0 / (vp_kms.powi(2) / 2.0 - 64.0 / 2.0 - 398600.0 / rp_km));

        // from the orbit equation we can get the true anomaly at 10am
        // acos is symmetric so we need to tell it which branch to choose, in this case we know the
        // spacecraft is travelling towards us so this should be negative
        let true_anomaly_at_10am =
            -((angular_momentum_km2s.powi(2) / (398600.0 * distance_at_10am_km) - 1.0)
                / eccentricity)
                .acos();

        let eccentric_anomaly_at_10am = 2.0
            * (((eccentricity - 1.0) / (eccentricity + 1.0)).sqrt()
                * (true_anomaly_at_10am / 2.0).tan())
            .atanh();

        let mean_anomaly_at_10am =
            eccentricity * eccentric_anomaly_at_10am.sinh() - eccentric_anomaly_at_10am;
        let time_to_perigee_from_10am = angular_momentum_km2s.powi(3) / 398600.0_f64.powi(2)
            * (eccentricity.powi(2) - 1.0).powf(-1.5)
            * mean_anomaly_at_10am;

        let time_s: f64 = (17.0 - 10.0) * 60.0 * 60.0 + time_to_perigee_from_10am;

        let hyperbolic_mean_anomaly = 398600.0_f64.powi(2) / angular_momentum_km2s.powi(3)
            * (eccentricity.powi(2) - 1.0).powf(1.5)
            * time_s;

        let eccentric_anomaly = newton(
            |x| eccentricity * x.sinh() - x - hyperbolic_mean_anomaly,
            |x| eccentricity * x.cosh() - 1.0,
            hyperbolic_mean_anomaly,
            1.0e-6,
        )
        .unwrap()
        .root;

        let true_anomaly = 2.0
            * (((eccentricity + 1.0) / (eccentricity - 1.0)).sqrt()
                * (eccentric_anomaly / 2.0).tanh())
            .atan();

        let distance: f64 = angular_momentum_km2s.powi(2) / 398600.0 * 1.0
            / (1.0 + eccentricity * true_anomaly.cos());

        assert_approx_eq(distance, 136250.0 + 6378.0, 1e-0, "Problem 3.16 solution");
    }
}
