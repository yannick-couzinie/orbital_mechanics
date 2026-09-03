//! Runge-Kutta-Fehlberg adaptive step integration.

use super::Integrator;
use super::errors::IntegrationError;
use super::rk_utils::{calculate_stages, combine_stages};
use super::structs::IntegrationResult;
use nalgebra;

pub struct AdaptiveRkParameters {
    pub step: f64,
    pub method: AdaptiveStepRkMethod,
    pub tolerance: f64,
}

#[derive(Debug, Clone, Copy)]
pub enum AdaptiveStepRkMethod {
    Rkf45,
}

struct ButcherTableau {
    a: &'static [f64],
    b: &'static [&'static [f64]],
    c1: &'static [f64],
    c2: &'static [f64],
}

static RK45_TABLEAU: ButcherTableau = ButcherTableau {
    a: &[0., 1. / 4., 3. / 8., 12. / 13., 1., 1. / 2.],
    b: &[
        &[0., 0., 0., 0., 0.],
        &[1. / 4., 0., 0., 0., 0.],
        &[3. / 32., 9. / 32., 0., 0., 0.],
        &[1932. / 2197., -7200. / 2197., 7296. / 2197., 0., 0.],
        &[439. / 216., -8., 3680. / 513., -845. / 4104., 0.],
        &[-8. / 27., 2., -3544. / 2565., 1859. / 4104., -11. / 40.],
    ],
    c1: &[25. / 216., 0., 1408. / 2565., 2197. / 4104., -1. / 5., 0.],
    c2: &[
        16. / 135.,
        0.,
        6656. / 12825.,
        28561. / 56430.,
        -9. / 50.,
        2. / 55.,
    ],
};

impl AdaptiveStepRkMethod {
    fn tableau(&self) -> &'static ButcherTableau {
        match self {
            Self::Rkf45 => &RK45_TABLEAU,
        }
    }
    fn power(&self) -> f64 {
        match self {
            Self::Rkf45 => 4.,
        }
    }
}

impl Integrator for AdaptiveRkParameters {
    fn integrate<const N: usize, F>(
        &self,
        // The function to integrate, the input and output vector sizes (y) need to be the same, the
        // input is t, y
        ode_function: F,
        // The timespan with the first entry being the start and the second the endpoint
        tspan: (f64, f64),
        // Column vector of the initial values of the vector y
        y0: nalgebra::SVector<f64, N>,
        // Time step
    ) -> Result<IntegrationResult<N>, IntegrationError<N>>
    where
        F: Fn(f64, &nalgebra::SVector<f64, N>) -> nalgebra::SVector<f64, N>,
    {
        let t0 = tspan.0;
        let tf = tspan.1;
        let mut t = t0;
        let mut y = y0;
        let mut tout = vec![t];
        let mut yout = vec![y];
        let tableau = self.method.tableau();
        let mut stages: Vec<nalgebra::SVector<f64, N>> = Vec::with_capacity(tableau.a.len());

        // only forward propagation
        if !self.step.is_finite() || self.step <= 0.0 {
            return Err(IntegrationError::InvalidStepSize { step: self.step });
        }

        // only forward propagation and start has to be before end
        if !t0.is_finite() || !tf.is_finite() || t0 > tf {
            return Err(IntegrationError::InvalidTimeSpan { start: t0, end: tf });
        }

        if !y0.iter().all(|x| x.is_finite()) {
            return Err(IntegrationError::NonFiniteState { state: y0 });
        }

        if !self.tolerance.is_finite() || self.tolerance == 0. || self.tolerance.is_sign_negative()
        {
            return Err(IntegrationError::AdaptiveParametersIncomplete {});
        }

        let mut step_h = self.step;

        while t < tf {
            let ti = t;
            step_h = step_h.min(tf - t);
            let mut error: f64;
            // make sure that t is not that big that the step gets eaten up in floating point precision
            if t + step_h == t {
                return Err(IntegrationError::StepDoesNotAdvanceTime {
                    time: t,
                    step: step_h,
                });
            }

            loop {
                calculate_stages(
                    &mut stages,
                    y,
                    ti,
                    step_h,
                    tableau.a,
                    tableau.b,
                    &ode_function,
                )?;

                let y_low_dim = combine_stages(&stages, y, step_h, tableau.c1);
                let y_high_dim = combine_stages(&stages, y, step_h, tableau.c2);

                error = (y_high_dim - y_low_dim).abs().max();

                if !error.is_finite() {
                    return Err(IntegrationError::NonFiniteError { error });
                }

                // when error = 0 use 5.0 as a replacement for infinity
                let factor = if error == 0.0 {
                    5.0
                } else {
                    0.8 * (self.tolerance / error).powf(1.0 / (self.method.power() + 1.0))
                };

                if error < self.tolerance {
                    t += step_h;
                    tout.push(t);
                    y = y_high_dim;
                    // this increases the step if error < tolerance
                    step_h *= factor;
                    break;
                } else {
                    // if you do not break update the step which effectively reduces it
                    step_h *= factor;
                }
            }

            if !y.iter().all(|x| x.is_finite()) {
                return Err(IntegrationError::NonFiniteState { state: y });
            }

            yout.push(y);
        }
        Ok(IntegrationResult {
            times: tout,
            states: yout,
        })
    }
}

#[cfg(test)]
mod tests {
    // we test some of the private components in this module
    use super::*;
    use crate::test_utils::assert_approx_eq;

    #[test]
    fn test_coefficient_dimensions() {
        assert_eq!(RK45_TABLEAU.a.len(), RK45_TABLEAU.b.len());
        assert_eq!(RK45_TABLEAU.a.len(), RK45_TABLEAU.c1.len());
        assert_eq!(RK45_TABLEAU.a.len(), RK45_TABLEAU.c2.len());

        for row in RK45_TABLEAU.b {
            assert!(
                row.len() >= RK45_TABLEAU.a.len().saturating_sub(1),
                "insufficient b coefficients for RK45"
            )
        }
    }

    #[test]
    fn test_c_coefficients_sum() {
        // see page 40 from Curtis on the top, each sum of c needs to be 1 and each row of b equals
        // the value of that row in a.
        let c_coefficient_sum: f64 = RK45_TABLEAU.c1.iter().sum();
        assert_approx_eq(
            c_coefficient_sum,
            1.0,
            1.0e-6,
            "c1 coefficients do not sum to one for RK45",
        );

        let c_coefficient_sum: f64 = RK45_TABLEAU.c2.iter().sum();
        assert_approx_eq(
            c_coefficient_sum,
            1.0,
            1.0e-6,
            "c2 coefficients do not sum to one for RK45",
        );
    }

    #[test]
    fn test_ab_coefficients_sum() {
        // see page 40 from Curtis on the top, each sum of c needs to be 1 and each row of b equals
        // the value of that row in a.
        for (i, a_row) in RK45_TABLEAU.a.iter().enumerate() {
            let b_row_coefficient_sum: f64 = RK45_TABLEAU.b[i].iter().sum();
            assert_approx_eq(
                b_row_coefficient_sum,
                *a_row,
                1.0e-6,
                format!("b row {i} does not sum to a[{i}] for RK45"),
            );
        }
    }

    #[test]
    fn test_large_step() {
        // exponential equation solved with large step will reduce the step first
        let integrator = AdaptiveRkParameters {
            step: 10.,
            method: AdaptiveStepRkMethod::Rkf45,
            tolerance: 1e-12,
        };
        let integration_result = integrator
            .integrate(
                |_t, y| nalgebra::SVector::<f64, 1>::new(y[0]),
                (0.0, 12.0),
                nalgebra::SVector::<f64, 1>::new(1.),
            )
            .expect("Rk run failed in test");

        // assert the step size was reduced
        assert!(integration_result.times[1] - integration_result.times[0] < 10.);

        for (i, t_entry) in integration_result.times.into_iter().enumerate() {
            // get 0 since the first entry in the vector is equal to y
            let expected = t_entry.exp();
            let actual = integration_result.states[i][0];

            assert_approx_eq(actual, expected, 1.0e-6, "Non autonomous ODE solve failed");
        }
    }

    #[test]
    fn test_non_autonomous() {
        // the solution to y' = t with y(0) = 0 is y(t) = 1/2 t^2
        let integrator = AdaptiveRkParameters {
            step: 0.1,
            method: AdaptiveStepRkMethod::Rkf45,
            tolerance: 1e-10,
        };
        let integration_result = integrator
            .integrate(
                |t, _y| nalgebra::SVector::<f64, 1>::new(t),
                (0.0, 2.0),
                nalgebra::SVector::<f64, 1>::new(0.),
            )
            .expect("Rk run failed in test");

        for (i, t_entry) in integration_result.times.into_iter().enumerate() {
            // get 0 since the first entry in the vector is equal to y
            let expected = 0.5 * t_entry.powi(2);
            let actual = integration_result.states[i][0];

            assert_approx_eq(actual, expected, 1.0e-6, "Non autonomous ODE solve failed");
        }
    }

    #[test]
    fn problem118_rksolvers() {
        // Run the complete solver and compare with the analytical result.
        // this is the first order ssytem for d4y/dt4 + 2d2y/dt2 + y = 0 with initial conditions y=1 and
        // dy/dt = d2y/dt2 = d3y/dt3 = 0 at t=0 solved for y(20) which should result in 9.545
        // the analytical solution is actually t.cos() + (t / 2.0) * t.sin() so we can compare against
        // that
        //
        // run the tests only for Rk4 as that should be proof enough that the solvers work, and getting
        // the same accuracy on Rk1 would require thousand-fold of steps
        let integrator = AdaptiveRkParameters {
            step: 0.01,
            tolerance: 1e-10,
            method: AdaptiveStepRkMethod::Rkf45,
        };
        let integration_result = integrator
            .integrate(
                |_t, y| nalgebra::SVector::<f64, 4>::new(y[1], y[2], y[3], -y[0] - 2.0 * y[2]),
                (0.0, 20.0),
                nalgebra::SVector::<f64, 4>::new(1., 0., 0., 0.),
            )
            .expect("Rk run failed in test");

        for (i, t_entry) in integration_result.times.into_iter().enumerate() {
            // get 0 since the first entry in the vector is equal to y
            let expected = t_entry.cos() + (t_entry * 0.5 * t_entry.sin());
            let actual = integration_result.states[i][0];

            assert_approx_eq(
                actual,
                expected,
                1.0e-6,
                format!("problem 1.18 failed using Rk4 at step {i}, t={t_entry}"),
            );
        }
    }

    #[test]
    fn harmonic_oscillator() {
        // The harmonic oscillator is d2x/dt2 = -x
        let integrator = AdaptiveRkParameters {
            step: 0.01,
            tolerance: 1e-10,
            method: AdaptiveStepRkMethod::Rkf45,
        };
        let integration_result = integrator
            .integrate(
                |_t, y| nalgebra::SVector::<f64, 2>::new(y[1], -y[0]),
                (0.0, 20.0),
                nalgebra::SVector::<f64, 2>::new(1., 0.),
            )
            .expect("Rk run failed in test");

        for (i, t_entry) in integration_result.times.into_iter().enumerate() {
            // we get cos since the initial condition is 1, i.e. the offset from sin is pi/2
            let expected = t_entry.cos();
            let actual = integration_result.states[i][0];

            assert_approx_eq(
                actual,
                expected,
                1.0e-6,
                "Could not solve harmonic oscillator correctly.",
            );
        }
    }
}
