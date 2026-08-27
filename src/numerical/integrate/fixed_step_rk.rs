//! One-dimensional Runge-Kutta numerical integration with methods from Rk1 to Rk4.

use super::errors::IntegrationError;
use super::structs::IntegrationResult;
use nalgebra;
use strum::EnumIter;

#[derive(EnumIter, Debug)]
pub enum FixedStepRkMethods {
    /// First-order Runge-Kutta numerical integration
    Rk1,
    /// Second-order Runge-Kutta numerical integration
    Rk2,
    /// Third-order Runge-Kutta numerical integration
    Rk3,
    /// Fourth-order Runge-Kutta numerical integration
    Rk4,
}

struct ButcherTableau {
    a: &'static [f64],
    b: &'static [&'static [f64]],
    c: &'static [f64],
}

static RK1_TABLEAU: ButcherTableau = ButcherTableau {
    a: &[0.],
    b: &[&[0.]],
    c: &[1.],
};
static RK2_TABLEAU: ButcherTableau = ButcherTableau {
    a: &[0., 1.],
    b: &[&[0.], &[1.0]],
    c: &[0.5, 0.5],
};
static RK3_TABLEAU: ButcherTableau = ButcherTableau {
    a: &[0., 0.5, 1.],
    b: &[&[0., 0.], &[0., 0.5], &[-1., 2.]],
    c: &[1. / 6., 2. / 3., 1. / 6.],
};
static RK4_TABLEAU: ButcherTableau = ButcherTableau {
    a: &[0., 0.5, 0.5, 1.],
    b: &[&[0., 0., 0.], &[0.5, 0., 0.], &[0., 0.5, 0.], &[0., 0., 1.]],
    c: &[1. / 6., 1. / 3., 1. / 3., 1. / 6.],
};

impl FixedStepRkMethods {
    fn tableau(&self) -> &'static ButcherTableau {
        match self {
            Self::Rk1 => &RK1_TABLEAU,
            Self::Rk2 => &RK2_TABLEAU,
            Self::Rk3 => &RK3_TABLEAU,
            Self::Rk4 => &RK4_TABLEAU,
        }
    }
}

pub fn rk1_4<const N: usize, F>(
    // The function to integrate, the input and output vector sizes (y) need to be the same, the
    // input is t, y
    ode_function: F,
    // The timespan with the first entry being the start and the second the endpoint
    tspan: (f64, f64),
    // Column vector of the initial values of the vector y
    y0: nalgebra::SVector<f64, N>,
    // Time step
    h: f64,
    rk: FixedStepRkMethods,
) -> Result<IntegrationResult<N>, IntegrationError<N>>
where
    F: Fn(f64, &nalgebra::SVector<f64, N>) -> nalgebra::SVector<f64, N>,
{
    let t0 = tspan.0;
    let tf = tspan.1;
    let mut step_h = h;
    let mut t = t0;
    let mut y = y0;
    let mut tout = vec![t];
    let mut yout = vec![y];
    let tableau = rk.tableau();

    // only forward propagation
    if !h.is_finite() || h <= 0.0 {
        return Err(IntegrationError::InvalidStepSize { step: h });
    }

    // only forward propagation and start has to be before end
    if !t0.is_finite() || !tf.is_finite() || t0 > tf {
        return Err(IntegrationError::InvalidTimeSpan { start: t0, end: tf });
    }

    if !y0.iter().all(|x| x.is_finite()) {
        return Err(IntegrationError::NonFiniteState { state: y0 });
    }

    while t < tf {
        let ti = t;
        let mut f: Vec<nalgebra::SVector<f64, N>> = Vec::new(); // the [n_stages, dim(y)] sized matrix with the function values
        step_h = step_h.min(tf - t);

        // make sure that t is not that big that the step gets eaten up in floating point precision
        if t + step_h == t {
            return Err(IntegrationError::StepDoesNotAdvanceTime {
                time: t,
                step: step_h,
            });
        }

        // evaluate the time derivates at the 'n_stages' points within the current interval
        for i in 0..tableau.a.len() {
            let t_inner = ti + tableau.a.get(i).unwrap() * step_h;
            let mut y_inner = y;

            for (j, f_entry) in f.iter().enumerate() {
                // this will not run on i=0
                y_inner.axpy(step_h * tableau.b[i][j], f_entry, 1.0);
            }

            let f_eval = ode_function(t_inner, &y_inner);

            if !f_eval.iter().all(|x| x.is_finite()) {
                return Err(IntegrationError::NonFiniteDerivative {
                    derivative: f_eval,
                    state: y_inner,
                    time: t,
                });
            }

            f.push(f_eval);
        }
        t += step_h;

        for (i, f_entry) in f.iter().enumerate() {
            y.axpy(step_h * tableau.c[i], f_entry, 1.0)
        }
        tout.push(t);
        yout.push(y);
    }
    Ok(IntegrationResult {
        times: tout,
        states: yout,
    })
}

#[cfg(test)]
mod tests {
    // we test some of the private components in this module (i.e. the FixedStepRkMethods)
    use super::*;
    use crate::test_utils::assert_approx_eq;
    use strum::IntoEnumIterator;

    #[test]
    fn test_coefficient_dimensions() {
        for rk_type in FixedStepRkMethods::iter() {
            assert_eq!(rk_type.tableau().a.len(), rk_type.tableau().b.len());
            assert_eq!(rk_type.tableau().a.len(), rk_type.tableau().c.len());

            for row in rk_type.tableau().b {
                assert!(
                    row.len() >= rk_type.tableau().a.len().saturating_sub(1),
                    "insufficient b coefficients for {rk_type:?}"
                )
            }
        }
    }

    #[test]
    fn test_c_coefficients_sum() {
        // see page 40 from Curtis on the top, each sum of c needs tob e 1 and each row of b equals
        // the value of that row in a.
        for rk_type in FixedStepRkMethods::iter() {
            let c_coefficient_sum: f64 = rk_type.tableau().c.into_iter().sum();
            assert_approx_eq(
                c_coefficient_sum,
                1.0,
                1.0e-6,
                format!("c coefficients do not sum to one for {rk_type:?}"),
            );
        }
    }

    #[test]
    fn test_ab_coefficients_sum() {
        // see page 40 from Curtis on the top, each sum of c needs tob e 1 and each row of b equals
        // the value of that row in a.
        for rk_type in FixedStepRkMethods::iter() {
            for (i, a_row) in rk_type.tableau().a.into_iter().enumerate() {
                let b_row_coefficient_sum: f64 = rk_type.tableau().b[i].iter().sum();
                assert_approx_eq(
                    b_row_coefficient_sum,
                    *a_row,
                    1.0e-6,
                    format!("b row {i} does not sum to a[{i}] for {rk_type:?}"),
                );
            }
        }
    }

    #[test]
    fn check_non_divisible_step_length() {
        // check that with a step length that is not divisible by the time range we correctly adapt
        // the last step
        let integration_result = rk1_4(
            |_t, y| nalgebra::SVector::<f64, 1>::new(y[0]),
            (0.0, 0.1),
            nalgebra::SVector::<f64, 1>::new(1.),
            0.3,
            FixedStepRkMethods::Rk2,
        )
        .expect("Rk run failed in test");

        assert_approx_eq(
            *integration_result.times.last().unwrap(),
            0.1,
            1e-6,
            "Step length 0.3 does not adapt correctly to time range (0.0, 0.1).",
        );

        assert_approx_eq(
            integration_result.states.last().unwrap()[0],
            1.105, // roughly exp(0.1)
            1e-3,
            "State is wrongly calculated with non-divisible step length",
        );
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
        let integration_result = rk1_4(
            |_t, y| nalgebra::SVector::<f64, 4>::new(y[1], y[2], y[3], -y[0] - 2.0 * y[2]),
            (0.0, 20.0),
            nalgebra::SVector::<f64, 4>::new(1., 0., 0., 0.),
            0.01,
            FixedStepRkMethods::Rk4,
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
}
