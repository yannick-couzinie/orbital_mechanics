//! One-dimensional Runge-Kutta numerical integration with methods from RK1 to RK4.

use super::errors::IntegrationError;
use strum::EnumIter;

/// The RkTypes we can use

#[derive(EnumIter, Debug)]
pub enum RkTypes {
    /// First-order Runge-Kutta numerical integration
    RK1,
    /// Second-order Runge-Kutta numerical integration
    RK2,
    /// Third-order Runge-Kutta numerical integration
    RK3,
    /// Fourth-order Runge-Kutta numerical integration
    RK4,
}

impl RkTypes {
    fn n_stages(&self) -> usize {
        match self {
            RkTypes::RK1 => 1,
            RkTypes::RK2 => 2,
            RkTypes::RK3 => 3,
            RkTypes::RK4 => 4,
        }
    }

    fn a(&self) -> Vec<f64> {
        match self {
            RkTypes::RK1 => vec![0.],
            RkTypes::RK2 => vec![0., 1.],
            RkTypes::RK3 => vec![0., 0.5, 1.],
            RkTypes::RK4 => vec![0., 0.5, 0.5, 1.],
        }
    }

    fn b(&self) -> Vec<Vec<f64>> {
        match self {
            RkTypes::RK1 => vec![vec![0.]],
            RkTypes::RK2 => vec![vec![0.], vec![1.]],
            RkTypes::RK3 => vec![vec![0., 0.], vec![0., 0.5], vec![-1., 2.]],
            RkTypes::RK4 => {
                vec![
                    vec![0., 0., 0.],
                    vec![0.5, 0., 0.],
                    vec![0., 0.5, 0.],
                    vec![0., 0., 1.],
                ]
            }
        }
    }

    fn c(&self) -> Vec<f64> {
        match self {
            RkTypes::RK1 => vec![1.],
            RkTypes::RK2 => vec![0.5, 0.5],
            RkTypes::RK3 => vec![1. / 6., 2. / 3., 1. / 6.],
            RkTypes::RK4 => vec![1. / 6., 1. / 3., 1. / 3., 1. / 6.],
        }
    }
}

pub fn rk1_4(
    // The function to integrate, the input and output vector sizes (y) need to be the same, the
    // input is t, y
    ode_function: &dyn Fn(f64, &Vec<f64>) -> Vec<f64>,
    // The timespan with the first entry being the start and the second the endpoint
    tspan: (f64, f64),
    // Column vector of the initial values of the vector y
    y0: Vec<f64>,
    // Time step
    h: f64,
    rk: RkTypes,
) -> Result<(Vec<f64>, Vec<Vec<f64>>), IntegrationError> {
    let t0 = tspan.0;
    let tf = tspan.1;
    let mut step_h = h;
    let mut t = t0;
    let mut y = y0.clone();
    let mut tout = vec![t];
    let mut yout = vec![y.clone()];

    // only forward propagation
    if !h.is_finite() || h <= 0.0 {
        return Err(IntegrationError::InvalidStepSize { step: h });
    }

    // only forward propagation and start has to be before end
    if !t0.is_finite() || !tf.is_finite() || t0 > tf {
        return Err(IntegrationError::InvalidTimeSpan { start: t0, end: tf });
    }

    // make sure that t is not that big that the step gets eaten up in floating point precision
    if t + step_h == t {
        return Err(IntegrationError::StepDoesNotAdvanceTime {
            time: t,
            step: step_h,
        });
    }

    while t < tf {
        let ti = t;
        let yi = y.clone();
        let mut f: Vec<Vec<f64>> = Vec::new(); // the [n_stages, dim(y)] sized matrix with the function values
        step_h = step_h.min(tf - t);

        // evaluate the time derivates at the 'n_stages' points within the current interval
        for i in 0..rk.n_stages() {
            let t_inner = ti + rk.a().get(i).unwrap() * step_h;
            let mut y_inner = yi.clone();
            if i == 0 {
                let f_eval = ode_function(t_inner, &y_inner);
                for f_entry in &f_eval {
                    // f is an empty vector at this point
                    // we want to push to f(:, i) where the : here is iterated in this for-loop
                    // and i = 0 (i.e. first entry)
                    f.push(vec![*f_entry])
                }
            } else {
                for j in 0..i {
                    // this will not run on i=1
                    for (k, y_entry) in y_inner.iter_mut().enumerate() {
                        *y_entry += step_h
                            * rk.b().get(i).unwrap().get(j).unwrap()
                            * f.get(k).unwrap().get(j).unwrap();
                    }
                }

                let f_eval = ode_function(t_inner, &y_inner);
                for (j, f_entry) in f_eval.into_iter().enumerate() {
                    // f is an empty vector at this point
                    // we want to push to f(:, i) where the : here is iterated in this for-loop
                    // and i = 0 (i.e. first entry)
                    f[j].push(f_entry)
                }
            }
        }
        t += step_h;
        let mut newy_vec = Vec::new();
        // loop over dimensions of y
        for (k1, y_entry) in yi.into_iter().enumerate() {
            // inner product  loop between f and c
            let mut new_y_entry = y_entry;
            for (k2, fentry) in f[k1].iter().enumerate() {
                new_y_entry += step_h * *fentry * rk.c()[k2];
            }
            newy_vec.push(new_y_entry)
        }
        y = newy_vec.clone();
        tout.push(t);
        yout.push(newy_vec);
    }
    Ok((tout, yout))
}

#[cfg(test)]
mod tests {
    // we test some of the private components in this module (i.e. the RkTypes)
    use super::*;
    use crate::test_utils::assert_approx_eq;
    use strum::IntoEnumIterator;

    #[test]
    fn test_coefficient_dimensions() {
        for rk_type in RkTypes::iter() {
            assert_eq!(rk_type.a().len(), rk_type.n_stages());
            assert_eq!(rk_type.b().len(), rk_type.n_stages());
            assert_eq!(rk_type.c().len(), rk_type.n_stages());
            for row in rk_type.b() {
                assert!(
                    row.len() >= rk_type.n_stages().saturating_sub(1),
                    "insufficient b coefficients for {rk_type:?}"
                )
            }
        }
    }

    #[test]
    fn test_c_coefficients_sum() {
        // see page 40 from Curtis on the top, each sum of c needs tob e 1 and each row of b equals
        // the value of that row in a.
        for rk_type in RkTypes::iter() {
            let c_coefficient_sum: f64 = rk_type.c().into_iter().sum();
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
        for rk_type in RkTypes::iter() {
            for (i, a_row) in rk_type.a().into_iter().enumerate() {
                let b_row_coefficient_sum: f64 = rk_type.b().get(i).unwrap().iter().sum();
                assert_approx_eq(
                    b_row_coefficient_sum,
                    a_row,
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
        let (times, states) = rk1_4(
            &|_t, y: &Vec<f64>| vec![y[0]],
            (0.0, 0.1),
            vec![1.],
            0.3,
            RkTypes::RK2,
        )
        .expect("Rk run failed in test");

        assert_approx_eq(
            *times.last().unwrap(),
            0.1,
            1e-6,
            "Step length 0.3 does not adapt correctly to time range (0.0, 0.1).",
        );

        assert_approx_eq(
            *states.last().unwrap().first().unwrap(),
            1.105, // roughly exp(0.1)
            1e-3,
            "State is wrongly calculated with non-divisible step length",
        );
    }
}
