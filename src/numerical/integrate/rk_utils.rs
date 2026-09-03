use super::errors::IntegrationError;
use nalgebra::SVector;

// Calculate the stages for RK integration
pub(super) fn calculate_stages<const N: usize, F>(
    stages: &mut Vec<SVector<f64, N>>,
    y: SVector<f64, N>,
    t: f64,
    step: f64,
    a: &'static [f64],
    b: &'static [&'static [f64]],
    ode_function: &F,
) -> Result<(), IntegrationError<N>>
where
    F: Fn(f64, &SVector<f64, N>) -> SVector<f64, N>,
{
    stages.clear();

    // evaluate the time derivates at the 'n_stages' points within the current interval
    for (i, a_entry) in a.iter().enumerate() {
        let t_inner = t + a_entry * step;
        let mut y_inner = y;

        for (j, f_entry) in stages.iter().enumerate() {
            // this will not run on i=0
            y_inner.axpy(step * b[i][j], f_entry, 1.0);
        }

        let f_eval = ode_function(t_inner, &y_inner);

        if !f_eval.iter().all(|x| x.is_finite()) {
            return Err(IntegrationError::NonFiniteDerivative {
                derivative: f_eval,
                state: y_inner,
                time: t_inner,
            });
        }

        stages.push(f_eval);
    }
    Ok(())
}

pub(super) fn combine_stages<const N: usize>(
    stages: &[SVector<f64, N>],
    y: SVector<f64, N>,
    step: f64,
    c: &'static [f64],
) -> SVector<f64, N> {
    let mut y_out = y;
    for (i, f_entry) in stages.iter().enumerate() {
        y_out.axpy(step * c[i], f_entry, 1.0)
    }
    y_out
}
