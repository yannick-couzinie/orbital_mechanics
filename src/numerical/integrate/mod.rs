pub mod adaptive_step_rk;
pub mod errors;
pub mod fixed_step_rk;
pub(crate) mod rk_utils;
pub mod structs;

use nalgebra::SVector;

pub trait Integrator {
    fn integrate<const N: usize, F>(
        &self,
        ode_function: F,
        tspan: (f64, f64),
        y0: SVector<f64, N>,
    ) -> Result<structs::IntegrationResult<N>, errors::IntegrationError<N>>
    where
        F: Fn(f64, &SVector<f64, N>) -> SVector<f64, N>;
}
