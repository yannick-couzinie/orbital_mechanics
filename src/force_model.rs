//! The force model for propagation with varying degrees of precision.

use nalgebra::SVector;
use physical_constants::NEWTONIAN_CONSTANT_OF_GRAVITATION;

pub trait ForceModel {
    // The ODE for this force-model
    fn ode(&self) -> impl Fn(f64, &SVector<f64, 6>) -> SVector<f64, 6>;
}

pub struct TwoBody {
    pub m1: f64,
    pub m2: f64,
}

impl ForceModel for TwoBody {
    fn ode(&self) -> impl Fn(f64, &SVector<f64, 6>) -> SVector<f64, 6> {
        let mu = NEWTONIAN_CONSTANT_OF_GRAVITATION * (self.m1 + self.m2);
        move |_t, y| {
            nalgebra::SVector::<f64, 6>::new(
                y[3],
                y[4],
                y[5],
                -mu * y[0] / (y[0].powi(2) + y[1].powi(2) + y[2].powi(2)).powf(3. / 2.),
                -mu * y[1] / (y[0].powi(2) + y[1].powi(2) + y[2].powi(2)).powf(3. / 2.),
                -mu * y[2] / (y[0].powi(2) + y[1].powi(2) + y[2].powi(2)).powf(3. / 2.),
            )
        }
    }
}
