//! Struct that contains a state

use nalgebra::SVector;
use nalgebra::Vector3;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReferenceFrame {
    GeocentricEquatorial,
}

pub struct State {
    // state_vector is usually given in x, y, z, vx, vy, vz
    pub state_vector: SVector<f64, 6>,
    pub frame: ReferenceFrame,
}

impl State {
    pub fn from_vector(state_vector: SVector<f64, 6>, frame: ReferenceFrame) -> Self {
        Self {
            state_vector,
            frame,
        }
    }

    pub fn from_array(values: [f64; 6], frame: ReferenceFrame) -> Self {
        Self {
            state_vector: SVector::from_row_slice(&values),
            frame,
        }
    }

    pub fn new(position: Vector3<f64>, velocity: Vector3<f64>, frame: ReferenceFrame) -> Self {
        Self::from_array(
            [
                position.x, position.y, position.z, velocity.x, velocity.y, velocity.z,
            ],
            frame,
        )
    }
    pub fn position(&self) -> Vector3<f64> {
        self.state_vector.fixed_rows::<3>(0).into_owned()
    }
    pub fn radius(&self) -> f64 {
        self.position().norm()
    }

    pub fn velocity(&self) -> Vector3<f64> {
        self.state_vector.fixed_rows::<3>(3).into_owned()
    }
    pub fn speed(&self) -> f64 {
        self.velocity().norm()
    }

    pub fn is_approx_eq(
        &self,
        other: &Self,
        position_tolerance: f64,
        velocity_tolerance: f64,
    ) -> bool {
        let position_error = (self.position() - other.position()).norm();
        let velocity_error = (self.velocity() - other.velocity()).norm();

        position_error <= position_tolerance && velocity_error <= velocity_tolerance
    }
}
