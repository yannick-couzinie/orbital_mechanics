//! Struct that contains a state

use nalgebra::SVector;

pub struct State {
    // state_vector is usually given in x, y, z, vx, vy, vz
    pub state_vector: SVector<f64, 6>,
}

impl From<SVector<f64, 6>> for State {
    fn from(state_vector: SVector<f64, 6>) -> Self {
        Self { state_vector }
    }
}

impl From<[f64; 6]> for State {
    fn from(values: [f64; 6]) -> Self {
        Self {
            state_vector: SVector::from_row_slice(&values),
        }
    }
}
