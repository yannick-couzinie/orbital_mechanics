//! The central struct of this library that contains information about the spacecraft.

use crate::state;

pub struct Spacecraft {
    pub state: state::State,
    pub mass: f64,
    pub time: f64,
}

impl Spacecraft {
    pub fn new(state: impl Into<state::State>, mass: f64, time: f64) -> Self {
        Self {
            state: state.into(),
            mass,
            time,
        }
    }
}
