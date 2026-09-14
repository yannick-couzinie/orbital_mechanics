//! Propagation code

use crate::force_model::ForceModel;
use crate::numerical::integrate::{Integrator, errors::IntegrationError};
use crate::spacecraft::Spacecraft;
use crate::state::State;

pub struct Trajectory {
    pub times: Vec<f64>,
    pub states: Vec<State>,
}

// For now the algorithm takes a central body mass as we only implement the two-body problem, as we
// move on to more general algorithms this will be changed. The idea would then be to have a sort of
// force-model (trait? struct?) in which the user can specify which forces to include and this force
// model would output a closure usable by the integrator
pub fn propagate(
    sc: Spacecraft,
    integrator: &impl Integrator,
    force_model: &impl ForceModel,
    end_time: f64,
) -> Result<Trajectory, IntegrationError<6>> {
    let int_result = integrator.integrate(
        force_model.ode(sc.mass),
        (sc.time, end_time),
        sc.state.state_vector,
    )?;

    Ok(Trajectory {
        times: int_result.times,
        states: int_result.states.into_iter().map(State::from).collect(),
    })
}
