//! At a given instant, a 1000-kg earth-orbiting satellite has the inertial position and velocity
//! vectors r0 = 6600i (km) and v0 = 12j (km/s). Solve Eq. (2.22) numerically to find the distance
//! of the spacecraft from the center of the earth and its speed 24h later.
//!
//! {Ans.: Using MATLAB's ode45, distance = 456,500 km, speed = 5 km/2}

const MASS_EARTH: f64 = 5.9722e24;

use curtis_orbital_mechanics::numerical::integrate::adaptive_step_rk::{
    AdaptiveRkParameters, AdaptiveStepRkMethod,
};
use curtis_orbital_mechanics::{force_model::TwoBody, propagate, spacecraft};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let sc = spacecraft::Spacecraft::new([6600000., 0., 0., 0., 12000., 0.], 1000., 0.0);

    let integrator = AdaptiveRkParameters {
        step: 1.0,
        absolute_tolerance: 1e-3,
        relative_tolerance: 1e-9,
        method: AdaptiveStepRkMethod::Rkf45,
    };

    let force_model = TwoBody {
        primary_mass: MASS_EARTH,
    };

    // propagate for 4 hours
    let traj = propagate::propagate(sc, &integrator, &force_model, 24. * 60. * 60.)?;

    let final_state = traj.states.last().unwrap().state_vector;
    let final_distance =
        (final_state[0].powi(2) + final_state[1].powi(2) + final_state[2].powi(2)).sqrt();
    let final_speed =
        (final_state[3].powi(2) + final_state[4].powi(2) + final_state[5].powi(2)).sqrt();

    println!(
        "Final distance {}km and speed {}km/s!",
        final_distance / 1000.,
        final_speed / 1000.
    );

    Ok(())
}
