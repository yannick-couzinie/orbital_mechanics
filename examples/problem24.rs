//! At a given instant t0, a 1000-kg earth-orbiting satellite has the inertial position and velocity
//! vectors r0 = 3207i+5459j+2714k (km) and v0 = -6.532i 0.7835j + 6.142k (km/s). Solve Eq. (2.22)
//! [the two-body problem] numerically to find the maximum altitude reachedby the satellite and the
//! time at which it occurs.
//!
//! {Ans.: Using MATLAB's ode45, the maximum altitude 9670km at 1.66 h after t0}

const MASS_EARTH: f64 = 5.9722e24;

use curtis_orbital_mechanics::numerical::integrate::adaptive_step_rk::{
    AdaptiveRkParameters, AdaptiveStepRkMethod,
};
use curtis_orbital_mechanics::state::{ReferenceFrame, State};
use curtis_orbital_mechanics::{force_model::TwoBody, propagate, spacecraft};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let sc = spacecraft::Spacecraft::new(
        State::from_array(
            [3207000., 5459000., 2714000., -6532., 783.5, 6142.],
            ReferenceFrame::GeocentricEquatorial,
        ),
        1000.,
        0.0,
    );

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
    let traj = propagate::propagate(sc, &integrator, &force_model, 4. * 60. * 60.)?;

    let mut max_val = 0.0;
    let mut time_of_max = 0.0;

    for (time, state) in std::iter::zip(traj.times, traj.states) {
        let norm = (state.state_vector[0].powi(2)
            + state.state_vector[1].powi(2)
            + state.state_vector[2].powi(2))
        .sqrt();
        if norm > max_val {
            max_val = dbg!(norm);
            time_of_max = dbg!(time);
        }
    }

    println!("Max achieved after {} hours!", time_of_max / 3600.);

    Ok(())
}
