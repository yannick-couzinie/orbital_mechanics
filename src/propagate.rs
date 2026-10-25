//! Propagation code

use crate::force_model::ForceModel;
use crate::numerical::integrate::{Integrator, errors::IntegrationError};
use crate::numerical::root::{errors::RootError, newton};
use crate::numerical::stumpff;
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

fn universal_kepler_equation_solver(
    delta_time: f64,
    r0: f64,
    vr0: f64,
    alpha: f64,
    mu: f64,
) -> Result<f64, RootError> {
    // Algorithm 3.3 from Curtis, i.e. solve the universal Kepler's equation for the universal
    // anomaly x given delta t, r0, vr0 and alpha.
    let init_estimate = mu.sqrt() * alpha.abs() * delta_time;

    // this is inefficient because we calculate C2 and C3 separately for botht he function eval and
    // derivative eval but it'll do for now
    let newton_result = newton::newton(
        |x| {
            r0 * vr0 / mu.sqrt() * x.powi(2) * stumpff::c2(alpha * x.powi(2))
                + (1.0 - alpha * r0) * x.powi(3) * stumpff::c3(alpha * x.powi(2))
                + r0 * x
                - mu.sqrt() * delta_time
        },
        |x| {
            r0 * vr0 / mu.sqrt() * x * (1.0 - alpha * x.powi(2) * stumpff::c3(alpha * x.powi(2)))
                + (1.0 - alpha * r0) * x.powi(2) * stumpff::c2(alpha * x.powi(2))
                + r0
        },
        init_estimate,
        1.0e-6,
    )?;

    Ok(newton_result.root)
}

pub fn propagate_universal_variable(
    initial: &State,
    delta_time: f64,
    mu: f64,
) -> Result<State, RootError> {
    // Algorithm 3.4 from Curtis, i.e. given r0 and v0 initial state vector find r, v at a time delta t later.

    let init_radius_km = initial.radius();
    let init_speed_kms = initial.speed();
    let vr0 = initial.position().dot(&initial.velocity()) / init_radius_km;
    let alpha = 2.0 / init_radius_km - init_speed_kms.powi(2) / mu;

    let universal_anomaly =
        universal_kepler_equation_solver(delta_time, init_radius_km, vr0, alpha, mu)?;

    // find the lagrange parameters
    let stumpff_arg = alpha * universal_anomaly.powi(2);
    let stumpff_s = stumpff::c3(stumpff_arg);
    let stumpff_c = stumpff::c2(stumpff_arg);
    let f = 1.0 - universal_anomaly.powi(2) / init_radius_km * stumpff_c;
    let g = delta_time - 1.0 / mu.sqrt() * universal_anomaly.powi(3) * stumpff_s;
    let new_position = f * initial.position() + g * initial.velocity();
    let new_radius_km = new_position.norm();
    let fdot = mu.sqrt() / (new_radius_km * init_radius_km)
        * (alpha * universal_anomaly.powi(3) * stumpff_s - universal_anomaly);
    let gdot = 1.0 - universal_anomaly.powi(2) / new_radius_km * stumpff_c;
    let new_velocity = fdot * initial.position() + gdot * initial.velocity();

    Ok(State::new(new_position, new_velocity))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_utils::{assert_approx_eq, assert_state_approx_eq};

    #[test]
    fn test_problem_319() {
        //! At a given instant, the radial position of an earth-orbiting satellite is 7200km and its
        //! radial speed is 1km/s. If the semimajor axis is 10000km, use Algorithm 3.3 to find the
        //! universal anomaly 60 min later. Check your result using Eq. (3.58).
        let time_s = 60.0 * 60.0;

        let calculated =
            universal_kepler_equation_solver(time_s, 7200.0, 1.0, 1.0 / 10_000.0, 398600.0)
                .unwrap();

        // Independently calculate the initial eccentric anomaly.
        let r0_km: f64 = 7200.0;
        let vr0_kms: f64 = 1.0;
        let sma_km: f64 = 10_000.0;
        let mu_km3s2: f64 = 398_600.0;

        // From:
        // r = a(1 - e cos(E))
        // vr = sqrt(mu * a) / r * e sin(E)
        let e_cos_e0: f64 = 1.0 - r0_km / sma_km;
        let e_sin_e0: f64 = r0_km * vr0_kms / (mu_km3s2 * sma_km).sqrt();

        let eccentricity = e_cos_e0.hypot(e_sin_e0);
        let eccentric_anomaly_0 = e_sin_e0.atan2(e_cos_e0);

        let mean_anomaly_0 = eccentric_anomaly_0 - eccentricity * eccentric_anomaly_0.sin();

        let mean_motion = (mu_km3s2 / sma_km.powi(3)).sqrt();
        let mean_anomaly_1 = mean_anomaly_0 + mean_motion * time_s;

        let eccentric_anomaly_1 = newton::newton(
            |x| x - eccentricity * x.sin() - mean_anomaly_1,
            |x| 1.0 - eccentricity * x.cos(),
            mean_anomaly_1,
            1.0e-10,
        )
        .unwrap()
        .root;

        // Equation 3.58: chi = sqrt(a) * change in eccentric anomaly.
        let expected = sma_km.sqrt() * (eccentric_anomaly_1 - eccentric_anomaly_0);

        assert_approx_eq(calculated, expected, 1e-6, "Problem 3.19 solution");
    }

    #[test]
    fn test_problem_320() {
        //! At a given instant, a space object has the following posiiton and velocity vectors
        //! relative to an earth-centered inertial frame of reference:
        //! r0 = 20_000i - 105_000j-19_000k (km)
        //! v0 = 0.9000i - 3.40000j - 1.5000k (km)
        //! Use Algorithm 3.4 to find r and v 2 h later.
        //!
        //! { Ans.:
        //!         r = 26338i - 128750j - 29656k (km)
        //!         v = 0.86280i - 3.2116j - 1.4613k (km/s)
        //! }

        let initial_state = State::from([20_000.0, -105_000.0, -19_000.0, 0.9, -3.4, -1.5]);

        let calculated =
            propagate_universal_variable(&initial_state, 2.0 * 60.0 * 60.0, 398600.0).unwrap();

        let expected = State::from([26_338.0, -128_750.0, -29_656.0, 0.8628, -3.2116, -1.4613]);

        assert_state_approx_eq(&expected, &calculated, 5.1, 1e-4, "Problem 3.20");
    }
}
