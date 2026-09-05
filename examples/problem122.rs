//! Use one of the numerical methods discussed in this section to solve Eq. (1.127) for the time
//! required for the moon to fall to the earth after it is somehow stopped in its orbit while the
//! earth remains fixed in space. This will require a trial-and-error procedure known formally as a
//! shooting method. It is not necessary for this problem to code the procedure. Simply guess a time
//! and let the solver compute the final radius. On the basis of the deviation of that result from
//! the earth's radius (6378km), revise your time estimate and rerun the problem to compute a new
//! final radius. Repeat this process in a logical fashion until your time estimate yields a final
//! radius that is accurate to at least three significant figures. Compare your answer with the
//! analytical solution
//!
//! t = sqrt(r_0 / 2g_0 R_E^2) (pi/4 r_0) sqrt(r(r_0 - r)) + r_0/2 sin^(-1) ((r_0 - 2r)/r_0)
//!
//! where t is the time, r_0 is the initial radius, r is the final radius (r<r_0), g_0 is the sea
//! level acceleration of earth's gravity and R_E is the radius of the earth.

const R_E: f64 = 6378000.;
const R_M: f64 = 1737000.;
const G_0: f64 = 9.807; // from Eq. 1.35
const R_EM: f64 = 384400000.; // distance earth - moon (i.e. SMA of moon orbit)

use curtis_orbital_mechanics::numerical::integrate::{
    Integrator,
    adaptive_step_rk::{AdaptiveRkParameters, AdaptiveStepRkMethod},
    errors,
};
use std::f64::consts::PI;

fn calculate_analytical_solution(r: f64) -> f64 {
    (R_EM / (2.0 * G_0 * R_E.powi(2))).sqrt()
        * (PI / 4.0 * R_EM
            + (r * (R_EM - r)).sqrt()
            + R_EM / 2.0 * ((R_EM - 2.0 * r) / R_EM).asin())
}

fn calculate_final_radius(
    time: f64,
    integrator: &impl Integrator,
) -> Result<f64, errors::IntegrationError<2>> {
    let integration_result = integrator.integrate(
        |_t, y| nalgebra::SVector::<f64, 2>::new(y[1], -G_0 * R_E.powi(2) / y[0].powi(2)),
        (0.0, time),
        // initial distance is R_EM and it is stoped in its tracks so velocity 0
        nalgebra::SVector::<f64, 2>::new(R_EM, 0.),
    )?;

    Ok(integration_result.states.last().unwrap()[0])
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let integrator = AdaptiveRkParameters {
        step: 20.0,
        tolerance: 1e-6,
        method: AdaptiveStepRkMethod::Rkf45,
    };
    // we use y1 = x, y2 = y, y3 = z

    let initial_time = R_EM / 100.; // assume that it will be quicker than 100m/s on average

    println!(
        "We target a radius of R_E+R_M = {:?}km",
        (R_E + R_M) / 1000.
    );

    let analytical_solution = calculate_analytical_solution(R_E + R_M);

    // the initial guess goes beyond the moon crashing into earth which is a singularity so we have
    // two checks to see whether the time is too high either we are below R_E + R_M or the
    // calculation fails in the first place. We have a failing time above so let's do a binary
    // search.
    let mut cur_lower = 0.0;
    let mut cur_higher = initial_time;

    let mut loop_counter = 0;
    loop {
        loop_counter += 1;
        if loop_counter == 1000 {
            break;
        }
        let cur_time = (cur_lower + cur_higher) / 2.;
        println!(
            "Current time guess is: {}s, distance from analytical: {}",
            cur_time,
            analytical_solution - cur_time
        );

        let current_distance_from_earth_center = calculate_final_radius(cur_time, &integrator);

        match current_distance_from_earth_center {
            Ok(value) => {
                if (value - R_E - R_M).abs() < 1e-3 {
                    break;
                }
                if value > R_E + R_M {
                    // time is too short
                    cur_lower = cur_time;
                } else {
                    // time is too long
                    cur_higher = cur_time;
                    continue;
                }
            }
            Err(_) => {
                // time is too long
                cur_higher = cur_time;
            }
        }
    }
    println!(
        "Found solution after {} iterations as {}s, distance from analytical is {}s",
        loop_counter,
        (cur_lower + cur_higher) / 2.,
        analytical_solution - (cur_lower + cur_higher) / 2.
    );
    Ok(())
}
