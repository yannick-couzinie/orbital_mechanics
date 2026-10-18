//! Write a program, like that for Example 2.18, to compute the trajectory of a spacecraft using the
//! restricted three-body equations of motion. Use the program to design a trajectory from the earth
//! to the earth-moon Lagrange point L4, starting at a 200-km altitude burnout point. The path
//! should take the coasting spacecraft to within 500km of L4 with a relative speed of not more than
//! 1 km/s
use std::{
    fs::{self, File},
    io::{BufWriter, Write},
    path::Path,
};

// We put the barycenter of the Earth-Moon system at the center, put the Earths center at negative x
// and the moon at positive x, then the Earth-L4 / Earth-Moon encloses a 60 degree angle and has a
// distance of 384.400km from both the earth and the moon

const EARTH_MOON_DISTANCE_KM: f64 = 3.844e5;
const EARTH_MASS_KG: f64 = 5.974e24;
const MOON_MASS_KG: f64 = 7.348e22;
const OMEGA_RAD_S: f64 = 2.66538e-6;
const MU_1_KM3_S3: f64 = 398620.0;
const MU_2_KM3_S3: f64 = 4903.02;
const PI1: f64 = EARTH_MASS_KG / (EARTH_MASS_KG + MOON_MASS_KG);
const PI2: f64 = MOON_MASS_KG / (EARTH_MASS_KG + MOON_MASS_KG);

// we don't move in the z plane so everything is two-dimensional
//
// We use the Earth–Moon barycentric rotating frame. If d is the Earth–Moon
// distance, the Earth is at (-PI2 * d, 0) and the Moon is at (PI1 * d, 0).
//
// Earth, Moon, and L4 form an equilateral triangle. Therefore, relative to
// Earth, L4 has coordinates
//
//     (d * cos(60°), d * sin(60°)) = (d / 2, sqrt(3) * d / 2).
//
// Translating this point into barycentric coordinates gives
//
//     x_L4 = -PI2 * d + d / 2 = (1/2 - PI2) * d
//     y_L4 = sqrt(3) * d / 2.

const TARGET_L4: nalgebra::SVector<f64, 2> = nalgebra::SVector::<f64, 2>::new(
    (0.5 - PI2) * EARTH_MOON_DISTANCE_KM,
    1.732_050_807_568_877_2_f64 * EARTH_MOON_DISTANCE_KM / 2.0,
);

use curtis_orbital_mechanics::numerical::integrate::{
    Integrator,
    adaptive_step_rk::{AdaptiveRkParameters, AdaptiveStepRkMethod},
    errors,
};

fn diffeq(_t: f64, y: &nalgebra::SVector<f64, 4>) -> nalgebra::SVector<f64, 4> {
    let r1: f64 = ((y[0] + PI2 * EARTH_MOON_DISTANCE_KM).powi(2) + y[1].powi(2)).sqrt();
    let r2: f64 = ((y[0] - PI1 * EARTH_MOON_DISTANCE_KM).powi(2) + y[1].powi(2)).sqrt();

    nalgebra::SVector::<f64, 4>::new(
        y[2],
        y[3],
        2.0 * OMEGA_RAD_S * y[3] + OMEGA_RAD_S.powi(2) * y[0]
            - MU_1_KM3_S3 * (y[0] + PI2 * EARTH_MOON_DISTANCE_KM) / r1.powi(3)
            - MU_2_KM3_S3 * (y[0] - PI1 * EARTH_MOON_DISTANCE_KM) / r2.powi(3),
        -2.0 * OMEGA_RAD_S * y[2] + OMEGA_RAD_S.powi(2) * y[1]
            - MU_1_KM3_S3 * y[1] / r1.powi(3)
            - MU_2_KM3_S3 * y[1] / r2.powi(3),
    )
}

fn distance_from_l4(state: &nalgebra::SVector<f64, 4>) -> f64 {
    ((state[0] - TARGET_L4[0]).powi(2) + (state[1] - TARGET_L4[1]).powi(2)).sqrt()
}

fn find_closest_point_to_l4(
    time: f64,
    integrator: &impl Integrator,
    initial_state: nalgebra::SVector<f64, 4>,
) -> Result<(f64, f64, nalgebra::SVector<f64, 4>), errors::IntegrationError<4>> {
    let integration_result = integrator.integrate(diffeq, (0.0, time), initial_state)?;

    let mut minimum_distance = f64::INFINITY;
    let mut minimum_index: usize = 0;
    for (index, state) in integration_result.states.iter().enumerate() {
        let distance = distance_from_l4(state);
        if distance < minimum_distance {
            minimum_distance = distance;
            minimum_index = index
        }
    }

    Ok((
        minimum_distance,
        integration_result.times[minimum_index],
        integration_result.states[minimum_index],
    ))
}

fn build_velocity_vector(angle_deg: f64, burnout_speed_km_s: f64) -> (f64, f64) {
    (
        -burnout_speed_km_s * angle_deg.to_radians().sin(),
        burnout_speed_km_s * angle_deg.to_radians().cos(),
    )
}

// calculate the distance to L4 and return distance time and state
fn evaluate(
    angle: f64,
    burnout_speed_km_s: f64,
    integrator: &impl Integrator,
    time: f64,
) -> (f64, f64, nalgebra::SVector<f64, 4>) {
    let (vx, vy) = build_velocity_vector(angle, burnout_speed_km_s);

    find_closest_point_to_l4(
        time,
        integrator,
        nalgebra::SVector::<f64, 4>::new(-PI2 * EARTH_MOON_DISTANCE_KM + 6578.0, 0.0, vx, vy),
    )
    .unwrap()
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let integrator = AdaptiveRkParameters {
        step: 20.0,
        absolute_tolerance: 1e-3,
        relative_tolerance: 1e-9,
        method: AdaptiveStepRkMethod::Rkf45,
    };

    // Given that we can reach the dark side of the moon in 3-4 days assume a time of 5 days is
    // enough to reach L4
    let time = 5.0 * 24.0 * 60.0 * 60.0; // assume that it will be quicker than 100m/s on average

    // let's use the same initial position as given in Example 2.17 and only alter the speed which
    // needs to be at least 10.86728km/s to reach L4 as per Example 2.17 and Fig. 2.37
    // we are basically doing the shooting method shooting for the best flight path angle at the
    // burnout point.
    //
    // We find the angle by bisection i.e we take gamma from -90 (all in positive x direction) over 0 (all in y direction) to 90 (all in
    // negative x direction) and bisect by targeting whichever is closer to L4.
    //
    // Given that the burnout velocity should be barely enough to get us to L4 in the first place we
    // will only check that we are within the 1km/s bound after finding a point that is close to
    // 500km

    let burnout_speed_km_s: f64 = 10.91271;

    let mut best_distance: f64 = f64::INFINITY;
    let mut best_angle_deg: f64 = f64::INFINITY;
    let mut best_coast_time: f64 = f64::INFINITY;

    // first perform a grid search
    for angle_step in -9000..=9000 {
        let angle_deg = angle_step as f64 / 100.0;
        let (vx, vy) = build_velocity_vector(angle_deg, burnout_speed_km_s);
        let (cur_dist, cur_time, _cur_state) = find_closest_point_to_l4(
            time,
            &integrator,
            nalgebra::SVector::<f64, 4>::new(-PI2 * EARTH_MOON_DISTANCE_KM + 6578.0, 0.0, vx, vy),
        )?;

        if cur_dist < best_distance {
            best_distance = cur_dist;
            best_angle_deg = angle_deg;
            best_coast_time = cur_time
        }
    }

    let mut angle_lower_deg: f64 = best_angle_deg - 0.1;
    let mut angle_upper_deg: f64 = best_angle_deg + 0.1;

    // do a ternary search, i.e. break the angle window into 4 points, 0, 1/3, 2/3 and 3/3. We
    // are trying to find the minimum of the distance. We roughly expect our window to be
    // monotonic so if the distance at 1/3 is lower than 2/3 then probably it is even higher at
    // 3/3 so continue the search from [0, 2/3].
    let angle_tolerance = 1e-4;
    while angle_upper_deg - angle_lower_deg > angle_tolerance {
        let third = (angle_upper_deg - angle_lower_deg) / 3.0;
        let left_angle_deg = angle_lower_deg + third;
        let right_angle_deg = angle_upper_deg - third;

        let (left_distance, _left_coast_time, _left_state) =
            evaluate(left_angle_deg, burnout_speed_km_s, &integrator, time);
        let (right_distance, _right_coast_time, _right_state) =
            evaluate(right_angle_deg, burnout_speed_km_s, &integrator, time);

        if left_distance < right_distance {
            angle_upper_deg = right_angle_deg;
        } else {
            angle_lower_deg = left_angle_deg;
        }
    }
    println!(
        "Grid best angle: {}°, distance: {} km, time: {} d",
        best_angle_deg,
        best_distance,
        best_coast_time / (24.0 * 60.0 * 60.0)
    );

    let refined_angle_deg = (angle_lower_deg + angle_upper_deg) / 2.0;

    let (refined_distance, refined_coast_time, refined_state) =
        evaluate(refined_angle_deg, burnout_speed_km_s, &integrator, time);

    println!(
        "Refined result: angle = {}°, distance = {} km, relative speed = {} km/s, coast time = {} d",
        refined_angle_deg,
        refined_distance,
        (refined_state[2].powi(2) + refined_state[3].powi(2)).sqrt(),
        refined_coast_time / (24.0 * 60.0 * 60.0),
    );

    let (vx, vy) = build_velocity_vector(refined_angle_deg, burnout_speed_km_s);

    let initial_state =
        nalgebra::SVector::<f64, 4>::new(-PI2 * EARTH_MOON_DISTANCE_KM + 6578.0, 0.0, vx, vy);

    let trajectory = integrator.integrate(diffeq, (0.0, refined_coast_time), initial_state)?;

    let output_directory = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("examples/problem247_trajectory_to_l4_artifacts");

    fs::create_dir_all(&output_directory)?;

    let file = File::create(output_directory.join("trajectory.csv"))?;
    let mut writer = BufWriter::new(file);

    writeln!(writer, "time,x,y,vx,vy")?;

    for (time, state) in trajectory.times.iter().zip(&trajectory.states) {
        writeln!(
            writer,
            "{time},{},{},{},{}",
            state[0], state[1], state[2], state[3],
        )?;
    }

    Ok(())
}
