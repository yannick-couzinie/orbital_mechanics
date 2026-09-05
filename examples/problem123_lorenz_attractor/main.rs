//! Use the RKF45 solver to solve the nonlinear Lorenz equations, due to the American meterologist
//! and mathematician E.N. Lorenz
//!
//! dot(x) = sigma*(y-x)
//! dot(y) = x*(rho-z) - y
//! dot(z) = x*y - beta * z
//!
//! Star off by using the values that Lorenz (1963) used in his paper (namely, sigma = 10, beta =
//! 8/3, and rho = 28). For initial conditions use x = 0, y = 1, and z = 0 at t=0. Let t range to a
//! value of at least 20. Plot the phase trajectory x = x(t), y = y(t), zz(t) in three dimensions to
//! see the now-famous "Lorenz attractor".

const SIGMA: f64 = 10.;
const RHO: f64 = 28.;
const BETA: f64 = 8. / 3.;

use curtis_orbital_mechanics::numerical::integrate::{
    Integrator,
    adaptive_step_rk::{AdaptiveRkParameters, AdaptiveStepRkMethod},
};

use std::{
    fs::{self, File},
    io::{BufWriter, Write},
    path::Path,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let integrator = AdaptiveRkParameters {
        step: 20.0,
        tolerance: 1e-10,
        method: AdaptiveStepRkMethod::Rkf45,
    };
    let integration_result = integrator.integrate(
        |_t, y| {
            nalgebra::SVector::<f64, 3>::new(
                SIGMA * (y[1] - y[0]),
                y[0] * (RHO - y[2]) - y[1],
                y[0] * y[1] - BETA * y[2],
            )
        },
        (0.0, 50.),
        nalgebra::SVector::<f64, 3>::new(0., 1., 0.),
    )?;

    let output_directory = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("examples/problem123_lorenz_attractor/artifacts");

    fs::create_dir_all(&output_directory)?;

    let file = File::create(output_directory.join("trajectory.csv"))?;
    let mut writer = BufWriter::new(file);

    writeln!(writer, "time,x,y,z")?;

    for (time, state) in integration_result
        .times
        .iter()
        .zip(&integration_result.states)
    {
        writeln!(writer, "{time},{},{},{}", state[0], state[1], state[2])?;
    }

    Ok(())
}
