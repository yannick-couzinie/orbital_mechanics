mod numerical;

use crate::numerical::rkf1_4::RkTypes;
use crate::numerical::rkf1_4::rk1_4;

fn main() {
    // this is the first order ssytem for d4y/dt4 + 2d2y/dt2 + y = 0 with initial conditions y=1 and
    // dy/dt = d2y/dt2 = d3y/dt3 = 0 at t=0 solved for y(20) which should result in 9.545
    let result = rk1_4(
        &|_t, y: &Vec<f64>| vec![y[1], y[2], y[3], -y[0] - 2.0 * y[2]],
        (0.0, 20.0),
        vec![1., 0., 0., 0.],
        0.0001,
        RkTypes::RK2,
    );
    println!("{:?}", result.1.last().unwrap()[0]);
}
