use curtis_orbital_mechanics::numerical::rkf1_4::{RkTypes, rk1_4};
use curtis_orbital_mechanics::test_utils::assert_approx_eq;

#[test]
fn problem118_rksolvers() {
    // Run the complete solver and compare with the analytical result.
    // this is the first order ssytem for d4y/dt4 + 2d2y/dt2 + y = 0 with initial conditions y=1 and
    // dy/dt = d2y/dt2 = d3y/dt3 = 0 at t=0 solved for y(20) which should result in 9.545
    // the analytical solution is actually t.cos() + (t / 2.0) * t.sin() so we can compare against
    // that
    //
    // run the tests only for RK4 as that should be proof enough that the solvers work, and getting
    // the same accuracy on RK1 would require thousand-fold of steps
    let (times, states) = rk1_4(
        &|_t, y: &Vec<f64>| vec![y[1], y[2], y[3], -y[0] - 2.0 * y[2]],
        (0.0, 20.0),
        vec![1., 0., 0., 0.],
        0.01,
        RkTypes::RK4,
    );

    for (i, t_entry) in times.into_iter().enumerate() {
        // get 0 since the first entry in the vector is equal to y
        let expected = t_entry.cos() + (t_entry * 0.5 * t_entry.sin());
        let actual = states.get(i).unwrap().first().copied().unwrap();

        assert_approx_eq(
            actual,
            expected,
            1.0e-6,
            format!("problem 1.18 failed using RK4 at step {i}, t={t_entry}"),
        );
    }
}
