//! One-dimensional Runge-Kutta numerical integration with methods from RK1 to RK4.

/// The RkTypes we can use
pub enum RkTypes {
    /// First-order Runge-Kutta numerical integration
    RK1,
    /// Second-order Runge-Kutta numerical integration
    RK2,
    /// Third-order Runge-Kutta numerical integration
    RK3,
    /// Fourth-order Runge-Kutta numerical integration
    RK4,
}

impl RkTypes {
    fn n_stages(&self) -> usize {
        match self {
            RkTypes::RK1 => 1,
            RkTypes::RK2 => 2,
            RkTypes::RK3 => 3,
            RkTypes::RK4 => 4,
        }
    }

    fn a(&self) -> Vec<f64> {
        match self {
            RkTypes::RK1 => vec![0.],
            RkTypes::RK2 => vec![0., 1.],
            RkTypes::RK3 => vec![0., 0.5, 1.],
            RkTypes::RK4 => vec![0., 0.5, 0.5, 1.],
        }
    }

    fn b(&self) -> Vec<Vec<f64>> {
        match self {
            RkTypes::RK1 => vec![vec![0.]],
            RkTypes::RK2 => vec![vec![0., 1.]],
            RkTypes::RK3 => vec![vec![0., 0.], vec![0., 0.5], vec![-1., 2.]],
            RkTypes::RK4 => {
                vec![
                    vec![0., 0., 0.],
                    vec![0.5, 0., 0.],
                    vec![0., 0.5, 0.],
                    vec![0., 0., 1.],
                ]
            }
        }
    }

    fn c(&self) -> Vec<f64> {
        match self {
            RkTypes::RK1 => vec![1.],
            RkTypes::RK2 => vec![0.5, 0.5],
            RkTypes::RK3 => vec![1. / 6., 2. / 3., 1. / 6.],
            RkTypes::RK4 => vec![1. / 6., 1. / 3., 1. / 3., 1. / 6.],
        }
    }
}

pub fn rk1_4(
    // The function to integrate, the input and output vector sizes (y) need to be the same, the
    // input is t, y
    ode_function: &dyn Fn(f64, &Vec<f64>) -> Vec<f64>,
    // The timespan with the first entry being the start and the second the endpoint
    tspan: (f64, f64),
    // Column vector of the initial values of the vector y
    y0: Vec<f64>,
    // Time step
    h: f64,
    rk: RkTypes,
) -> (Vec<f64>, Vec<Vec<f64>>) {
    let t0 = tspan.0;
    let tf = tspan.1;
    let mut step_h = h;
    let mut t = t0;
    let mut y = y0.clone();
    let mut tout = vec![t];
    let mut yout = vec![y.clone()];

    while t < tf {
        let ti = t;
        let yi = y.clone();
        let mut f: Vec<Vec<f64>> = Vec::new(); // the [n_stages, dim(y)] sized matrix with the function values

        // evaluate the time derivates at the 'n_stages' points within the current interval
        for i in 0..rk.n_stages() {
            let t_inner = ti + rk.a().get(i).unwrap() * h;
            let mut y_inner = yi.clone();
            if i == 0 {
                let f_eval = ode_function(t_inner, &y_inner);
                for f_entry in &f_eval {
                    // f is an empty vector at this point
                    // we want to push to f(:, i) where the : here is iterated in this for-loop
                    // and i = 0 (i.e. first entry)
                    f.push(vec![*f_entry])
                }
            } else {
                for j in 0..(i - 1) {
                    // this will not run on i=1
                    for (k, y_entry) in y_inner.iter_mut().enumerate() {
                        *y_entry += step_h
                            * rk.b().get(i).unwrap().get(j).unwrap()
                            * f.get(k).unwrap().get(j).unwrap();
                    }
                }

                let f_eval = ode_function(t_inner, &y_inner);
                for (j, f_entry) in f_eval.into_iter().enumerate() {
                    // f is an empty vector at this point
                    // we want to push to f(:, i) where the : here is iterated in this for-loop
                    // and i = 0 (i.e. first entry)
                    f[j].push(f_entry)
                }
            }
        }
        step_h = step_h.min(tf - t);
        t += step_h;
        let mut newy_vec = Vec::new();
        // loop over dimensions of y
        for (k1, y_entry) in yi.into_iter().enumerate() {
            // inner product  loop between f and c
            let mut new_y_entry = y_entry;
            for (k2, fentry) in f[k1].iter().enumerate() {
                new_y_entry += step_h * *fentry * rk.c()[k2];
            }
            newy_vec.push(new_y_entry)
        }
        y = newy_vec.clone();
        tout.push(t);
        yout.push(newy_vec);
    }
    (tout, yout)
}
