//! For the sun-earth system, find the distance of the collinear Lagrange points L1, L2, and L3 from
//! the barycenter.

const M1_KG: f64 = 5.974e24; // Earth
const M2_KG: f64 = 1.989e30; // Sun
const PI1: f64 = M1_KG / (M1_KG + M2_KG);
const PI2: f64 = 1.0 - PI1;
const AU_KM: f64 = 149.6e6;
const TOLERANCE: f64 = 1.0 / AU_KM;

fn function_to_root(xi: f64) -> f64 {
    PI1 * (xi + PI2) / ((xi + PI2).abs().powi(3)) + PI2 * (xi - PI1) / ((xi - PI1).abs().powi(3))
        - xi
}

// find the root of function_to_root starting with the guesses as input
fn bisect_for_root(lower: f64, upper: f64) -> (f64, f64) {
    let middle: f64 = (lower + upper) / 2.0;
    if function_to_root(lower) * function_to_root(middle) > 0.0 {
        (middle, upper)
    } else {
        (lower, middle)
    }
}

fn find_root(initial_lower: f64, initial_upper: f64) -> (f64, f64) {
    let iterations: usize =
        (((initial_upper - initial_lower).abs() / TOLERANCE).log2()).ceil() as usize;
    let mut lower: f64 = initial_lower;
    let mut upper: f64 = initial_upper;
    for _ in 0..iterations {
        (lower, upper) = bisect_for_root(lower, upper);
    }
    (lower, upper)
}

fn main() {
    let mut root: f64;
    (root, _) = find_root(-1.1, -1.001);
    println!("L2 (beyond earth) is at {}km", root * AU_KM);

    (root, _) = find_root(-0.999, -0.9);
    println!("L1 (between earth and sun) is at {}km", root * AU_KM);

    // L3
    (root, _) = find_root(0.9, 1.1);
    println!("L3 (beyond sun) is at {}km", root * AU_KM);
}
