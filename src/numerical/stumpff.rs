//! Small helper library that contains the Stumpff function formulation necessary for the universal
//! variable approach to orbits.

pub fn c2(z: f64) -> f64 {
    //! Curtis' C function
    if z > 1e-6 {
        (1.0 - z.sqrt().cos()) / z
    } else if z < -1e-6 {
        (1.0 - (-z).sqrt().cosh()) / z
    } else if z == 0.0 {
        0.5
    } else {
        // small z branch
        0.5 - z / 24.0 + z.powi(2) / 720.0 - z.powi(3) / 40320.0
    }
}

pub fn c3(z: f64) -> f64 {
    //! Curtis' Z function
    if z > 1e-6 {
        (z.sqrt() - z.sqrt().sin()) / z.sqrt().powi(3)
    } else if z < -1e-6 {
        (-(-z).sqrt() + (-z).sqrt().sinh()) / (-z).sqrt().powi(3)
    } else if z == 0.0 {
        1.0 / 6.0
    } else {
        // small z branch
        1.0 / 6.0 - z / 120.0 + z.powi(2) / 5040.0 - z.powi(3) / 362880.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_utils::assert_approx_eq;

    #[test]
    fn test_stumpff_at_zero() {
        assert_approx_eq(c2(0.0), 0.5, 0.0, "C(0)");
        assert_approx_eq(c3(0.0), 1.0 / 6.0, 0.0, "S(0)");
    }

    #[test]
    fn test_stumpff_positive_argument() {
        let z = std::f64::consts::PI.powi(2);

        assert_approx_eq(c2(z), 0.20264236728467555, 1.0e-14, "C(pi^2)");
        assert_approx_eq(c3(z), 0.10132118364233778, 1.0e-14, "S(pi^2)");
    }

    #[test]
    fn test_stumpff_negative_argument() {
        assert_approx_eq(c2(-1.0), 0.5430806348152437, 1.0e-14, "C(-1)");
        assert_approx_eq(c3(-1.0), 0.17520119364380138, 1.0e-14, "S(-1)");
    }

    #[test]
    fn test_stumpff_small_positive_argument() {
        let z = 1.0e-8;

        assert_approx_eq(c2(z), 0.49999999958333335, 1.0e-14, "C(1e-8)");
        assert_approx_eq(c3(z), 0.16666666658333332, 1.0e-14, "S(1e-8)");
    }

    #[test]
    fn test_stumpff_small_negative_argument() {
        let z = -1.0e-8;

        assert_approx_eq(c2(z), 0.5000000004166667, 1.0e-14, "C(-1e-8)");
        assert_approx_eq(c3(z), 0.16666666675, 1.0e-14, "S(-1e-8)");
    }
}
