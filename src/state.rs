//! Struct that contains a state

use nalgebra::SVector;
use nalgebra::Vector3;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReferenceFrame {
    GeocentricEquatorial,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct State {
    // state_vector is usually given in x, y, z, vx, vy, vz
    pub state_vector: SVector<f64, 6>,
    pub frame: ReferenceFrame,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EquatorialDirection {
    pub right_ascension_rad: f64,
    pub declination_rad: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EquatorialDirectionError {
    WrongReferenceFrame,
    ZeroPosition,
}

impl State {
    pub fn from_vector(state_vector: SVector<f64, 6>, frame: ReferenceFrame) -> Self {
        Self {
            state_vector,
            frame,
        }
    }

    pub fn from_array(values: [f64; 6], frame: ReferenceFrame) -> Self {
        Self {
            state_vector: SVector::from_row_slice(&values),
            frame,
        }
    }

    pub fn new(position: Vector3<f64>, velocity: Vector3<f64>, frame: ReferenceFrame) -> Self {
        Self::from_array(
            [
                position.x, position.y, position.z, velocity.x, velocity.y, velocity.z,
            ],
            frame,
        )
    }

    pub fn position(&self) -> Vector3<f64> {
        self.state_vector.fixed_rows::<3>(0).into_owned()
    }

    pub fn radius(&self) -> f64 {
        self.position().norm()
    }

    pub fn velocity(&self) -> Vector3<f64> {
        self.state_vector.fixed_rows::<3>(3).into_owned()
    }
    pub fn speed(&self) -> f64 {
        self.velocity().norm()
    }

    pub fn equatorial_direction(&self) -> Result<EquatorialDirection, EquatorialDirectionError> {
        if self.frame != ReferenceFrame::GeocentricEquatorial {
            return Err(EquatorialDirectionError::WrongReferenceFrame);
        }

        let position = self.position();
        let radius = position.norm();

        if radius == 0.0 {
            return Err(EquatorialDirectionError::ZeroPosition);
        }

        Ok(EquatorialDirection {
            right_ascension_rad: position
                .y
                .atan2(position.x)
                .rem_euclid(std::f64::consts::TAU),
            declination_rad: (position.z / radius).clamp(-1.0, 1.0).asin(),
        })
    }

    pub fn is_approx_eq(
        &self,
        other: &Self,
        position_tolerance: f64,
        velocity_tolerance: f64,
    ) -> bool {
        let position_error = (self.position() - other.position()).norm();
        let velocity_error = (self.velocity() - other.velocity()).norm();

        position_error <= position_tolerance && velocity_error <= velocity_tolerance
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_utils::assert_approx_eq;

    #[test]
    fn test_equatorial_direction_problem41b() {
        // Right  ascension of r = (-3000, -6000, -9000) is 243.4 deg
        let expected_alpha_deg: f64 = 243.4;
        let expected_delta_deg: f64 = -53.30;
        let calculated = State::from_array(
            [-3000.0, -6000.0, -9000.0, 0.0, 0.0, 0.0],
            ReferenceFrame::GeocentricEquatorial,
        )
        .equatorial_direction()
        .unwrap();
        assert_approx_eq(
            expected_alpha_deg.to_radians(),
            calculated.right_ascension_rad,
            1e-3,
            "Equatorial Direction RA test",
        );

        assert_approx_eq(
            expected_delta_deg.to_radians(),
            calculated.declination_rad,
            1e-4,
            "Equatorial Direction decl test",
        );
    }
}
