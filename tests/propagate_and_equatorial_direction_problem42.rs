//! At a given instant, a spacecraft is 500 km above the earth, with a right ascension of 300 deg
//! and a declination of -2 0 deg relative to the geocentric equatorial frame. Its velocity is 10
//! km/s directly north, nomral to the equatorial plane. Find alpha and delta 30 min later.
//! { Ans.: Alpha = 120 deg, delta = -29.98 deg }

use curtis_orbital_mechanics::propagate::propagate_universal_variable;
use curtis_orbital_mechanics::state::{EquatorialDirection, ReferenceFrame, State};
use nalgebra::Vector3;

#[test]
fn test_propagation_from_equatorial_to_equatorial() {
    let start_equatorial_direction = EquatorialDirection {
        right_ascension_rad: 300.0_f64.to_radians(),
        declination_rad: (-60.0_f64).to_radians(),
    };
    let initial_position = start_equatorial_direction.to_vector(500.0 + 6378.0);

    let initial_velocity = Vector3::new(0.0, 0.0, 10.0);

    let initial_state = State::new(
        initial_position,
        initial_velocity,
        ReferenceFrame::GeocentricEquatorial,
    );

    let final_state = propagate_universal_variable(&initial_state, 30.0 * 60.0, 398600.0)
        .unwrap()
        .equatorial_direction()
        .unwrap();

    let expected_ra_deg: f64 = 120.0;
    let expected_decl_deg: f64 = -29.98;
    let tolerance = 1e-2;

    assert!(
        (final_state.right_ascension_rad - expected_ra_deg.to_radians()).abs() <= tolerance,
        "expected RA {expected_ra_deg}, got {}",
        final_state.right_ascension_rad.to_degrees(),
    );

    assert!(
        (final_state.declination_rad - expected_decl_deg.to_radians()).abs() <= tolerance,
        "expected declination {expected_decl_deg}, got {}",
        final_state.declination_rad.to_degrees(),
    );
}
