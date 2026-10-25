use crate::state::State;
use std::fmt::Display;

#[track_caller]
pub(crate) fn assert_approx_eq(actual: f64, expected: f64, tolerance: f64, context: impl Display) {
    let error = (actual - expected).abs();

    assert!(
        error <= tolerance,
        "{context}: expected {expected}, got {actual}, absolute error {error}, tolerance {tolerance}"
    );
}

#[track_caller]
pub(crate) fn assert_state_approx_eq(
    actual: &State,
    expected: &State,
    position_tolerance: f64,
    velocity_tolerance: f64,
    context: &str,
) {
    let position_error = (actual.position() - expected.position()).norm();
    let velocity_error = (actual.velocity() - expected.velocity()).norm();

    assert!(
        position_error <= position_tolerance && velocity_error <= velocity_tolerance,
        "{context}: position error {position_error}, \
         velocity error {velocity_error}"
    );
}
