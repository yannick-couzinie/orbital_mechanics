use std::fmt::Display;

#[track_caller]
pub(crate) fn assert_approx_eq(actual: f64, expected: f64, tolerance: f64, context: impl Display) {
    let error = (actual - expected).abs();

    assert!(
        error <= tolerance,
        "{context}: expected {expected}, got {actual}, absolute error {error}, tolerance {tolerance}"
    );
}
