#[derive(Debug)]
pub enum IntegrationError {
    InvalidStepSize {
        step: f64,
    },
    InvalidTimeSpan {
        start: f64,
        end: f64,
    },
    StepDoesNotAdvanceTime {
        time: f64,
        step: f64,
    },
    MismatchedDerivativeDimension {
        state_dimension: usize,
        derivative_dimension: usize,
    },
    NonFiniteDerivative {
        derivative: Vec<f64>,
        state: Vec<f64>,
        time: f64,
    },
    NonFiniteState {
        state: Vec<f64>,
    },
}
