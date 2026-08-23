#[derive(Debug)]
pub enum IntegrateError {
    InvalidStepSize { step: f64 },
    InvalidTimeSpan { start: f64, end: f64 },
    StepDoesNotAdvanceTime { time: f64, step: f64 },
    MismatchedDerivativeDimension,
    NonFiniteDerivative,
}
