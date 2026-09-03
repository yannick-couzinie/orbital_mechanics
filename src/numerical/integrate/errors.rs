use nalgebra::SVector;
use snafu::Snafu;

#[derive(Debug, Snafu)]
pub enum IntegrationError<const N: usize> {
    #[snafu(display("invalid integration step size: {step}"))]
    InvalidStepSize { step: f64 },
    #[snafu(display("Supply non-zero, positive, finite tolerance RKF45."))]
    AdaptiveParametersIncomplete {},
    #[snafu(display("The error has become non finite during integration {error}"))]
    NonFiniteError { error: f64 },
    #[snafu(display("invalid integration time span: {start} to {end}"))]
    InvalidTimeSpan { start: f64, end: f64 },
    #[snafu(display("invalid integration step (does not advance time): {time} + {step}"))]
    StepDoesNotAdvanceTime { time: f64, step: f64 },
    #[snafu(display(
        "derivative is non-finite: derivative={derivative}, state={state}, time={time}"
    ))]
    NonFiniteDerivative {
        derivative: SVector<f64, N>,
        state: SVector<f64, N>,
        time: f64,
    },
    #[snafu(display("state is non-finite: state={state}"))]
    NonFiniteState { state: SVector<f64, N> },
}
