use nalgebra::SVector;

#[derive(Debug)]
pub enum IntegrationError<const N: usize> {
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
    NonFiniteDerivative {
        derivative: SVector<f64, N>,
        state: SVector<f64, N>,
        time: f64,
    },
    NonFiniteState {
        state: SVector<f64, N>,
    },
}
