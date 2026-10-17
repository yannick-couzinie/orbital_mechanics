use snafu::Snafu;

#[derive(Debug, Snafu)]
pub enum RootError {
    #[snafu(display(
        "Function evaluated to something non-finite with input {input}, output {output}."
    ))]
    NonFiniteEvaluation { input: f64, output: f64 },

    #[snafu(display(
        "Function derivative evaluated to something non-finite with input {input}, output {output}."
    ))]
    NonFiniteDerivative { input: f64, output: f64 },

    #[snafu(display("The root is not correctly bracketed between {lower} to {upper}."))]
    RootNotBracketed { lower: f64, upper: f64 },

    #[snafu(display(
        "The given interval interval {lower} to {upper} contains invalid numbers, both ends must be finite."
    ))]
    InvalidSearchInterval { lower: f64, upper: f64 },

    #[snafu(display(
        "The tolerance needs to be positive and finite, the given tolerance {tolerance} is invalid."
    ))]
    InvalidTolerance { tolerance: f64 },
}
