use nalgebra;

pub struct IntegrationResult<const N: usize> {
    pub times: Vec<f64>,
    pub states: Vec<nalgebra::SVector<f64, N>>,
}
