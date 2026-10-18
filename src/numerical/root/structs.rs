use std::fmt;

pub struct RootResult {
    pub root: f64,
    pub residual: f64,
    pub iterations: usize,
}

impl fmt::Display for RootResult {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "root: {}, residual: {}, iterations: {}",
            self.root, self.residual, self.iterations
        )
    }
}
