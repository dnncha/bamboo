use thiserror::Error;

#[derive(Debug, Error)]
pub enum HtslibError {
    #[error(transparent)]
    Io(#[from] std::io::Error),
    #[error(transparent)]
    Htslib(#[from] rust_htslib::errors::Error),
    #[error("{0}")]
    Message(String),
}