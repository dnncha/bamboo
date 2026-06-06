use std::fmt;

#[derive(Debug)]
pub enum IoError {
    InvalidUri(String),
    ObjectStore(object_store::Error),
    Local(std::io::Error),
    Message(String),
}

impl fmt::Display for IoError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidUri(uri) => write!(f, "invalid URI '{uri}'"),
            Self::ObjectStore(err) => write!(f, "{err}"),
            Self::Local(err) => write!(f, "{err}"),
            Self::Message(message) => write!(f, "{message}"),
        }
    }
}

impl std::error::Error for IoError {}

impl From<std::io::Error> for IoError {
    fn from(value: std::io::Error) -> Self {
        Self::Local(value)
    }
}

impl From<object_store::Error> for IoError {
    fn from(value: object_store::Error) -> Self {
        Self::ObjectStore(value)
    }
}