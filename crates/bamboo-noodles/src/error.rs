use bamboo_core::RegionParseError;
use bamboo_io::IoError;
use std::fmt;
use std::io;

#[derive(Debug)]
pub enum NoodlesError {
    Io(io::Error),
    ObjectStore(IoError),
    Region(RegionParseError),
    MissingIndex { path: String },
    MissingReference { name: String },
    Message(String),
}

impl fmt::Display for NoodlesError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(err) => write!(f, "{err}"),
            Self::ObjectStore(err) => write!(f, "{err}"),
            Self::Region(err) => write!(f, "{err}"),
            Self::MissingIndex { path } => {
                write!(f, "missing BAM index for indexed fetch: {path}")
            }
            Self::MissingReference { name } => {
                write!(f, "reference sequence not found in header: {name}")
            }
            Self::Message(message) => write!(f, "{message}"),
        }
    }
}

impl std::error::Error for NoodlesError {}

impl From<io::Error> for NoodlesError {
    fn from(value: io::Error) -> Self {
        Self::Io(value)
    }
}

impl From<RegionParseError> for NoodlesError {
    fn from(value: RegionParseError) -> Self {
        Self::Region(value)
    }
}

impl From<IoError> for NoodlesError {
    fn from(value: IoError) -> Self {
        Self::ObjectStore(value)
    }
}