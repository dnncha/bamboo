use std::fmt;

/// A genomic fetch region.
///
/// Coordinates follow the Python/pysam convention: 0-based, half-open intervals.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FetchRegion {
    pub reference_name: String,
    pub start: Option<u32>,
    pub end: Option<u32>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RegionParseError {
    Empty,
    InvalidFormat(String),
    InvalidCoordinate(String),
}

impl fmt::Display for RegionParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => write!(f, "region string is empty"),
            Self::InvalidFormat(value) => write!(
                f,
                "invalid region '{value}': expected 'chrom', 'chrom:start-end', or 'chrom:start'"
            ),
            Self::InvalidCoordinate(value) => write!(f, "invalid coordinate in region: {value}"),
        }
    }
}

impl std::error::Error for RegionParseError {}

impl FetchRegion {
    /// Parse a samtools-style region string (1-based, inclusive start) into a
    /// Python/pysam-style 0-based half-open interval.
    pub fn from_samtools_region(region: &str) -> Result<Self, RegionParseError> {
        let region = region.trim();
        if region.is_empty() {
            return Err(RegionParseError::Empty);
        }

        let (reference_name, interval) = match region.split_once(':') {
            Some((reference_name, interval)) => (reference_name.to_string(), Some(interval)),
            None => (region.to_string(), None),
        };

        if reference_name.is_empty() {
            return Err(RegionParseError::InvalidFormat(region.to_string()));
        }

        let (start, end) = match interval {
            None => (None, None),
            Some(interval) => {
                if let Some((start, end)) = interval.split_once('-') {
                    let start = parse_samtools_start(start)?;
                    let end = parse_samtools_end(end)?;
                    (Some(start), Some(end))
                } else {
                    let start = parse_samtools_start(interval)?;
                    (Some(start), None)
                }
            }
        };

        Ok(Self {
            reference_name,
            start,
            end,
        })
    }

    /// Build a samtools/noodles region string from Python coordinates.
    pub fn to_samtools_region(&self) -> String {
        match (self.start, self.end) {
            (None, None) => self.reference_name.clone(),
            (Some(start), None) => format!("{}:{}", self.reference_name, start + 1),
            (Some(start), Some(end)) => {
                format!("{}:{}-{}", self.reference_name, start + 1, end)
            }
            (None, Some(end)) => format!("{}:1-{}", self.reference_name, end),
        }
    }
}

fn parse_samtools_start(value: &str) -> Result<u32, RegionParseError> {
    let one_based: u32 = value
        .parse()
        .map_err(|_| RegionParseError::InvalidCoordinate(value.to_string()))?;
    if one_based == 0 {
        return Err(RegionParseError::InvalidCoordinate(value.to_string()));
    }
    Ok(one_based - 1)
}

fn parse_samtools_end(value: &str) -> Result<u32, RegionParseError> {
    value
        .parse()
        .map_err(|_| RegionParseError::InvalidCoordinate(value.to_string()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_whole_contig() {
        let region = FetchRegion::from_samtools_region("chr1").unwrap();
        assert_eq!(region.reference_name, "chr1");
        assert_eq!(region.start, None);
        assert_eq!(region.end, None);
    }

    #[test]
    fn parses_interval() {
        let region = FetchRegion::from_samtools_region("chr1:1000-2000").unwrap();
        assert_eq!(region.reference_name, "chr1");
        assert_eq!(region.start, Some(999));
        assert_eq!(region.end, Some(2000));
    }
}