use crate::error::NoodlesError;
use noodles::sam::Header;
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::ReferenceSequence;
use std::num::NonZeroUsize;

/// Build a SAM header from reference name → length pairs.
pub fn header_from_references(refs: &[(String, u32)]) -> Result<Header, NoodlesError> {
    if refs.is_empty() {
        return Err(NoodlesError::Message(
            "BAM header requires at least one reference sequence".to_string(),
        ));
    }

    let mut builder = Header::builder();
    for (name, length) in refs {
        let length = *length as usize;
        let non_zero = NonZeroUsize::new(length).ok_or_else(|| {
            NoodlesError::Message(format!("reference '{name}' has invalid length {length}"))
        })?;
        builder = builder.add_reference_sequence(
            name.as_str(),
            Map::<ReferenceSequence>::new(non_zero),
        );
    }

    Ok(builder.build())
}