use crate::error::NoodlesError;
use bamboo_core::{BamColumn, BamScanOptions, BamTable, TagValue};
use noodles::bam as bam;
use noodles::sam::Header;
use noodles::sam::alignment::RecordBuf;

use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::record::MappingQuality;
use noodles::sam::alignment::record::cigar::Op;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::Cigar;
use noodles::sam::alignment::record_buf::Data;
use noodles::sam::alignment::record_buf::QualityScores;
use noodles::sam::alignment::record_buf::Sequence;
use noodles::sam::alignment::record_buf::data::field::Value;

/// A single parsed alignment record exposed to Python.
#[derive(Debug, Clone)]
pub struct AlignedRecord {
    pub query_name: Option<String>,
    pub flag: u16,
    pub reference_name: Option<String>,
    pub reference_start: Option<i64>,
    pub mapping_quality: Option<u8>,
    pub cigar: String,
    pub mate_reference_name: Option<String>,
    pub mate_reference_start: Option<i64>,
    pub template_length: Option<i32>,
    pub query_sequence: Option<String>,
    pub query_qualities: Option<String>,
    pub tags: Vec<(String, TagValue)>,
}

impl AlignedRecord {
    pub fn from_record_buf(header: &Header, record: &RecordBuf, options: &BamScanOptions) -> Self {
        let needs = FieldNeeds::from_options(options);

        let reference_name = if needs.reference_name {
            record
                .reference_sequence_id()
                .and_then(|id| header.reference_sequences().get_index(id))
                .map(|(name, _)| name.to_string())
        } else {
            None
        };

        let mate_reference_name = if needs.mate_reference_name {
            record
                .mate_reference_sequence_id()
                .and_then(|id| header.reference_sequences().get_index(id))
                .map(|(name, _)| name.to_string())
        } else {
            None
        };

        let tags = if needs.tags {
            options
                .tags
                .iter()
                .map(|name| {
                    let value = record
                        .data()
                        .get(&Tag::from([name.as_bytes()[0], name.as_bytes()[1]]))
                        .map(convert_tag_value)
                        .unwrap_or(TagValue::Missing);
                    (name.clone(), value)
                })
                .collect()
        } else {
            Vec::new()
        };

        Self {
            query_name: if needs.query_name {
                record.name().map(|name| name.to_string())
            } else {
                None
            },
            flag: if needs.flag {
                record.flags().bits()
            } else {
                0
            },
            reference_name,
            reference_start: if needs.reference_start {
                record
                    .alignment_start()
                    .map(|position| position.get() as i64 - 1)
            } else {
                None
            },
            mapping_quality: if needs.mapping_quality {
                record.mapping_quality().map(|quality| quality.get())
            } else {
                None
            },
            cigar: if needs.cigar {
                cigar_to_string(record.cigar())
            } else {
                String::new()
            },
            mate_reference_name,
            mate_reference_start: if needs.mate_reference_start {
                record
                    .mate_alignment_start()
                    .map(|position| position.get() as i64 - 1)
            } else {
                None
            },
            template_length: if needs.template_length {
                Some(record.template_length())
            } else {
                None
            },
            query_sequence: if needs.query_sequence {
                bytes_to_optional_string(record.sequence().as_ref())
            } else {
                None
            },
            query_qualities: if needs.query_qualities {
                bytes_to_optional_string(record.quality_scores().as_ref())
            } else {
                None
            },
            tags,
        }
    }

    pub fn from_bam_record(header: &Header, record: &bam::Record, options: &BamScanOptions) -> Self {
        let needs = FieldNeeds::from_options(options);

        let reference_name = if needs.reference_name {
            record
                .reference_sequence_id()
                .transpose()
                .ok()
                .flatten()
                .and_then(|id| header.reference_sequences().get_index(id))
                .map(|(name, _)| name.to_string())
        } else {
            None
        };

        let mate_reference_name = if needs.mate_reference_name {
            record
                .mate_reference_sequence_id()
                .transpose()
                .ok()
                .flatten()
                .and_then(|id| header.reference_sequences().get_index(id))
                .map(|(name, _)| name.to_string())
        } else {
            None
        };

        Self {
            query_name: if needs.query_name {
                record.name().map(|name| name.to_string())
            } else {
                None
            },
            flag: if needs.flag {
                record.flags().bits()
            } else {
                0
            },
            reference_name,
            reference_start: if needs.reference_start {
                record
                    .alignment_start()
                    .transpose()
                    .ok()
                    .flatten()
                    .map(|position| position.get() as i64 - 1)
            } else {
                None
            },
            mapping_quality: if needs.mapping_quality {
                record.mapping_quality().map(|quality| quality.get())
            } else {
                None
            },
            cigar: if needs.cigar {
                bam_cigar_to_string(record.cigar())
            } else {
                String::new()
            },
            mate_reference_name,
            mate_reference_start: if needs.mate_reference_start {
                record
                    .mate_alignment_start()
                    .transpose()
                    .ok()
                    .flatten()
                    .map(|position| position.get() as i64 - 1)
            } else {
                None
            },
            template_length: if needs.template_length {
                Some(record.template_length())
            } else {
                None
            },
            query_sequence: if needs.query_sequence {
                bam_sequence_to_optional_string(record.sequence())
            } else {
                None
            },
            query_qualities: if needs.query_qualities {
                bam_quality_to_optional_string(record.quality_scores())
            } else {
                None
            },
            tags: Vec::new(),
        }
    }

    pub fn passes_filters(&self, options: &BamScanOptions) -> bool {
        if let Some(min_mapq) = options.min_mapq {
            match self.mapping_quality {
                Some(mapq) if mapq >= min_mapq => {}
                _ => return false,
            }
        }

        if let Some(reference_name) = &options.reference_name {
            match &self.reference_name {
                Some(name) if name == reference_name => {}
                _ => return false,
            }
        }

        if let Some(region) = &options.region {
            let record_ref = self.reference_name.as_deref().unwrap_or("");
            if record_ref != region.reference_name {
                return false;
            }
            if let Some(start) = region.start {
                let pos = self.reference_start.unwrap_or(-1);
                if pos < start as i64 {
                    return false;
                }
            }
            if let Some(end) = region.end {
                let pos = self.reference_start.unwrap_or(i64::MAX);
                if pos >= end as i64 {
                    return false;
                }
            }
        }

        true
    }

    pub fn to_record_buf(&self, header: &Header) -> Result<RecordBuf, NoodlesError> {
        let mut builder = RecordBuf::builder().set_flags(Flags::from_bits_truncate(self.flag));

        if let Some(name) = &self.query_name {
            builder = builder.set_name(name.as_str());
        }

        if let Some(reference_name) = &self.reference_name {
            let id = header
                .reference_sequences()
                .get_index_of(reference_name.as_bytes())
                .ok_or_else(|| NoodlesError::MissingReference {
                    name: reference_name.clone(),
                })?;
            builder = builder.set_reference_sequence_id(id);
        }

        if let Some(start) = self.reference_start {
            let position = noodles::core::Position::try_from((start + 1) as usize).map_err(|err| {
                NoodlesError::Message(format!("invalid alignment start {start}: {err}"))
            })?;
            builder = builder.set_alignment_start(position);
        }

        if let Some(mapq) = self.mapping_quality {
            let quality = MappingQuality::try_from(mapq).map_err(|err| {
                NoodlesError::Message(format!("invalid mapping quality {mapq}: {err}"))
            })?;
            builder = builder.set_mapping_quality(quality);
        }

        if !self.cigar.is_empty() && self.cigar != "*" {
            builder = builder.set_cigar(parse_cigar(&self.cigar)?);
        }

        if let Some(mate_reference_name) = &self.mate_reference_name {
            if let Some(id) = header
                .reference_sequences()
                .get_index_of(mate_reference_name.as_bytes())
            {
                builder = builder.set_mate_reference_sequence_id(id);
            }
        }

        if let Some(start) = self.mate_reference_start {
            let position = noodles::core::Position::try_from((start + 1) as usize).map_err(|err| {
                NoodlesError::Message(format!("invalid mate alignment start {start}: {err}"))
            })?;
            builder = builder.set_mate_alignment_start(position);
        }

        if let Some(template_length) = self.template_length {
            builder = builder.set_template_length(template_length);
        }

        if let Some(sequence) = &self.query_sequence {
            builder = builder.set_sequence(Sequence::from(sequence.as_bytes()));
        }

        if let Some(qualities) = &self.query_qualities {
            builder = builder.set_quality_scores(QualityScores::from(qualities.as_bytes().to_vec()));
        }

        let mut data = Data::default();
        for (name, value) in &self.tags {
            if matches!(value, TagValue::Missing) {
                continue;
            }
            let tag = tag_from_name(name)?;
            data.insert(tag, tag_value_to_noodles(value)?);
        }
        builder = builder.set_data(data);

        Ok(builder.build())
    }

    /// Reconstruct a record from a columnar table row (used by fetch iterators).
    pub fn from_table_row(table: &BamTable, row: usize) -> Self {
        let has = |column: BamColumn| table.columns.iter().any(|c| *c == column);

        Self {
            query_name: if has(BamColumn::QueryName) {
                table.qname.get(row).cloned().flatten()
            } else {
                None
            },
            flag: if has(BamColumn::Flag) {
                table.flag.get(row).copied().unwrap_or(0)
            } else {
                0
            },
            reference_name: if has(BamColumn::ReferenceName) {
                table.rname.get(row).cloned().flatten()
            } else {
                None
            },
            reference_start: if has(BamColumn::Position) {
                table
                    .pos
                    .get(row)
                    .cloned()
                    .flatten()
                    .map(|value| value as i64)
            } else {
                None
            },
            mapping_quality: if has(BamColumn::MappingQuality) {
                table.mapq.get(row).copied().flatten()
            } else {
                None
            },
            cigar: if has(BamColumn::Cigar) {
                table.cigar.get(row).cloned().unwrap_or_default()
            } else {
                String::new()
            },
            mate_reference_name: if has(BamColumn::MateReferenceName) {
                table.rnext.get(row).cloned().flatten()
            } else {
                None
            },
            mate_reference_start: if has(BamColumn::MatePosition) {
                table
                    .pnext
                    .get(row)
                    .cloned()
                    .flatten()
                    .map(|value| value as i64)
            } else {
                None
            },
            template_length: if has(BamColumn::TemplateLength) {
                table.tlen.get(row).copied().flatten()
            } else {
                None
            },
            query_sequence: if has(BamColumn::Sequence) {
                table.seq.get(row).cloned().flatten()
            } else {
                None
            },
            query_qualities: if has(BamColumn::Quality) {
                table.qual.get(row).cloned().flatten()
            } else {
                None
            },
            tags: Vec::new(),
        }
    }

    pub fn append_to_table(&self, table: &mut BamTable, options: &BamScanOptions) {
        for column in &options.columns {
            match column {
                BamColumn::QueryName => table.qname.push(self.query_name.clone()),
                BamColumn::Flag => table.flag.push(self.flag),
                BamColumn::ReferenceName => table.rname.push(self.reference_name.clone()),
                BamColumn::Position => table
                    .pos
                    .push(self.reference_start.map(|value| value as i32)),
                BamColumn::MappingQuality => table.mapq.push(self.mapping_quality),
                BamColumn::Cigar => table.cigar.push(self.cigar.clone()),
                BamColumn::MateReferenceName => {
                    table.rnext.push(self.mate_reference_name.clone())
                }
                BamColumn::MatePosition => table
                    .pnext
                    .push(self.mate_reference_start.map(|value| value as i32)),
                BamColumn::TemplateLength => table.tlen.push(self.template_length),
                BamColumn::Sequence => table.seq.push(self.query_sequence.clone()),
                BamColumn::Quality => table.qual.push(self.query_qualities.clone()),
            }
        }

        for (index, tag_name) in options.tags.iter().enumerate() {
            let value = self
                .tags
                .iter()
                .find(|(name, _)| name == tag_name)
                .map(|(_, value)| value.clone())
                .unwrap_or(TagValue::Missing);
            table.tags[index].values.push(value);
        }
    }
}

struct FieldNeeds {
    query_name: bool,
    flag: bool,
    reference_name: bool,
    reference_start: bool,
    mapping_quality: bool,
    cigar: bool,
    mate_reference_name: bool,
    mate_reference_start: bool,
    template_length: bool,
    query_sequence: bool,
    query_qualities: bool,
    tags: bool,
}

impl FieldNeeds {
    fn from_options(options: &BamScanOptions) -> Self {
        let filter_needs_reference = options.region.is_some() || options.reference_name.is_some();
        let filter_needs_mapq = options.min_mapq.is_some();
        let filter_needs_position = options.region.is_some();

        Self {
            query_name: options.wants_column(BamColumn::QueryName),
            flag: options.wants_column(BamColumn::Flag),
            reference_name: options.wants_column(BamColumn::ReferenceName) || filter_needs_reference,
            reference_start: options.wants_column(BamColumn::Position) || filter_needs_position,
            mapping_quality: options.wants_column(BamColumn::MappingQuality) || filter_needs_mapq,
            cigar: options.wants_column(BamColumn::Cigar),
            mate_reference_name: options.wants_column(BamColumn::MateReferenceName),
            mate_reference_start: options.wants_column(BamColumn::MatePosition),
            template_length: options.wants_column(BamColumn::TemplateLength),
            query_sequence: options.wants_column(BamColumn::Sequence),
            query_qualities: options.wants_column(BamColumn::Quality),
            tags: !options.tags.is_empty(),
        }
    }
}

pub(crate) fn bam_sequence_to_optional_string(sequence: bam::record::Sequence<'_>) -> Option<String> {
    if sequence.is_empty() {
        None
    } else {
        Some(sequence.iter().map(|base| base as char).collect())
    }
}

pub(crate) fn bam_quality_to_optional_string(qualities: bam::record::QualityScores<'_>) -> Option<String> {
    if qualities.is_empty() {
        None
    } else {
        Some(String::from_utf8_lossy(qualities.as_ref()).into_owned())
    }
}

pub(crate) fn bam_cigar_to_string(cigar: bam::record::Cigar<'_>) -> String {
    cigar
        .iter()
        .filter_map(Result::ok)
        .map(|op| format!("{}{}", op.len(), kind_to_char(op.kind())))
        .collect()
}

fn cigar_to_string(cigar: &noodles::sam::alignment::record_buf::Cigar) -> String {
    cigar
        .as_ref()
        .iter()
        .map(|op| format!("{}{}", op.len(), kind_to_char(op.kind())))
        .collect()
}

fn kind_to_char(kind: Kind) -> char {
    match kind {
        Kind::Match => 'M',
        Kind::Insertion => 'I',
        Kind::Deletion => 'D',
        Kind::Skip => 'N',
        Kind::SoftClip => 'S',
        Kind::HardClip => 'H',
        Kind::Pad => 'P',
        Kind::SequenceMatch => '=',
        Kind::SequenceMismatch => 'X',
    }
}

fn bytes_to_optional_string(bytes: &[u8]) -> Option<String> {
    if bytes.is_empty() {
        None
    } else {
        Some(String::from_utf8_lossy(bytes).into_owned())
    }
}

fn parse_cigar(cigar: &str) -> Result<Cigar, NoodlesError> {
    let mut ops = Vec::new();
    let mut digits = String::new();

    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            digits.push(ch);
            continue;
        }

        let len = if digits.is_empty() {
            1
        } else {
            digits
                .parse::<usize>()
                .map_err(|err| NoodlesError::Message(format!("invalid CIGAR '{cigar}': {err}")))?
        };
        digits.clear();
        ops.push(Op::new(char_to_kind(ch)?, len));
    }

    if !digits.is_empty() {
        return Err(NoodlesError::Message(format!(
            "invalid CIGAR '{cigar}': trailing digits without operation"
        )));
    }

    Ok(ops.into_iter().collect())
}

fn char_to_kind(ch: char) -> Result<Kind, NoodlesError> {
    Ok(match ch {
        'M' => Kind::Match,
        'I' => Kind::Insertion,
        'D' => Kind::Deletion,
        'N' => Kind::Skip,
        'S' => Kind::SoftClip,
        'H' => Kind::HardClip,
        'P' => Kind::Pad,
        '=' => Kind::SequenceMatch,
        'X' => Kind::SequenceMismatch,
        other => {
            return Err(NoodlesError::Message(format!(
                "invalid CIGAR operation '{other}'"
            )));
        }
    })
}

fn tag_from_name(name: &str) -> Result<Tag, NoodlesError> {
    let bytes = name.as_bytes();
    if bytes.len() != 2 {
        return Err(NoodlesError::Message(format!(
            "auxiliary tag names must be two characters, got '{name}'"
        )));
    }
    Ok(Tag::from([bytes[0], bytes[1]]))
}

fn tag_value_to_noodles(value: &TagValue) -> Result<Value, NoodlesError> {
    Ok(match value {
        TagValue::Int(v) => Value::from(*v as i32),
        TagValue::Float(v) => Value::from(*v as f32),
        TagValue::String(v) => Value::from(v.as_str()),
        TagValue::Missing => Value::from(""),
    })
}

fn convert_tag_value(value: &Value) -> TagValue {
    match value {
        Value::Int8(v) => TagValue::Int(*v as i64),
        Value::UInt8(v) => TagValue::Int(*v as i64),
        Value::Int16(v) => TagValue::Int(*v as i64),
        Value::UInt16(v) => TagValue::Int(*v as i64),
        Value::Int32(v) => TagValue::Int(*v as i64),
        Value::UInt32(v) => TagValue::Int(*v as i64),
        Value::Float(v) => TagValue::Float(*v as f64),
        Value::String(v) => TagValue::String(v.to_string()),
        Value::Character(v) => TagValue::String((*v as char).to_string()),
        Value::Hex(v) => TagValue::String(v.to_string()),
        Value::Array(array) => TagValue::String(format!("{array:?}")),
    }
}