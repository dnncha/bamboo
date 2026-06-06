use bamboo_core::{BamColumn, BamScanOptions, BamTable, TagValue};
use noodles::sam::Header;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::Tag;
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
    pub fn from_record_buf(header: &Header, record: RecordBuf, tag_names: &[String]) -> Self {
        let reference_name = record
            .reference_sequence_id()
            .and_then(|id| header.reference_sequences().get_index(id))
            .map(|(name, _)| name.to_string());

        let mate_reference_name = record
            .mate_reference_sequence_id()
            .and_then(|id| header.reference_sequences().get_index(id))
            .map(|(name, _)| name.to_string());

        let tags = tag_names
            .iter()
            .map(|name| {
                let value = record
                    .data()
                    .get(&Tag::from([name.as_bytes()[0], name.as_bytes()[1]]))
                    .map(convert_tag_value)
                    .unwrap_or(TagValue::Missing);
                (name.clone(), value)
            })
            .collect();

        Self {
            query_name: record.name().map(|name| name.to_string()),
            flag: record.flags().bits(),
            reference_name,
            reference_start: record
                .alignment_start()
                .map(|position| position.get() as i64 - 1),
            mapping_quality: record.mapping_quality().map(|quality| quality.get()),
            cigar: cigar_to_string(record.cigar()),
            mate_reference_name,
            mate_reference_start: record
                .mate_alignment_start()
                .map(|position| position.get() as i64 - 1),
            template_length: Some(record.template_length()),
            query_sequence: bytes_to_optional_string(record.sequence().as_ref()),
            query_qualities: bytes_to_optional_string(record.quality_scores().as_ref()),
            tags,
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