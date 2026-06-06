use crate::error::NoodlesError;
use noodles::bam as bam;
use noodles::bgzf as bgzf;
use noodles::core::region::Interval;
use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;
use noodles::sam::Header;
use noodles::sam::alignment::Record as _;

enum ChunkState {
    SeekNextChunk,
    ReadUntil(bgzf::VirtualPosition),
    Exhausted,
}

/// Lazy indexed region fetch that streams records chunk-by-chunk.
pub struct LazyIndexedFetch<R>
where
    R: bgzf::io::BufRead + bgzf::io::Seek,
{
    reader: bam::io::Reader<R>,
    header: Header,
    reference_sequence_id: usize,
    interval: Interval,
    chunks: std::vec::IntoIter<Chunk>,
    chunk_state: ChunkState,
    record: bam::Record,
}

impl<R> LazyIndexedFetch<R>
where
    R: bgzf::io::BufRead + bgzf::io::Seek,
{
    pub fn open(
        indexed: bam::io::IndexedReader<R>,
        header: Header,
        region: &noodles::core::Region,
    ) -> Result<Self, NoodlesError> {
        let reference_sequence_id = header
            .reference_sequences()
            .get_index_of(region.name())
            .ok_or_else(|| {
                NoodlesError::Message(format!(
                    "region reference sequence does not exist in header: {}",
                    region.name()
                ))
            })?;

        let chunks = indexed
            .index()
            .query(reference_sequence_id, region.interval())
            .map_err(NoodlesError::from)?;

        let reader = bam::io::Reader::from(indexed.into_inner());

        Ok(Self {
            reader,
            header,
            reference_sequence_id,
            interval: region.interval(),
            chunks: chunks.into_iter(),
            chunk_state: ChunkState::SeekNextChunk,
            record: bam::Record::default(),
        })
    }

    pub fn header(&self) -> &Header {
        &self.header
    }

    pub fn next_record(&mut self) -> Result<Option<bam::Record>, NoodlesError> {
        loop {
            match self.chunk_state {
                ChunkState::Exhausted => return Ok(None),
                ChunkState::SeekNextChunk => match self.chunks.next() {
                    Some(chunk) => {
                        self.reader
                            .get_mut()
                            .seek_to_virtual_position(chunk.start())
                            .map_err(NoodlesError::from)?;
                        self.chunk_state = ChunkState::ReadUntil(chunk.end());
                    }
                    None => {
                        self.chunk_state = ChunkState::Exhausted;
                        return Ok(None);
                    }
                },
                ChunkState::ReadUntil(chunk_end) => {
                    if self.reader.get_mut().virtual_position() >= chunk_end {
                        self.chunk_state = ChunkState::SeekNextChunk;
                        continue;
                    }

                    let block_size = self
                        .reader
                        .read_record(&mut self.record)
                        .map_err(NoodlesError::from)?;
                    if block_size == 0 {
                        self.chunk_state = ChunkState::SeekNextChunk;
                        continue;
                    }

                    if record_intersects_region(
                        &self.record,
                        self.reference_sequence_id,
                        self.interval,
                    )? {
                        return Ok(Some(std::mem::take(&mut self.record)));
                    }
                }
            }
        }
    }
}

fn record_intersects_region(
    record: &bam::Record,
    reference_sequence_id: usize,
    region_interval: Interval,
) -> Result<bool, NoodlesError> {
    match (
        record.reference_sequence_id().transpose().map_err(NoodlesError::from)?,
        record
            .alignment_start()
            .transpose()
            .map_err(NoodlesError::from)?,
        record
            .alignment_end()
            .transpose()
            .map_err(NoodlesError::from)?,
    ) {
        (Some(id), Some(start), Some(end)) => {
            let alignment_interval = (start..=end).into();
            Ok(id == reference_sequence_id && region_interval.intersects(alignment_interval))
        }
        _ => Ok(false),
    }
}