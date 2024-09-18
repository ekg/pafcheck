use anyhow::{Result, Context};
use rust_htslib::bgzf::Reader;
use rust_htslib::tbx::Index;

pub struct FastaReader {
    reader: Reader,
    index: Index,
}

impl FastaReader {
    pub fn new(fasta_path: &str) -> Result<Self> {
        let reader = Reader::from_path(fasta_path).context("Failed to open FASTA file")?;
        let index = Index::load(fasta_path).context("Failed to load FASTA index")?;
        Ok(FastaReader { reader, index })
    }

    pub fn fetch_sequence(
        &mut self,
        seq_name: &str,
        start: u64,
        end: u64,
    ) -> Result<String> {
        let tid = self.index.tid(seq_name).context("Failed to get sequence ID")?;
        let intervals = self.index.fetch(tid, start, end).context("Failed to fetch sequence")?;

        let mut sequence = String::new();
        for interval in intervals {
            let data = self.reader.read_interval(interval).context("Failed to read sequence interval")?;
            sequence.push_str(&String::from_utf8_lossy(&data));
        }
        Ok(sequence)
    }
}
