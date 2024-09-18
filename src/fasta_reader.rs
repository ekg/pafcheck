use anyhow::{Context, Result};
use rust_htslib::faidx;
use std::path::Path;

pub struct FastaReader {
    reader: faidx::Reader,
}

impl FastaReader {
    pub fn new<P: AsRef<Path>>(fasta_path: P) -> Result<Self> {
        let reader = faidx::Reader::from_path(fasta_path)
            .context("Failed to open FASTA file and load index")?;
        Ok(FastaReader { reader })
    }

    pub fn fetch_sequence(&self, seq_name: &str, start: usize, end: usize) -> Result<String> {
        self.reader
            .fetch_seq_string(seq_name, start, end - 1) // Adjust for 0-based indexing
            .context("Failed to fetch sequence")
    }
}
