use anyhow::{Context, Result};
use rust_htslib::faidx;
use std::path::Path;

pub struct MultiFastaReader {
    query_reader: faidx::Reader,
    target_reader: faidx::Reader,
}

impl MultiFastaReader {
    pub fn new<P: AsRef<Path>>(query_fasta: P, target_fasta: P) -> Result<Self> {
        let query_reader = faidx::Reader::from_path(&query_fasta)
            .context(format!("Failed to open query FASTA file: {:?}", query_fasta.as_ref()))?;
        let target_reader = faidx::Reader::from_path(&target_fasta)
            .context(format!("Failed to open target FASTA file: {:?}", target_fasta.as_ref()))?;
        Ok(MultiFastaReader {
            query_reader,
            target_reader,
        })
    }

    pub fn fetch_query_sequence(&self, seq_name: &str, start: usize, end: usize) -> Result<String> {
        self.fetch_sequence(&self.query_reader, seq_name, start, end)
            .context("Failed to fetch query sequence")
    }

    pub fn fetch_target_sequence(&self, seq_name: &str, start: usize, end: usize) -> Result<String> {
        self.fetch_sequence(&self.target_reader, seq_name, start, end)
            .context("Failed to fetch target sequence")
    }

    fn fetch_sequence(&self, reader: &faidx::Reader, seq_name: &str, start: usize, end: usize) -> Result<String> {
        reader
            .fetch_seq_string(seq_name, start, end - 1) // Adjust for 0-based indexing
            .context(format!("Failed to fetch sequence: {}", seq_name))
    }
}
