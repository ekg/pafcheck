use anyhow::{Context, Result};
use rust_htslib::faidx;
use std::path::Path;
use std::collections::HashMap;

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

    pub fn from_strings(query_fasta: &str, target_fasta: &str) -> Result<Self> {
        let query_map = parse_fasta(query_fasta)?;
        let target_map = parse_fasta(target_fasta)?;

        let query_reader = create_in_memory_reader(&query_map)?;
        let target_reader = create_in_memory_reader(&target_map)?;

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

fn parse_fasta(fasta_content: &str) -> Result<HashMap<String, String>> {
    let mut sequences = HashMap::new();
    let mut current_seq_name = String::new();
    let mut current_seq = String::new();

    for line in fasta_content.lines() {
        if line.starts_with('>') {
            if !current_seq_name.is_empty() {
                sequences.insert(current_seq_name, current_seq);
                current_seq = String::new();
            }
            current_seq_name = line[1..].trim().to_string();
        } else {
            current_seq.push_str(line.trim());
        }
    }

    if !current_seq_name.is_empty() {
        sequences.insert(current_seq_name, current_seq);
    }

    Ok(sequences)
}

fn create_in_memory_reader(sequences: &HashMap<String, String>) -> Result<faidx::Reader> {
    let mut fasta_content = String::new();
    for (name, seq) in sequences {
        fasta_content.push('>');
        fasta_content.push_str(name);
        fasta_content.push('\n');
        fasta_content.push_str(seq);
        fasta_content.push('\n');
    }

    let temp_file = tempfile::NamedTempFile::new().context("Failed to create temporary file")?;
    std::fs::write(temp_file.path(), fasta_content).context("Failed to write temporary FASTA file")?;
    faidx::Reader::from_path(temp_file.path()).context("Failed to create FASTA reader from temporary file")
}
