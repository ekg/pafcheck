use std::collections::HashSet;
use std::fs;

use super::agc_index::AgcIndex;
use super::faidx::FastaIndex;

pub fn collect_sequence_paths(
    mut files: Vec<String>,
    sequence_list: Option<String>,
) -> Result<Vec<String>, String> {
    if let Some(list_path) = sequence_list {
        let contents = fs::read_to_string(&list_path)
            .map_err(|e| format!("Failed to read --sequence-list '{list_path}': {e}"))?;
        for line in contents.lines() {
            let entry = line.trim();
            if entry.is_empty() || entry.starts_with('#') {
                continue;
            }
            files.push(entry.to_string());
        }
    }

    let mut seen = HashSet::new();
    files.retain(|path| seen.insert(path.clone()));
    Ok(files)
}

// Helper to extract file extension
fn get_extension(filename: &str) -> &str {
    if let Some(dot_pos) = filename.rfind('.') {
        &filename[dot_pos + 1..]
    } else {
        ""
    }
}

/// Unified sequence index that supports both FASTA and AGC formats
#[derive(Debug)]
pub enum SequenceIndex {
    Fasta(FastaIndex),
    Agc(AgcIndex),
}

impl SequenceIndex {
    pub fn build(files: &[String]) -> Result<Self, String> {
        if files.is_empty() {
            return Err("No sequence files provided".to_string());
        }

        // Check file extensions to determine type
        let first_ext = get_extension(&files[0]);
        let all_same_type = files.iter().all(|f| get_extension(f) == first_ext);

        if !all_same_type {
            return Err(
                "Mixed file types not supported. All files must be the same format (all .fa/.fasta/.fna or all .agc)"
                    .to_string(),
            );
        }

        match first_ext {
            "fa" | "fasta" | "fna" => {
                let index = FastaIndex::build(files)?;
                Ok(SequenceIndex::Fasta(index))
            }
            "agc" => {
                let index = AgcIndex::build(files)?;
                Ok(SequenceIndex::Agc(index))
            }
            _ => {
                // Default to FASTA for unknown extensions or no extension
                let index = FastaIndex::build(files)?;
                Ok(SequenceIndex::Fasta(index))
            }
        }
    }

    pub fn fetch_sequence(
        &self,
        seq_name: &str,
        start: usize,
        end: usize,
    ) -> Result<Vec<u8>, String> {
        match self {
            SequenceIndex::Fasta(index) => index.fetch_sequence(seq_name, start, end),
            SequenceIndex::Agc(index) => index.fetch_sequence(seq_name, start, end),
        }
    }

    /// Create a SequenceIndex from in-memory FASTA content (for testing)
    pub fn from_fasta_content(fasta_content: &str) -> Result<Self, String> {
        let index = FastaIndex::from_fasta_content(fasta_content)?;
        Ok(SequenceIndex::Fasta(index))
    }
}
