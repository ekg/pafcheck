use crate::cigar_parser::{parse_cigar, CigarOp};
use crate::fasta_reader::MultiFastaReader;
use crate::paf_parser::PafRecord;
use anyhow::{Context, Result};
use thiserror::Error;

#[derive(Debug, Clone)]
pub enum ErrorType {
    Mismatch,
    LengthMismatch,
    // Add other error types as needed
}

#[derive(Error, Debug)]
pub struct ValidationError {
    pub errors: Vec<(ErrorType, String)>,
}

impl std::fmt::Display for ValidationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Validation errors: {}", self.errors.iter().map(|(_, msg)| msg).collect::<Vec<_>>().join("; "))
    }
}

pub fn validate_record(
    record: &PafRecord,
    fasta_reader: &mut MultiFastaReader,
    error_mode: &str,
) -> Result<()> {
    // Fetch sequences
    let mut query_seq = fasta_reader
        .fetch_query_sequence(&record.query_name, record.query_start, record.query_end)
        .context(format!(
            "Failed to fetch query sequence: {} ({}:{})",
            record.query_name, record.query_start, record.query_end
        ))?;
    let mut target_seq = fasta_reader
        .fetch_target_sequence(&record.target_name, record.target_start, record.target_end)
        .context(format!(
            "Failed to fetch target sequence: {} ({}:{})",
            record.target_name, record.target_start, record.target_end
        ))?;

    // Handle reverse strand
    if record.strand == '-' {
        query_seq = reverse_complement(&query_seq);
    }

    // Convert sequences to uppercase and to byte arrays
    let query_seq = query_seq.to_uppercase().into_bytes();
    let target_seq = target_seq.to_uppercase().into_bytes();

    // Debug statements
    println!(
        "Expected query length: {}, fetched length: {}",
        record.query_end - record.query_start,
        query_seq.len()
    );
    println!(
        "Expected target length: {}, fetched length: {}",
        record.target_end - record.target_start,
        target_seq.len()
    );

    // Parse CIGAR string
    let cigar_ops = parse_cigar(&record.cigar).context("Failed to parse CIGAR string")?;

    let mut q_idx: usize = 0;
    let mut t_idx: usize = 0;
    let mut errors: Vec<(ErrorType, String)> = Vec::new();

    for (op_idx, op) in cigar_ops.iter().enumerate() {
        match op {
            CigarOp::Match(len) => {
                let len = *len as usize;
                let q_slice = &query_seq[q_idx..q_idx + len];
                let t_slice = &target_seq[t_idx..t_idx + len];
                for i in 0..len {
                    let q = q_slice[i];
                    let t = t_slice[i];
                    if q != t {
                        let error_message = format!(
                            "Mismatch in Match operation at CIGAR op {}, position {}: query {} vs target {}",
                            op_idx,
                            record.query_start + q_idx + i,
                            q as char,
                            t as char
                        );
                        errors.push((ErrorType::Mismatch, error_message));
                    }
                }
                q_idx += len;
                t_idx += len;
            }
            CigarOp::Mismatch(len) => {
                let len = *len as usize;
                let q_slice = &query_seq[q_idx..q_idx + len];
                let t_slice = &target_seq[t_idx..t_idx + len];
                for i in 0..len {
                    let q = q_slice[i];
                    let t = t_slice[i];
                    if q == t {
                        let error_message = format!(
                            "Match in Mismatch operation at CIGAR op {}, position {}: query {} vs target {}",
                            op_idx,
                            record.query_start + q_idx + i,
                            q as char,
                            t as char
                        );
                        errors.push((ErrorType::Mismatch, error_message));
                    }
                }
                q_idx += len;
                t_idx += len;
            }
            CigarOp::Insertion(len) => {
                q_idx += *len as usize;
            }
            CigarOp::Deletion(len) => {
                t_idx += *len as usize;
            }
        }
    }

    // Verify endpoints
    if q_idx != (record.query_end - record.query_start) as usize {
        errors.push((
            ErrorType::LengthMismatch,
            format!(
                "Query sequence length mismatch after CIGAR operations: expected {}, got {}",
                record.query_end - record.query_start,
                q_idx
            )
        ));
    }
    if t_idx != (record.target_end - record.target_start) as usize {
        errors.push((
            ErrorType::LengthMismatch,
            format!(
                "Target sequence length mismatch after CIGAR operations: expected {}, got {}",
                record.target_end - record.target_start,
                t_idx
            )
        ));
    }

    if !errors.is_empty() {
        for (_, error) in &errors {
            println!("Error: {}", error);
        }
        if error_mode == "report" {
            let messages: Vec<String> = errors.iter().map(|(_, msg)| msg.clone()).collect();
            anyhow::bail!("{}", messages.join("\n"));
        } else {
            return Err(ValidationError { errors }.into());
        }
    }

    Ok(())
}

fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|base| match base {
            'A' | 'a' => 'T',
            'T' | 't' => 'A',
            'G' | 'g' => 'C',
            'C' | 'c' => 'G',
            'N' | 'n' => 'N',
            _ => 'N',
        })
        .collect()
}
