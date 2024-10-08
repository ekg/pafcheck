use crate::cigar_parser::{parse_cigar, CigarOp};
use crate::fasta_reader::MultiFastaReader;
use crate::paf_parser::PafRecord;
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::io::Write;
use thiserror::Error;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum ErrorType {
    Mismatch,
    LengthMismatch,
    CigarMismatch,
}

#[derive(Error, Debug)]
pub struct ValidationError {
    pub errors: HashMap<ErrorType, ErrorInfo>,
}

#[derive(Debug)]
pub struct ErrorInfo {
    pub first_message: String,
    pub count: usize,
}

impl std::fmt::Display for ValidationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Validation errors:")?;
        for (error_type, error_info) in &self.errors {
            writeln!(f, "  - {:?}: {}", error_type, error_info.first_message)?;
            if error_info.count > 1 {
                writeln!(
                    f,
                    "  - {:?}: Total occurrences: {}",
                    error_type, error_info.count
                )?;
            }
        }
        Ok(())
    }
}

pub fn validate_record<W: Write>(
    record: &PafRecord,
    fasta_reader: &mut MultiFastaReader,
    error_mode: &str,
    output: &mut W,
) -> Result<()> {
    let query_seq = fasta_reader
        .fetch_query_sequence(&record.query_name, record.query_start, record.query_end)
        .context(format!(
            "Failed to fetch query sequence: {} ({}:{})",
            record.query_name, record.query_start, record.query_end
        ))?;
    let target_seq = fasta_reader
        .fetch_target_sequence(&record.target_name, record.target_start, record.target_end)
        .context(format!(
            "Failed to fetch target sequence: {} ({}:{})",
            record.target_name, record.target_start, record.target_end
        ))?;

    let query_seq = if record.strand == '-' {
        reverse_complement(&query_seq)
    } else {
        query_seq
    };

    let query_seq = query_seq.to_uppercase().into_bytes();
    let target_seq = target_seq.to_uppercase().into_bytes();

    let cigar_ops = parse_cigar(&record.cigar).context("Failed to parse CIGAR string")?;

    let mut q_idx: usize = 0;
    let mut t_idx: usize = 0;
    let mut errors: HashMap<ErrorType, ErrorInfo> = HashMap::new();

    for (op_idx, op) in cigar_ops.iter().enumerate() {
        match op {
            CigarOp::Match(len) | CigarOp::Mismatch(len) => {
                let len = *len as usize;
                let q_slice = query_seq
                    .get(q_idx..q_idx + len)
                    .ok_or_else(|| anyhow::anyhow!("Query sequence index out of range"))?;
                let t_slice = target_seq
                    .get(t_idx..t_idx + len)
                    .ok_or_else(|| anyhow::anyhow!("Target sequence index out of range"))?;

                for i in 0..len {
                    let q = q_slice[i];
                    let t = t_slice[i];
                    let is_match = q == t;
                    let expected_match = matches!(op, CigarOp::Match(_));

                    if is_match != expected_match {
                        let error_type = if expected_match {
                            ErrorType::Mismatch
                        } else {
                            ErrorType::CigarMismatch
                        };

                        let error_message = format!(
                            "CIGAR mismatch at operation {}: query char '{}' at pos {} vs target char '{}' at pos {}",
                            op_idx, q as char, record.query_start + q_idx + i, t as char, record.target_start + t_idx + i
                        );

                        errors
                            .entry(error_type)
                            .and_modify(|e| e.count += 1)
                            .or_insert(ErrorInfo {
                                first_message: error_message,
                                count: 1,
                            });
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

    if q_idx != query_seq.len() {
        let error_type = ErrorType::LengthMismatch;
        let error_message = format!(
            "Query sequence length mismatch: CIGAR implies {}, actual length {}",
            q_idx,
            query_seq.len()
        );
        errors
            .entry(error_type)
            .and_modify(|e| e.count += 1)
            .or_insert(ErrorInfo {
                first_message: error_message,
                count: 1,
            });
    }
    if t_idx != target_seq.len() {
        let error_type = ErrorType::LengthMismatch;
        let error_message = format!(
            "Target sequence length mismatch: CIGAR implies {}, actual length {}",
            t_idx,
            target_seq.len()
        );
        errors
            .entry(error_type)
            .and_modify(|e| e.count += 1)
            .or_insert(ErrorInfo {
                first_message: error_message,
                count: 1,
            });
    }

    if !errors.is_empty() {
        if error_mode == "report" {
            for (error_type, error_info) in &errors {
                writeln!(output, "{:?}: {}", error_type, error_info.first_message)?;
                if error_info.count > 1 {
                    writeln!(
                        output,
                        "{:?}: Total occurrences: {}",
                        error_type, error_info.count
                    )?;
                }
            }
            Ok(())
        } else {
            Err(anyhow::anyhow!(ValidationError { errors }))
        }
    } else {
        Ok(())
    }
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::paf_parser::PafRecord;
    use std::io::BufWriter;

    #[test]
    fn test_false_mismatch_detection() {
        let query_fasta_content = ">query\nACTGACTGACTG";
        let target_fasta_content = ">target\nACTGACTGACTG";
        let cigar = "5=1X6=";

        let paf_record = PafRecord {
            query_name: "query".to_string(),
            query_length: 12,
            query_start: 0,
            query_end: 12,
            strand: '+',
            target_name: "target".to_string(),
            target_length: 12,
            target_start: 0,
            target_end: 12,
            cigar: cigar.to_string(),
        };

        let mut fasta_reader =
            MultiFastaReader::from_strings(query_fasta_content, target_fasta_content).unwrap();
        let mut output = BufWriter::new(Vec::new());

        let result = validate_record(&paf_record, &mut fasta_reader, "omit", &mut output);

        assert!(
            result.is_err(),
            "Expected an error due to CIGAR mismatch, but got success"
        );

        if let Err(e) = result {
            if let Some(validation_error) = e.downcast_ref::<ValidationError>() {
                assert!(
                    validation_error
                        .errors
                        .iter()
                        .any(|(error_type, error_info)| {
                            *error_type == ErrorType::CigarMismatch
                                && error_info
                                    .first_message
                                    .contains("CIGAR mismatch at operation 1")
                        }),
                    "Expected CIGAR mismatch error at operation 1, but got: {:?}",
                    validation_error.errors
                );
            } else {
                panic!("Expected ValidationError, but got: {}", e);
            }
        }
    }

    #[test]
    fn test_false_match_detection() {
        let query_fasta_content = ">query\nACTGACCGACTG";
        let target_fasta_content = ">target\nACTGACTGACTG";
        let cigar = "12=";

        let paf_record = PafRecord {
            query_name: "query".to_string(),
            query_length: 12,
            query_start: 0,
            query_end: 12,
            strand: '+',
            target_name: "target".to_string(),
            target_length: 12,
            target_start: 0,
            target_end: 12,
            cigar: cigar.to_string(),
        };

        let mut fasta_reader =
            MultiFastaReader::from_strings(query_fasta_content, target_fasta_content).unwrap();
        let mut output = BufWriter::new(Vec::new());

        let result = validate_record(&paf_record, &mut fasta_reader, "omit", &mut output);

        assert!(
            result.is_err(),
            "Expected an error due to sequence mismatch, but got success"
        );

        if let Err(e) = result {
            if let Some(validation_error) = e.downcast_ref::<ValidationError>() {
                assert!(
                    validation_error
                        .errors
                        .iter()
                        .any(|(error_type, error_info)| {
                            *error_type == ErrorType::Mismatch
                                && error_info
                                    .first_message
                                    .contains("CIGAR mismatch at operation 0")
                        }),
                    "Expected Mismatch error at operation 0, but got: {:?}",
                    validation_error.errors
                );
            } else {
                panic!("Expected ValidationError, but got: {}", e);
            }
        }
    }
}
