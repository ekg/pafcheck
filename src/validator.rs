use crate::cigar_parser::{parse_cigar, CigarOp};
use crate::fasta_reader::FastaReader;
use crate::paf_parser::PafRecord;
use anyhow::{Context, Result};

pub fn validate_record(
    record: &PafRecord,
    fasta_reader: &mut FastaReader,
    error_mode: &str,
) -> Result<()> {
    let query_seq = fasta_reader
        .fetch_sequence(&record.query_name, record.query_start, record.query_end)
        .context("Failed to fetch query sequence")?;
    let target_seq = fasta_reader
        .fetch_sequence(&record.target_name, record.target_start, record.target_end)
        .context("Failed to fetch target sequence")?;

    // Debug print
    println!("Query sequence length: {}", query_seq.len());
    println!("Target sequence length: {}", target_seq.len());

    let cigar_ops = parse_cigar(&record.cigar).context("Failed to parse CIGAR string")?;

    let mut q_idx: usize = 0;
    let mut t_idx: usize = 0;
    for op in cigar_ops {
        match op {
            CigarOp::Match(len) => {
                let q_slice = &query_seq[q_idx..q_idx + len as usize];
                let t_slice = &target_seq[t_idx..t_idx + len as usize];
                if q_slice != t_slice {
                    report_error(
                        error_mode,
                        &format!(
                            "Mismatch in Match operation at position {}: query {} vs target {}",
                            record.query_start + q_idx,
                            q_slice,
                            t_slice
                        ),
                        record,
                    )?;
                }
                q_idx += len as usize;
                t_idx += len as usize;
            }
            CigarOp::Mismatch(len) => {
                let q_slice = &query_seq[q_idx..q_idx + len as usize];
                let t_slice = &target_seq[t_idx..t_idx + len as usize];
                if q_slice == t_slice {
                    report_error(
                        error_mode,
                        &format!(
                            "Match in Mismatch operation at position {}: query {} vs target {}",
                            record.query_start + q_idx,
                            q_slice,
                            t_slice
                        ),
                        record,
                    )?;
                }
                q_idx += len as usize;
                t_idx += len as usize;
            }
            CigarOp::Insertion(len) => {
                q_idx += len as usize;
            }
            CigarOp::Deletion(len) => {
                t_idx += len as usize;
            }
        }
    }

    // Verify endpoints
    if q_idx != (record.query_end - record.query_start) {
        let _ = report_error(
            error_mode,
            "Query sequence length mismatch after CIGAR operations",
            record,
        );
    }
    if t_idx != (record.target_end - record.target_start) {
        let _ = report_error(
            error_mode,
            "Target sequence length mismatch after CIGAR operations",
            record,
        );
    }

    Ok(())
}

fn report_error(error_mode: &str, message: &str, record: &PafRecord) -> Result<()> {
    if error_mode == "report" {
        println!("[PAF_CHECK] {}: {:?}", message, record);
        anyhow::bail!("{}", message);
    }
    Ok(())
}
