use crate::cigar_parser::{parse_cigar, CigarOp};
use crate::fasta_reader::MultiFastaReader;
use crate::paf_parser::PafRecord;
use anyhow::{Context, Result};

pub fn validate_record(
    record: &PafRecord,
    fasta_reader: &mut MultiFastaReader,
    error_mode: &str,
) -> Result<()> {
    let query_seq = fasta_reader
        .fetch_query_sequence(&record.query_name, record.query_start, record.query_end)
        .context(format!("Failed to fetch query sequence: {} ({}:{})", record.query_name, record.query_start, record.query_end))?;
    let target_seq = fasta_reader
        .fetch_target_sequence(&record.target_name, record.target_start, record.target_end)
        .context(format!("Failed to fetch target sequence: {} ({}:{})", record.target_name, record.target_start, record.target_end))?;

    println!("Query sequence: {} ({}:{}) length: {}", record.query_name, record.query_start, record.query_end, query_seq.len());
    println!("Target sequence: {} ({}:{}) length: {}", record.target_name, record.target_start, record.target_end, target_seq.len());

    let cigar_ops = parse_cigar(&record.cigar).context("Failed to parse CIGAR string")?;

    let mut q_idx: usize = 0;
    let mut t_idx: usize = 0;
    for (op_idx, op) in cigar_ops.iter().enumerate() {
        match op {
            CigarOp::Match(len) => {
                let q_slice = &query_seq[q_idx..q_idx + *len as usize];
                let t_slice = &target_seq[t_idx..t_idx + *len as usize];
                for (i, (q, t)) in q_slice.chars().zip(t_slice.chars()).enumerate() {
                    if q != t {
                        let error_message = format!(
                            "Mismatch in Match operation at CIGAR op {}, position {}: query {} vs target {}",
                            op_idx,
                            record.query_start + q_idx + i,
                            q,
                            t
                        );
                        println!("Warning: {}", error_message);
                        report_error(error_mode, &error_message, record)?;
                    }
                }
                q_idx += *len as usize;
                t_idx += *len as usize;
            }
            CigarOp::Mismatch(len) => {
                let q_slice = &query_seq[q_idx..q_idx + *len as usize];
                let t_slice = &target_seq[t_idx..t_idx + *len as usize];
                for (i, (q, t)) in q_slice.chars().zip(t_slice.chars()).enumerate() {
                    if q == t {
                        println!(
                            "Warning: Match in Mismatch operation at CIGAR op {}, position {}: query {} vs target {}",
                            op_idx,
                            record.query_start + q_idx + i,
                            q,
                            t
                        );
                    } else {
                        println!(
                            "Mismatch confirmed at CIGAR op {}, position {}: query {} vs target {}",
                            op_idx,
                            record.query_start + q_idx + i,
                            q,
                            t
                        );
                    }
                }
                q_idx += *len as usize;
                t_idx += *len as usize;
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
    if q_idx != (record.query_end - record.query_start) {
        report_error(
            error_mode,
            &format!(
                "Query sequence length mismatch after CIGAR operations: expected {}, got {}",
                record.query_end - record.query_start,
                q_idx
            ),
            record,
        )?;
    }
    if t_idx != (record.target_end - record.target_start) {
        report_error(
            error_mode,
            &format!(
                "Target sequence length mismatch after CIGAR operations: expected {}, got {}",
                record.target_end - record.target_start,
                t_idx
            ),
            record,
        )?;
    }

    Ok(())
}

fn report_error(error_mode: &str, message: &str, record: &PafRecord) -> Result<()> {
    println!("[PAF_CHECK] {}: {:?}", message, record);
    if error_mode == "report" {
        anyhow::bail!("{}", message);
    }
    Ok(())
}
