use anyhow::Result;
use pafcheck::fasta_reader::MultiFastaReader;
use pafcheck::paf_parser::PafRecord;
use pafcheck::validator::validate_record;
use std::io::{Write, BufWriter};
use tempfile::NamedTempFile;

fn create_temp_fasta(sequences: &[(&str, &str)]) -> Result<NamedTempFile> {
    let mut temp_file = NamedTempFile::new()?;
    for (name, seq) in sequences {
        writeln!(temp_file, ">{}", name)?;
        writeln!(temp_file, "{}", seq)?;
    }
    Ok(temp_file)
}

#[test]
fn test_mismatch_detection() -> Result<()> {
    // Create temporary FASTA files
    let query_fasta = create_temp_fasta(&[("query1", "ACGTC")])?;
    let target_fasta = create_temp_fasta(&[("target1", "ACGTT")])?;

    // Create a PAF record with an intentional mismatch
    let paf_record = PafRecord {
        query_name: "query1".to_string(),
        query_length: 5,
        query_start: 0,
        query_end: 5,
        strand: '+',
        target_name: "target1".to_string(),
        target_length: 5,
        target_start: 0,
        target_end: 5,
        cigar: "4=1X".to_string(),
    };

    // Create MultiFastaReader
    let mut fasta_reader = MultiFastaReader::new(query_fasta.path(), target_fasta.path())?;

    // Capture output
    let mut output = Vec::new();
    {
        let mut writer = BufWriter::new(&mut output);
        let result = validate_record(&paf_record, &mut fasta_reader, "report", &mut writer);
        assert!(result.is_ok(), "Expected validation to succeed, but it failed");
        writer.flush()?;
    }

    // Convert captured output to string
    let output_str = String::from_utf8_lossy(&output);
    assert!(output_str.contains("Mismatch in Match operation at CIGAR op 0, position 4: query C vs target T"),
        "Expected mismatch was not reported in the output");

    Ok(())
}

#[test]
fn test_false_match_detection() -> Result<()> {
    let query_fasta = create_temp_fasta(&[("query1", "ACGTC")])?;
    let target_fasta = create_temp_fasta(&[("target1", "ACGTT")])?;

    let paf_record = PafRecord {
        query_name: "query1".to_string(),
        query_length: 5,
        query_start: 0,
        query_end: 5,
        strand: '+',
        target_name: "target1".to_string(),
        target_length: 5,
        target_start: 0,
        target_end: 5,
        cigar: "5=".to_string(),
    };

    let mut fasta_reader = MultiFastaReader::new(query_fasta.path(), target_fasta.path())?;

    let mut output = BufWriter::new(Vec::new());
    let result = validate_record(&paf_record, &mut fasta_reader, "report", &mut output);

    assert!(result.is_err(), "Expected validation to fail, but it succeeded");

    if let Err(e) = result {
        let error_message = e.to_string();
        assert!(
            error_message.contains("Mismatch in Match operation at CIGAR op 0, position 4: query C vs target T"),
            "Unexpected error message: {}",
            error_message
        );
    }

    Ok(())
}

#[test]
fn test_false_mismatch_detection() -> Result<()> {
    let query_fasta = create_temp_fasta(&[("query1", "ACGTT")])?;
    let target_fasta = create_temp_fasta(&[("target1", "ACGTT")])?;

    let paf_record = PafRecord {
        query_name: "query1".to_string(),
        query_length: 5,
        query_start: 0,
        query_end: 5,
        strand: '+',
        target_name: "target1".to_string(),
        target_length: 5,
        target_start: 0,
        target_end: 5,
        cigar: "4=1X".to_string(),
    };

    let mut fasta_reader = MultiFastaReader::new(query_fasta.path(), target_fasta.path())?;

    let mut output = BufWriter::new(Vec::new());
    let result = validate_record(&paf_record, &mut fasta_reader, "report", &mut output);

    assert!(result.is_err(), "Expected validation to fail, but it succeeded");

    if let Err(e) = result {
        let error_message = e.to_string();
        assert!(
            error_message.contains("Match in Mismatch operation at CIGAR op 1, position 4: query T vs target T"),
            "Unexpected error message: {}",
            error_message
        );
    }

    Ok(())
}
#[test]
fn test_mixed_match_mismatch_errors() -> Result<()> {
    let query_fasta = create_temp_fasta(&[("query1", "ACGTACGT")])?;
    let target_fasta = create_temp_fasta(&[("target1", "ACGTTCGT")])?;

    let paf_record = PafRecord {
        query_name: "query1".to_string(),
        query_length: 8,
        query_start: 0,
        query_end: 8,
        strand: '+',
        target_name: "target1".to_string(),
        target_length: 8,
        target_start: 0,
        target_end: 8,
        cigar: "4=1X3=".to_string(),
    };

    let mut fasta_reader = MultiFastaReader::new(query_fasta.path(), target_fasta.path())?;

    // Capture output
    let mut output = Vec::new();
    {
        let mut writer = BufWriter::new(&mut output);
        let result = validate_record(&paf_record, &mut fasta_reader, "report", &mut writer);
        assert!(result.is_ok(), "Expected validation to succeed in report mode");
        writer.flush()?;
    }

    // Convert captured output to string
    let output_str = String::from_utf8_lossy(&output);
    assert!(
        output_str.contains("Mismatch in Match operation at CIGAR op 0, position 4: query A vs target T") &&
        output_str.contains("Match in Mismatch operation at CIGAR op 1, position 5: query C vs target C"),
        "Expected mismatches were not reported in the output: {}",
        output_str
    );

    Ok(())
}
