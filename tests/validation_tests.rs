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
        assert!(result.is_ok(), "Expected validation to succeed in report mode");
        writer.flush()?;
    }

    // Convert captured output to string
    let output_str = String::from_utf8_lossy(&output);
    assert!(
        output_str.is_empty(),
        "Expected no errors, but got: {}",
        output_str
    );

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

    let mut output = Vec::new();
    let result = {
        let mut writer = BufWriter::new(&mut output);
        let result = validate_record(&paf_record, &mut fasta_reader, "report", &mut writer);
        writer.flush()?;
        result
    };

    assert!(result.is_ok(), "Expected validation to succeed in report mode");

    let output_str = String::from_utf8_lossy(&output);
    assert!(
        output_str.contains("Mismatch: CIGAR mismatch at operation 0: query char 'C' at pos 4 vs target char 'T' at pos 4"),
        "Expected mismatch was not reported in the output: {}",
        output_str
    );

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

    let mut output = Vec::new();
    let result = {
        let mut writer = BufWriter::new(&mut output);
        let result = validate_record(&paf_record, &mut fasta_reader, "report", &mut writer);
        writer.flush()?;
        result
    };
    
    assert!(result.is_ok(), "Expected validation to succeed in report mode");

    let output_str = String::from_utf8_lossy(&output);
    assert!(
        output_str.contains("CigarMismatch: CIGAR mismatch at operation 1: query char 'T' at pos 4 vs target char 'T' at pos 4"),
        "Expected mismatch was not reported in the output: {}",
        output_str
    );

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
        output_str.is_empty(),
        "Expected no errors, but got: {}",
        output_str
    );

    Ok(())
}
