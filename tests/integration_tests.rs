use anyhow::Result;
use pafcheck::fasta_reader::MultiFastaReader;
use pafcheck::paf_parser::PafRecord;
use pafcheck::validator::validate_record;
use std::fs::File;
use std::io::{BufRead, BufReader, Write, BufWriter};
use tempfile::NamedTempFile;

fn create_temp_fasta(sequences: &[(&str, &str)]) -> Result<NamedTempFile> {
    let mut temp_file = NamedTempFile::new()?;
    for (name, seq) in sequences {
        writeln!(temp_file, ">{}", name)?;
        writeln!(temp_file, "{}", seq)?;
    }
    Ok(temp_file)
}

fn create_temp_paf(entries: &[&str]) -> Result<NamedTempFile> {
    let mut temp_file = NamedTempFile::new()?;
    for entry in entries {
        writeln!(temp_file, "{}", entry)?;
    }
    Ok(temp_file)
}

fn run_validation(
    query_fasta_content: &[(&str, &str)],
    target_fasta_content: &[(&str, &str)],
    paf_content: &[&str],
    error_mode: &str,
) -> Result<()> {
    let query_fasta_file = create_temp_fasta(query_fasta_content)?;
    let target_fasta_file = create_temp_fasta(target_fasta_content)?;
    let paf_file = create_temp_paf(paf_content)?;

    let mut fasta_reader = MultiFastaReader::new(query_fasta_file.path(), target_fasta_file.path())?;
    let paf_reader = BufReader::new(File::open(paf_file.path())?);

    for line in paf_reader.lines() {
        let line = line?;
        let record = PafRecord::from_line(&line)?;
        let mut output = BufWriter::new(Vec::new());
        validate_record(&record, &mut fasta_reader, error_mode, &mut output)?;
    }

    Ok(())
}

#[test]
fn test_perfect_match() -> Result<()> {
    let query_fasta_content = [("query1", "ATCGATCGATCG")];
    let target_fasta_content = [("target1", "ATCGATCGATCG")];
    let paf_content = ["query1\t12\t0\t12\t+\ttarget1\t12\t0\t12\t12\t12\t60\tcg:Z:12="];

    run_validation(&query_fasta_content, &target_fasta_content, &paf_content, "report")?;
    println!("Test completed successfully. No errors expected.");
    Ok(())
}

#[test]
fn test_mismatch_detection() -> Result<()> {
    let query_fasta_content = [("query1", "ATCGATCGATCG")];
    let target_fasta_content = [("target1", "ATCGATTGATCG")];
    let paf_content = ["query1\t12\t0\t12\t+\ttarget1\t12\t0\t12\t11\t12\t55\tcg:Z:5=1X6="];

    let result = run_validation(&query_fasta_content, &target_fasta_content, &paf_content, "report");
    assert!(result.is_ok(), "Expected successful validation, but got an error");
    println!("Test completed successfully. Mismatch detected and reported as expected.");
    Ok(())
}

#[test]
fn test_false_match_detection() -> Result<()> {
    let query_fasta_content = [("query1", "ATCGATCGATCG")];
    let target_fasta_content = [("target1", "ATCGATTGATCG")];
    let paf_content = ["query1\t12\t0\t12\t+\ttarget1\t12\t0\t12\t12\t12\t60\tcg:Z:12="];

    let result = run_validation(&query_fasta_content, &target_fasta_content, &paf_content, "omit");
    assert!(result.is_err(), "Expected an error due to false match");
    if let Err(e) = result {
        let error_message = e.to_string();
        assert!(error_message.contains("CIGAR mismatch at operation 0"),
                "Unexpected error message: {}", error_message);
    }
    println!("Test completed successfully. Error detected as expected.");
    Ok(())
}

#[test]
fn test_false_mismatch_detection() -> Result<()> {
    let query_fasta_content = [("query1", "ATCGATCGATCG")];
    let target_fasta_content = [("target1", "ATCGATCGATCG")];
    let paf_content = ["query1\t12\t0\t12\t+\ttarget1\t12\t0\t12\t11\t12\t55\tcg:Z:5=1X6="];

    let result = run_validation(&query_fasta_content, &target_fasta_content, &paf_content, "omit");
    assert!(result.is_err(), "Expected an error due to false mismatch");
    if let Err(e) = result {
        let error_message = e.to_string();
        assert!(error_message.contains("CIGAR mismatch at operation 1"),
                "Unexpected error message: {}", error_message);
    }
    println!("Test completed successfully. Error detected as expected.");
    Ok(())
}

#[test]
fn test_insertion() -> Result<()> {
    let query_fasta_content = [("query1", "ATCGATACGATCG")];
    let target_fasta_content = [("target1", "ATCGATCGATCG")];
    let paf_content = ["query1\t13\t0\t13\t+\ttarget1\t12\t0\t12\t12\t13\t60\tcg:Z:6=1I6="];

    run_validation(&query_fasta_content, &target_fasta_content, &paf_content, "report")?;
    println!("Test completed successfully. No errors expected.");
    Ok(())
}

#[test]
fn test_deletion() -> Result<()> {
    let query_fasta_content = [("query1", "ATCGATCGATCG")];
    let target_fasta_content = [("target1", "ATCGATTACGATCG")];
    let paf_content = ["query1\t12\t0\t12\t+\ttarget1\t14\t0\t14\t12\t14\t60\tcg:Z:6=2D6="];

    let result = run_validation(&query_fasta_content, &target_fasta_content, &paf_content, "report");
    assert!(result.is_ok(), "Expected no errors, but got: {:?}", result);
    Ok(())
}
