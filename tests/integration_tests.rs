use anyhow::Result;
use pafcheck::fasta_reader::MultiFastaReader;
use pafcheck::paf_parser::PafRecord;
use pafcheck::validator::validate_record;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
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
) -> Result<Vec<String>> {
    let query_fasta_file = create_temp_fasta(query_fasta_content)?;
    let target_fasta_file = create_temp_fasta(target_fasta_content)?;
    let paf_file = create_temp_paf(paf_content)?;

    let mut fasta_reader = MultiFastaReader::new(query_fasta_file.path(), target_fasta_file.path())?;
    let paf_reader = BufReader::new(File::open(paf_file.path())?);

    let mut errors = Vec::new();

    for line in paf_reader.lines() {
        let line = line?;
        let record = PafRecord::from_line(&line)?;
        if let Err(e) = validate_record(&record, &mut fasta_reader, error_mode) {
            errors.push(e.to_string());
        }
    }

    Ok(errors)
}

#[test]
fn test_perfect_match() -> Result<()> {
    let query_fasta_content = [("query1", "ATCGATCGATCG")];
    let target_fasta_content = [("target1", "ATCGATCGATCG")];
    let paf_content = ["query1\t12\t0\t12\t+\ttarget1\t12\t0\t12\t12\t12\t60\tcg:Z:12="];

    let errors = run_validation(&query_fasta_content, &target_fasta_content, &paf_content, "report")?;
    assert!(errors.is_empty(), "Expected no errors, but got: {:?}", errors);
    Ok(())
}

#[test]
fn test_mismatch_detection() -> Result<()> {
    let query_fasta_content = [("query1", "ATCGATCGATCG")];
    let target_fasta_content = [("target1", "ATCGATTGATCG")];
    let paf_content = ["query1\t12\t0\t12\t+\ttarget1\t12\t0\t12\t11\t12\t55\tcg:Z:5=1X6="];

    let errors = run_validation(&query_fasta_content, &target_fasta_content, &paf_content, "report")?;
    assert!(errors.is_empty(), "Expected no errors, but got: {:?}", errors);
    println!("Test completed without errors. Check console output for any warnings.");
    Ok(())
}

#[test]
fn test_false_match_detection() -> Result<()> {
    let query_fasta_content = [("query1", "ATCGATCGATCG")];
    let target_fasta_content = [("target1", "ATCGATTGATCG")];
    let paf_content = ["query1\t12\t0\t12\t+\ttarget1\t12\t0\t12\t12\t12\t60\tcg:Z:12="];

    let errors = run_validation(&query_fasta_content, &target_fasta_content, &paf_content, "report")?;
    assert_eq!(errors.len(), 1, "Expected one error, but got: {:?}", errors);
    assert!(errors[0].contains("Mismatch in Match operation"), "Unexpected error message: {}", errors[0]);
    Ok(())
}

#[test]
fn test_false_mismatch_detection() -> Result<()> {
    let query_fasta_content = [("query1", "ATCGATCGATCG")];
    let target_fasta_content = [("target1", "ATCGATCGATCG")];
    let paf_content = ["query1\t12\t0\t12\t+\ttarget1\t12\t0\t12\t11\t12\t55\tcg:Z:5=1X6="];

    let errors = run_validation(&query_fasta_content, &target_fasta_content, &paf_content, "report")?;
    assert_eq!(errors.len(), 1, "Expected one error, but got: {:?}", errors);
    assert!(errors[0].contains("Match in Mismatch operation"), "Unexpected error message: {}", errors[0]);
    Ok(())
}

#[test]
fn test_insertion() -> Result<()> {
    let query_fasta_content = [("query1", "ATCGATACGATCG")];
    let target_fasta_content = [("target1", "ATCGATCGATCG")];
    let paf_content = ["query1\t13\t0\t13\t+\ttarget1\t12\t0\t12\t12\t13\t60\tcg:Z:6=1I6="];

    let errors = run_validation(&query_fasta_content, &target_fasta_content, &paf_content, "report")?;
    assert!(errors.is_empty(), "Expected no errors, but got: {:?}", errors);
    Ok(())
}

#[test]
fn test_deletion() -> Result<()> {
    let query_fasta_content = [("query1", "ATCGATCGATCG")];
    let target_fasta_content = [("target1", "ATCGATTACGATCG")];
    let paf_content = ["query1\t12\t0\t12\t+\ttarget1\t14\t0\t14\t12\t14\t60\tcg:Z:6=2D6="];

    let errors = run_validation(&query_fasta_content, &target_fasta_content, &paf_content, "report")?;
    assert!(errors.is_empty(), "Expected no errors, but got: {:?}", errors);
    Ok(())
}
