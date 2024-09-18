use anyhow::Result;
use pafcheck::fasta_reader::FastaReader;
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

fn run_validation(fasta_content: &[(&str, &str)], paf_content: &[&str], error_mode: &str) -> Result<()> {
    let fasta_file = create_temp_fasta(fasta_content)?;
    let paf_file = create_temp_paf(paf_content)?;

    let mut fasta_reader = FastaReader::new(fasta_file.path())?;
    let paf_reader = BufReader::new(File::open(paf_file.path())?);

    for line in paf_reader.lines() {
        let line = line?;
        let record = PafRecord::from_line(&line)?;
        validate_record(&record, &mut fasta_reader, error_mode)?;
    }

    Ok(())
}

#[test]
fn test_perfect_match() -> Result<()> {
    let fasta_content: [(&str, &str); 2] = [
        ("query1", &"A".repeat(1000)),
        ("target1", &"A".repeat(1000)),
    ];
    let paf_content = ["query1\t1000\t100\t200\t+\ttarget1\t1000\t150\t250\t100\t100\t255\tcg:Z:100="];
    
    let errors = run_validation(&fasta_content, &paf_content, "report")?;
    assert!(errors.is_empty(), "Expected no errors, but got: {:?}", errors);
    Ok(())
}

#[test]
fn test_false_mismatch_detection() -> Result<()> {
    let fasta_content: [(&str, &str); 2] = [
        ("query1", &"C".repeat(1000)),
        ("target1", &"C".repeat(1000)),
    ];
    let paf_content = ["query1\t1000\t100\t200\t+\ttarget1\t1000\t150\t250\t100\t100\t255\tcg:Z:50=1X49="];
    
    let result = run_validation(&fasta_content, &paf_content, "report");
    assert!(result.is_err(), "Expected an error, but got: {:?}", result);
    let error = result.unwrap_err();
    assert!(error.to_string().contains("Match in Mismatch operation"), "Unexpected error message: {}", error);
    Ok(())
}

// Add more test functions for the remaining test cases...
