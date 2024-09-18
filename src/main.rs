use anyhow::{Context, Result};
use clap::{App, Arg};
use std::fs::File;
use std::io::{BufRead, BufReader};

use pafcheck::fasta_reader::MultiFastaReader;
use pafcheck::paf_parser::PafRecord;
use pafcheck::validator::validate_record;

fn main() -> Result<()> {
    let matches = App::new("PAF Validator")
        .version("1.0")
        .author("Your Name")
        .about("Validates PAF CIGAR strings against FASTA files")
        .arg(
            Arg::with_name("query_fasta")
                .short('q')
                .long("query-fasta")
                .value_name("QUERY_FASTA")
                .help("Path to the bgzip-compressed and tabix-indexed query FASTA file")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("target_fasta")
                .short('t')
                .long("target-fasta")
                .value_name("TARGET_FASTA")
                .help("Path to the bgzip-compressed and tabix-indexed target FASTA file")
                .takes_value(true)
                .required(false),
        )
        .arg(
            Arg::with_name("paf")
                .short('p')
                .long("paf")
                .value_name("PAF")
                .help("Path to the PAF file")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("error-mode")
                .short('e')
                .long("error-mode")
                .value_name("MODE")
                .help("Error handling mode: omit, report")
                .takes_value(true)
                .required(false)
                .default_value("omit"),
        )
        .get_matches();

    let query_fasta_path = matches.value_of("query_fasta").unwrap();
    let target_fasta_path = matches.value_of("target_fasta").unwrap_or(query_fasta_path);
    let paf_path = matches.value_of("paf").unwrap();
    let error_mode = matches.value_of("error-mode").unwrap();

    validate_paf(query_fasta_path, target_fasta_path, paf_path, error_mode)
        .context("Failed to validate PAF file")?;

    Ok(())
}

fn validate_paf(query_fasta: &str, target_fasta: &str, paf_path: &str, error_mode: &str) -> Result<()> {
    let mut fasta_reader = MultiFastaReader::new(query_fasta, target_fasta)
        .context("Failed to create FASTA readers")?;
    println!("Using query FASTA: {}", query_fasta);
    println!("Using target FASTA: {}", target_fasta);
    let paf_file = File::open(paf_path).context("Failed to open PAF file")?;
    let reader = BufReader::new(paf_file);

    for (line_number, line) in reader.lines().enumerate() {
        let line = line.context("Failed to read PAF line")?;
        let record = PafRecord::from_line(&line).context(format!(
            "Failed to parse PAF record at line {}",
            line_number + 1
        ))?;
        validate_record(&record, &mut fasta_reader, error_mode).context(format!(
            "Failed to validate PAF record at line {}",
            line_number + 1
        ))?;
    }

    println!("PAF validation completed successfully.");
    Ok(())
}
