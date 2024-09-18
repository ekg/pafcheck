use clap::{App, Arg};
use anyhow::{Result, Context};
use std::fs::File;
use std::io::{BufRead, BufReader};

mod paf_parser;
mod fasta_reader;
mod cigar_parser;
mod validator;

use crate::fasta_reader::FastaReader;
use crate::paf_parser::PafRecord;
use crate::validator::validate_record;

fn main() -> Result<()> {
    let matches = App::new("PAF Validator")
        .version("1.0")
        .author("Your Name")
        .about("Validates PAF CIGAR strings against a FASTA file")
        .arg(
            Arg::with_name("fasta")
                .short('f')
                .long("fasta")
                .value_name("FASTA")
                .help("Path to the bgzip-compressed and tabix-indexed FASTA file")
                .takes_value(true)
                .required(true),
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

    let fasta_path = matches.value_of("fasta").unwrap();
    let paf_path = matches.value_of("paf").unwrap();
    let error_mode = matches.value_of("error-mode").unwrap();

    validate_paf(fasta_path, paf_path, error_mode)
        .context("Failed to validate PAF file")?;

    Ok(())
}

fn validate_paf(fasta_path: &str, paf_path: &str, error_mode: &str) -> Result<()> {
    let mut fasta_reader = FastaReader::new(fasta_path)
        .context("Failed to create FASTA reader")?;
    let paf_file = File::open(paf_path)
        .context("Failed to open PAF file")?;
    let reader = BufReader::new(paf_file);

    for (line_number, line) in reader.lines().enumerate() {
        let line = line.context("Failed to read PAF line")?;
        let record = PafRecord::from_line(&line)
            .context(format!("Failed to parse PAF record at line {}", line_number + 1))?;
        validate_record(&record, &mut fasta_reader, error_mode)
            .context(format!("Failed to validate PAF record at line {}", line_number + 1))?;
    }

    println!("PAF validation completed successfully.");
    Ok(())
}
