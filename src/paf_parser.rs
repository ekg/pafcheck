use anyhow::{Result, Context};

#[derive(Debug)]
pub struct PafRecord {
    pub query_name: String,
    pub query_length: u64,
    pub query_start: u64,
    pub query_end: u64,
    pub strand: char,
    pub target_name: String,
    pub target_length: u64,
    pub target_start: u64,
    pub target_end: u64,
    pub cigar: String,
}

impl PafRecord {
    pub fn from_line(line: &str) -> Result<Self> {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 {
            anyhow::bail!("PAF line does not have enough fields");
        }

        let mut cigar = String::new();
        for field in &fields[12..] {
            if field.starts_with("cg:Z:") {
                cigar = field[5..].to_string();
                break;
            }
        }

        Ok(PafRecord {
            query_name: fields[0].to_string(),
            query_length: fields[1].parse().context("Failed to parse query length")?,
            query_start: fields[2].parse().context("Failed to parse query start")?,
            query_end: fields[3].parse().context("Failed to parse query end")?,
            strand: fields[4].chars().next().unwrap(),
            target_name: fields[5].to_string(),
            target_length: fields[6].parse().context("Failed to parse target length")?,
            target_start: fields[7].parse().context("Failed to parse target start")?,
            target_end: fields[8].parse().context("Failed to parse target end")?,
            cigar,
        })
    }
}
