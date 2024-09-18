use anyhow::{Result, Context};

#[derive(Debug)]
pub enum CigarOp {
    Match(u64),
    Mismatch(u64),
    Insertion(u64),
    Deletion(u64),
}

pub fn parse_cigar(cigar: &str) -> Result<Vec<CigarOp>> {
    let mut ops = Vec::new();
    let mut num = String::new();

    for c in cigar.chars() {
        if c.is_ascii_digit() {
            num.push(c);
        } else {
            let count = num.parse::<u64>().context("Failed to parse CIGAR operation count")?;
            match c {
                '=' => ops.push(CigarOp::Match(count)),
                'X' => ops.push(CigarOp::Mismatch(count)),
                'I' => ops.push(CigarOp::Insertion(count)),
                'D' => ops.push(CigarOp::Deletion(count)),
                _ => anyhow::bail!("Unknown CIGAR operation: {}", c),
            }
            num.clear();
        }
    }
    Ok(ops)
}
