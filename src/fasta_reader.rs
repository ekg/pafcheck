use anyhow::{Result, Context};
use rust_htslib::bgzf::Reader;
use rust_htslib::faidx;

pub struct FastaReader {
    reader: Reader,
    index: faidx::Index,
}

impl FastaReader {
    pub fn new(fasta_path: &str) -> Result<Self> {
        let reader = Reader::from_path(fasta_path).context("Failed to open FASTA file")?;
        let index = faidx::Index::from_file(fasta_path).context("Failed to load FASTA index")?;
        Ok(FastaReader { reader, index })
    }

    pub fn fetch_sequence(
        &mut self,
        seq_name: &str,
        start: u64,
        end: u64,
    ) -> Result<String> {
        let tid = self.index.tid(seq_name).context("Failed to get sequence ID")?;
        let mut sequence = Vec::new();
        self.index.fetch_seq(&self.reader, tid, start, end, &mut sequence)
            .context("Failed to fetch sequence")?;
        Ok(String::from_utf8_lossy(&sequence).into_owned())
    }
}
