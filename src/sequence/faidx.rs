use rust_htslib::faidx::Reader as FastaReader;
use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;

#[derive(Debug)]
pub struct FastaIndex {
    fasta_paths: Vec<PathBuf>,
    sequence_to_file: HashMap<String, usize>,
}

impl FastaIndex {
    pub fn build(fasta_files: &[String]) -> Result<Self, String> {
        let mut fasta_paths = Vec::new();
        let mut sequence_to_file = HashMap::new();

        for (idx, entry) in fasta_files.iter().enumerate() {
            let path = PathBuf::from(entry);
            if !path.exists() {
                return Err(format!("FASTA file '{entry}' not found"));
            }

            let path_str = path
                .to_str()
                .ok_or_else(|| format!("FASTA path contains invalid UTF-8: {entry}"))?
                .to_string();
            let fai_path = PathBuf::from(format!("{path_str}.fai"));

            if !fai_path.exists() {
                if let Err(e) = FastaReader::from_path(&path) {
                    return Err(format!("Failed to prepare FASTA index for '{entry}': {e}"));
                }
            }

            let fai_content = fs::read_to_string(&fai_path).map_err(|e| {
                format!("Failed to read FASTA index '{}': {}", fai_path.display(), e)
            })?;

            for line in fai_content.lines() {
                let fields: Vec<&str> = line.split('\t').collect();
                if let Some(seq_name) = fields.first().map(|s| s.trim()).filter(|s| !s.is_empty()) {
                    sequence_to_file.entry(seq_name.to_string()).or_insert(idx);
                }
            }

            fasta_paths.push(path);
        }

        Ok(Self {
            fasta_paths,
            sequence_to_file,
        })
    }

    pub fn fetch_sequence(
        &self,
        seq_name: &str,
        start: usize,
        end: usize,
    ) -> Result<Vec<u8>, String> {
        if start >= end {
            return Ok(Vec::new());
        }

        let fasta_idx = self
            .sequence_to_file
            .get(seq_name)
            .ok_or_else(|| format!("Sequence '{seq_name}' not found in supplied FASTA files"))?;

        let fasta_path = &self.fasta_paths[*fasta_idx];
        let reader = FastaReader::from_path(fasta_path)
            .map_err(|e| format!("Failed to open FASTA '{}': {}", fasta_path.display(), e))?;

        let raw_seq = reader.fetch_seq(seq_name, start, end - 1).map_err(|e| {
            format!(
                "Failed to fetch {seq_name}:{start}-{end} from '{}': {e}",
                fasta_path.display()
            )
        })?;
        let mut seq_vec = raw_seq.to_vec();
        // Apply the fix for the rust_htslib memory leak bug
        // https://github.com/rust-bio/rust-htslib/issues/401#issuecomment-1704290171
        unsafe { libc::free(raw_seq.as_ptr() as *mut std::ffi::c_void) };
        seq_vec.iter_mut().for_each(|b| *b = b.to_ascii_uppercase());
        Ok(seq_vec)
    }

    /// Create a FastaIndex from in-memory FASTA content (for testing)
    pub fn from_fasta_content(fasta_content: &str) -> Result<Self, String> {
        let temp_file = tempfile::NamedTempFile::new()
            .map_err(|e| format!("Failed to create temp file: {e}"))?;
        std::fs::write(temp_file.path(), fasta_content)
            .map_err(|e| format!("Failed to write temp FASTA: {e}"))?;

        let path_str = temp_file
            .path()
            .to_str()
            .ok_or("Temp path not valid UTF-8")?
            .to_string();

        // Create index
        FastaReader::from_path(temp_file.path())
            .map_err(|e| format!("Failed to index temp FASTA: {e}"))?;

        // We need to keep the temp file alive, so leak it
        let path = temp_file.into_temp_path();
        let path_buf = path.to_path_buf();
        std::mem::forget(path); // Leak to keep file alive

        let fai_path = PathBuf::from(format!("{path_str}.fai"));
        let fai_content = fs::read_to_string(&fai_path)
            .map_err(|e| format!("Failed to read FASTA index: {e}"))?;

        let mut sequence_to_file = HashMap::new();
        for line in fai_content.lines() {
            let fields: Vec<&str> = line.split('\t').collect();
            if let Some(seq_name) = fields.first().map(|s| s.trim()).filter(|s| !s.is_empty()) {
                sequence_to_file.insert(seq_name.to_string(), 0);
            }
        }

        Ok(Self {
            fasta_paths: vec![path_buf],
            sequence_to_file,
        })
    }
}
