use ragc_reader::{Decompressor, DecompressorConfig};
use rayon::prelude::*;
use std::collections::HashMap;
use std::fmt;
use std::sync::{Arc, Mutex};

/// Structure to manage AGC archives using ragc-reader
pub struct AgcIndex {
    decompressors: Arc<Mutex<Vec<Decompressor>>>,
    pub agc_paths: Vec<String>,
    sample_contig_to_agc: HashMap<String, usize>,
    // Precomputed mapping from contig name to (sample, full_contig_name, agc_idx)
    contig_to_sample_info: HashMap<String, (String, String, usize)>,
}

impl fmt::Debug for AgcIndex {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("AgcIndex")
            .field("agc_paths", &self.agc_paths)
            .field(
                "num_decompressors",
                &self.decompressors.lock().unwrap().len(),
            )
            .finish_non_exhaustive()
    }
}

impl AgcIndex {
    fn new() -> Self {
        AgcIndex {
            decompressors: Arc::new(Mutex::new(Vec::new())),
            agc_paths: Vec::new(),
            sample_contig_to_agc: HashMap::default(),
            contig_to_sample_info: HashMap::default(),
        }
    }

    fn extract_short_contig_name(full_name: &str) -> &str {
        full_name.split_whitespace().next().unwrap_or(full_name)
    }

    pub fn build(agc_files: &[String]) -> Result<Self, String> {
        let mut index = AgcIndex::new();

        // Parallel metadata extraction phase
        let metadata_results: Vec<_> =
            agc_files
                .par_iter()
                .enumerate()
                .map(
                    |(agc_idx, agc_path)| -> Result<
                        (usize, String, Decompressor, Vec<(String, Vec<String>)>),
                        String,
                    > {
                        let config = DecompressorConfig {
                            verbosity: 0,
                            max_segment_cache_entries: 1, // ~1MB cache (16 x 60KB segments)
                        };
                        let mut decompressor = Decompressor::open(agc_path, config)
                            .map_err(|e| format!("Failed to open AGC file: {agc_path}: {e}"))?;

                        // Get all samples and their contigs
                        let samples = decompressor.list_samples();
                        let sample_contigs: Vec<_> = samples
                            .into_iter()
                            .map(|sample| {
                                let contigs =
                                    decompressor.list_contigs(&sample).unwrap_or_default();
                                (sample, contigs)
                            })
                            .collect();

                        Ok((agc_idx, agc_path.clone(), decompressor, sample_contigs))
                    },
                )
                .collect::<Result<Vec<_>, String>>()?;

        // Sequential assembly phase to maintain order and avoid shared mutable state issues
        for (agc_idx, agc_path, decompressor, sample_contigs) in metadata_results {
            index.agc_paths.push(agc_path);

            for (sample, contigs) in sample_contigs {
                for contig in contigs {
                    // Create a key that combines full contig name and sample name
                    let key = format!("{contig}@{sample}");
                    index.sample_contig_to_agc.insert(key, agc_idx);

                    // Also insert just the full contig name if it's unique
                    index
                        .sample_contig_to_agc
                        .entry(contig.clone())
                        .or_insert(agc_idx);

                    // Precompute contig-to-sample mappings for fast lookup
                    let sample_info = (sample.clone(), contig.clone(), agc_idx);

                    // Map full contig name to sample info
                    index
                        .contig_to_sample_info
                        .entry(contig.clone())
                        .or_insert(sample_info.clone());

                    // Extract short contig name and create mappings
                    let short_contig = Self::extract_short_contig_name(&contig);

                    // If short name differs from full name, also create mappings for short name
                    if short_contig != contig {
                        // Create key with short contig name and sample
                        let short_key = format!("{short_contig}@{sample}");
                        index
                            .sample_contig_to_agc
                            .entry(short_key)
                            .or_insert(agc_idx);

                        // Also insert just the short contig name if it's unique
                        index
                            .sample_contig_to_agc
                            .entry(short_contig.to_string())
                            .or_insert(agc_idx);

                        // Map short contig name to sample info
                        index
                            .contig_to_sample_info
                            .entry(short_contig.to_string())
                            .or_insert(sample_info);
                    }
                }
            }

            index.decompressors.lock().unwrap().push(decompressor);
        }

        // Shrink to fit after building
        index.sample_contig_to_agc.shrink_to_fit();
        index.contig_to_sample_info.shrink_to_fit();

        Ok(index)
    }

    fn parse_query(&self, seq_name: &str) -> (String, String, Option<usize>) {
        // Parse queries in the format:
        // - "contig@sample" -> (sample, contig, agc_idx)
        // - "contig" -> (sample, contig, agc_idx) if contig is unique

        if let Some((contig, sample)) = seq_name.split_once('@') {
            let agc_idx = self.sample_contig_to_agc.get(seq_name).copied();
            (sample.to_string(), contig.to_string(), agc_idx)
        } else if let Some((sample, full_contig, agc_idx)) =
            self.contig_to_sample_info.get(seq_name)
        {
            (sample.clone(), full_contig.clone(), Some(*agc_idx))
        } else {
            (String::new(), seq_name.to_string(), None)
        }
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

        let (sample, contig, agc_idx) = self.parse_query(seq_name);

        let agc_idx =
            agc_idx.ok_or_else(|| format!("Sequence '{seq_name}' not found in any AGC file"))?;

        // ragc uses 0-based coordinates with exclusive end
        let mut decompressors = self.decompressors.lock().unwrap();
        let sequence = decompressors[agc_idx]
            .get_contig_range(&sample, &contig, start, end)
            .map_err(|e| {
                format!("Failed to fetch sequence '{contig}@{sample}:{start}:{end}': {e}")
            })?;

        // Convert from numeric encoding (0-3) to ASCII (A,C,G,T)
        Ok(sequence
            .into_iter()
            .map(|b| match b {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                3 => b'T',
                _ => b'N',
            })
            .collect())
    }
}
