use anyhow::{Context, Result};
use clap::{App, Arg};
use indicatif::{ProgressBar, ProgressStyle};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::sync::Arc;
use std::thread;

use pafcheck::cigar_parser::{parse_cigar, CigarOp};
use pafcheck::fasta_reader::MultiFastaReader;
use pafcheck::paf_parser::PafRecord;
use pafcheck::validator::{validate_record, ErrorType, ValidationError};

fn main() {
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
                .required_unless_one(["info", "coverage"]),
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
        .arg(
            Arg::with_name("info")
                .long("info")
                .help("Show PAF statistics (mapping lengths, identities, CIGAR presence)")
                .takes_value(false)
                .required(false),
        )
        .arg(
            Arg::with_name("coverage")
                .long("coverage")
                .help("Show detailed coverage table for each query-target pair")
                .takes_value(false)
                .required(false),
        )
        .arg(
            Arg::with_name("threads")
                .short('j')
                .long("threads")
                .value_name("THREADS")
                .help("Number of threads to use for validation")
                .takes_value(true)
                .required(false),
        )
        .get_matches();

    let query_fasta_path = matches.value_of("query_fasta");
    let target_fasta_path = matches.value_of("target_fasta");
    let paf_path = matches.value_of("paf").unwrap();
    let error_mode = matches.value_of("error-mode").unwrap();

    let show_info = matches.is_present("info");
    let show_coverage = matches.is_present("coverage");

    if show_info {
        if let Err(e) = show_paf_info(paf_path) {
            eprintln!("[pafcheck] Error: {e}");
            std::process::exit(1);
        }
    } else if show_coverage {
        if let Err(e) = show_paf_coverage(paf_path) {
            eprintln!("[pafcheck] Error: {e}");
            std::process::exit(1);
        }
    } else {
        let num_threads = if let Some(threads_str) = matches.value_of("threads") {
            threads_str
                .parse::<usize>()
                .context("Invalid number of threads")
                .unwrap_or_else(|e| {
                    eprintln!("[pafcheck] Error: {e}");
                    std::process::exit(1);
                })
        } else if std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1)
            >= 4
        {
            4
        } else {
            1
        };

        let query_fasta: &str = query_fasta_path.unwrap();
        let target_fasta = target_fasta_path.unwrap_or(query_fasta);

        if let Err(e) = validate_paf(
            query_fasta,
            target_fasta,
            paf_path,
            error_mode,
            num_threads,
        ) {
            eprintln!("[pafcheck] Error: {e}");
            std::process::exit(1);
        }
    }
}

struct ThreadResult {
    error_count: usize,
    error_type_counts: HashMap<ErrorType, usize>,
    error_messages: Vec<String>,
}

#[derive(Debug)]
struct AlignmentPairStats {
    query_name: String,
    target_name: String,
    alignments: Vec<AlignmentInfo>,
    query_length: usize,
    target_length: usize,
}

#[derive(Debug, Clone)]
struct AlignmentInfo {
    query_start: usize,
    query_end: usize,
    target_start: usize,
    target_end: usize,
    matching_bases: usize,
    block_length: usize,
    cigar: String,
}

fn show_paf_coverage(paf_path: &str) -> Result<()> {
    let paf_file = File::open(paf_path).context("Failed to open PAF file")?;
    let reader = BufReader::new(paf_file);

    let mut pair_stats: HashMap<(String, String), AlignmentPairStats> = HashMap::new();

    for (line_number, line) in reader.lines().enumerate() {
        let line = line.context("Failed to read PAF line")?;
        let record = PafRecord::from_line(&line).context(format!(
            "Failed to parse PAF record at line {}",
            line_number + 1
        ))?;

        let key = (record.query_name.clone(), record.target_name.clone());

        let alignment_info = AlignmentInfo {
            query_start: record.query_start,
            query_end: record.query_end,
            target_start: record.target_start,
            target_end: record.target_end,
            matching_bases: record.matching_bases,
            block_length: record.block_length,
            cigar: record.cigar.clone(),
        };

        pair_stats
            .entry(key.clone())
            .and_modify(|stats| stats.alignments.push(alignment_info.clone()))
            .or_insert_with(|| AlignmentPairStats {
                query_name: record.query_name.clone(),
                target_name: record.target_name.clone(),
                alignments: vec![alignment_info],
                query_length: record.query_length,
                target_length: record.target_length,
            });
    }

    print_coverage_table(&pair_stats)?;
    Ok(())
}

fn print_coverage_table(pair_stats: &HashMap<(String, String), AlignmentPairStats>) -> Result<()> {
    // Print header
    println!(
        "{:<15} {:<15} {:<8} {:<10} {:<10} {:<10} {:<12} {:<12} {:<8} {:<12} {:<10} {:<12} {:<12}",
        "Query",
        "Target",
        "NumAlign",
        "MedianLen",
        "MeanLen",
        "StdLen",
        "QueryCov",
        "TargetCov",
        "Jaccard",
        "BPJaccard",
        "Identity",
        "GapCompId",
        "BlastId"
    );

    for ((_query, _target), stats) in pair_stats {
        let coverage_stats = calculate_pair_coverage_stats(stats)?;

        println!("{:<15} {:<15} {:<8} {:<10.1} {:<10.1} {:<10.1} {:<12.2} {:<12.2} {:<8.3} {:<12} {:<10.2} {:<12.2} {:<12.2}",
                 stats.query_name,
                 stats.target_name,
                 coverage_stats.num_alignments,
                 coverage_stats.median_length,
                 coverage_stats.mean_length,
                 coverage_stats.std_length,
                 coverage_stats.query_coverage * 100.0,
                 coverage_stats.target_coverage * 100.0,
                 coverage_stats.jaccard,
                 coverage_stats.bp_jaccard.map_or("NA".to_string(), |j| format!("{j:.3}")),
                 coverage_stats.identity * 100.0,
                 coverage_stats.gap_compressed_identity * 100.0,
                 coverage_stats.blast_identity * 100.0);
    }

    Ok(())
}

#[derive(Debug)]
struct PairCoverageStats {
    num_alignments: usize,
    median_length: f64,
    mean_length: f64,
    std_length: f64,
    query_coverage: f64,
    target_coverage: f64,
    jaccard: f64,
    bp_jaccard: Option<f64>,
    identity: f64,
    gap_compressed_identity: f64,
    blast_identity: f64,
}

fn calculate_pair_coverage_stats(stats: &AlignmentPairStats) -> Result<PairCoverageStats> {
    let num_alignments = stats.alignments.len();

    // Calculate length statistics
    let lengths: Vec<usize> = stats
        .alignments
        .iter()
        .map(|a| a.query_end - a.query_start)
        .collect();

    let mean_length = lengths.iter().sum::<usize>() as f64 / lengths.len() as f64;

    let mut sorted_lengths = lengths.clone();
    sorted_lengths.sort_unstable();
    let median_length = if sorted_lengths.len() % 2 == 0 {
        (sorted_lengths[sorted_lengths.len() / 2 - 1] + sorted_lengths[sorted_lengths.len() / 2])
            as f64
            / 2.0
    } else {
        sorted_lengths[sorted_lengths.len() / 2] as f64
    };

    let std_length = if lengths.len() > 1 {
        let variance = lengths
            .iter()
            .map(|&l| (l as f64 - mean_length).powi(2))
            .sum::<f64>()
            / (lengths.len() - 1) as f64;
        variance.sqrt()
    } else {
        0.0
    };

    // Calculate coverages
    let query_covered = calculate_covered_bases(
        &stats
            .alignments
            .iter()
            .map(|a| (a.query_start, a.query_end))
            .collect::<Vec<_>>(),
    );
    let target_covered = calculate_covered_bases(
        &stats
            .alignments
            .iter()
            .map(|a| (a.target_start, a.target_end))
            .collect::<Vec<_>>(),
    );

    let query_coverage = query_covered as f64 / stats.query_length as f64;
    let target_coverage = target_covered as f64 / stats.target_length as f64;

    // Calculate Jaccard index (intersection / union)
    // For simplicity, we approximate this as min(query_cov, target_cov) / max(query_cov, target_cov)
    let jaccard = if query_covered.max(target_covered) > 0 {
        query_covered.min(target_covered) as f64 / query_covered.max(target_covered) as f64
    } else {
        0.0
    };

    // Calculate base pair Jaccard using CIGAR strings
    let bp_jaccard = calculate_bp_jaccard(&stats.alignments)?;

    // Calculate identity metrics
    let total_matching: usize = stats.alignments.iter().map(|a| a.matching_bases).sum();
    let total_block_length: usize = stats.alignments.iter().map(|a| a.block_length).sum();
    let total_query_aligned: usize = stats
        .alignments
        .iter()
        .map(|a| a.query_end - a.query_start)
        .sum();

    let identity = total_matching as f64 / total_query_aligned as f64;
    let gap_compressed_identity = if total_block_length > 0 {
        total_matching as f64 / total_block_length as f64
    } else {
        0.0
    };

    // BLAST identity is typically matches / alignment_length including gaps
    let blast_identity = gap_compressed_identity; // Same as gap-compressed for now

    Ok(PairCoverageStats {
        num_alignments,
        median_length,
        mean_length,
        std_length,
        query_coverage,
        target_coverage,
        jaccard,
        bp_jaccard,
        identity,
        gap_compressed_identity,
        blast_identity,
    })
}

fn calculate_bp_jaccard(alignments: &[AlignmentInfo]) -> Result<Option<f64>> {
    let mut total_matches = 0u64;
    let mut total_mismatches = 0u64;
    let mut total_insertions = 0u64;
    let mut total_deletions = 0u64;

    let mut has_cigars = false;

    for alignment in alignments {
        if !alignment.cigar.is_empty() {
            has_cigars = true;
            let ops = parse_cigar(&alignment.cigar)?;

            for op in ops {
                match op {
                    CigarOp::Match(count) => total_matches += count,
                    CigarOp::Mismatch(count) => total_mismatches += count,
                    CigarOp::Insertion(count) => total_insertions += count,
                    CigarOp::Deletion(count) => total_deletions += count,
                }
            }
        }
    }

    if !has_cigars {
        return Ok(None);
    }

    // Base pair Jaccard: matches / (matches + mismatches + insertions + deletions)
    let total_operations = total_matches + total_mismatches + total_insertions + total_deletions;
    if total_operations > 0 {
        Ok(Some(total_matches as f64 / total_operations as f64))
    } else {
        Ok(Some(0.0))
    }
}

fn show_paf_info(paf_path: &str) -> Result<()> {
    let paf_file = File::open(paf_path).context("Failed to open PAF file")?;
    let reader = BufReader::new(paf_file);

    let mut total_records = 0;
    let mut records_with_cigar = 0;
    let mut query_lengths = Vec::new();
    let mut target_lengths = Vec::new();
    let mut gap_compressed_identities = Vec::new();
    let mut block_identities = Vec::new();

    // Coverage tracking
    let mut query_coverage: HashMap<String, Vec<(usize, usize)>> = HashMap::new();
    let mut target_coverage: HashMap<String, Vec<(usize, usize)>> = HashMap::new();
    let mut query_sequence_lengths: HashMap<String, usize> = HashMap::new();
    let mut target_sequence_lengths: HashMap<String, usize> = HashMap::new();

    for (line_number, line) in reader.lines().enumerate() {
        let line = line.context("Failed to read PAF line")?;
        let record = PafRecord::from_line(&line).context(format!(
            "Failed to parse PAF record at line {}",
            line_number + 1
        ))?;

        total_records += 1;

        if !record.cigar.is_empty() {
            records_with_cigar += 1;
        }

        let query_mapping_length = record.query_end - record.query_start;
        let target_mapping_length = record.target_end - record.target_start;

        query_lengths.push(query_mapping_length);
        target_lengths.push(target_mapping_length);

        let gap_compressed_identity = if record.block_length > 0 {
            (record.matching_bases as f64 / record.block_length as f64) * 100.0
        } else {
            0.0
        };
        gap_compressed_identities.push(gap_compressed_identity);

        let query_block_identity = if query_mapping_length > 0 {
            (record.matching_bases as f64 / query_mapping_length as f64) * 100.0
        } else {
            0.0
        };
        block_identities.push(query_block_identity);

        // Track coverage
        query_coverage
            .entry(record.query_name.clone())
            .or_default()
            .push((record.query_start, record.query_end));
        target_coverage
            .entry(record.target_name.clone())
            .or_default()
            .push((record.target_start, record.target_end));

        // Track sequence lengths
        query_sequence_lengths.insert(record.query_name, record.query_length);
        target_sequence_lengths.insert(record.target_name, record.target_length);
    }

    if total_records == 0 {
        println!("No PAF records found");
        return Ok(());
    }

    query_lengths.sort_unstable();
    target_lengths.sort_unstable();
    gap_compressed_identities.sort_by(|a, b| a.partial_cmp(b).unwrap());
    block_identities.sort_by(|a, b| a.partial_cmp(b).unwrap());

    println!("PAF Statistics:");
    println!("===============");
    println!("Total records: {total_records}");
    println!(
        "Records with CIGAR: {} ({:.1}%)",
        records_with_cigar,
        (records_with_cigar as f64 / total_records as f64) * 100.0
    );
    println!();

    println!("Query mapping lengths:");
    print_length_stats(&query_lengths);
    println!();

    println!("Target mapping lengths:");
    print_length_stats(&target_lengths);
    println!();

    println!("Gap-compressed identity:");
    print_identity_stats(&gap_compressed_identities);
    println!();

    println!("Block identity (matches/query_length):");
    print_identity_stats(&block_identities);
    println!();

    // Calculate coverage statistics
    let query_coverage_stats = calculate_coverage_stats(&query_coverage, &query_sequence_lengths);
    let target_coverage_stats =
        calculate_coverage_stats(&target_coverage, &target_sequence_lengths);

    println!("Query coverage:");
    print_coverage_stats(&query_coverage_stats);
    println!();

    println!("Target coverage:");
    print_coverage_stats(&target_coverage_stats);

    Ok(())
}

fn print_length_stats(lengths: &[usize]) {
    if lengths.is_empty() {
        println!("  No data");
        return;
    }

    let total: usize = lengths.iter().sum();
    let mean = total as f64 / lengths.len() as f64;
    let median = if lengths.len() % 2 == 0 {
        (lengths[lengths.len() / 2 - 1] + lengths[lengths.len() / 2]) as f64 / 2.0
    } else {
        lengths[lengths.len() / 2] as f64
    };

    println!("  Min: {}", lengths[0]);
    println!("  Max: {}", lengths[lengths.len() - 1]);
    println!("  Mean: {mean:.1}");
    println!("  Median: {median:.1}");
}

fn print_identity_stats(identities: &[f64]) {
    if identities.is_empty() {
        println!("  No data");
        return;
    }

    let total: f64 = identities.iter().sum();
    let mean = total / identities.len() as f64;
    let median = if identities.len() % 2 == 0 {
        (identities[identities.len() / 2 - 1] + identities[identities.len() / 2]) / 2.0
    } else {
        identities[identities.len() / 2]
    };

    println!("  Min: {:.2}%", identities[0]);
    println!("  Max: {:.2}%", identities[identities.len() - 1]);
    println!("  Mean: {mean:.2}%");
    println!("  Median: {median:.2}%");
}

#[derive(Debug)]
struct CoverageStats {
    coverage_percentages: Vec<f64>,
    mean: f64,
    variance: f64,
    skewness: f64,
    kurtosis: f64,
    total_sequences: usize,
}

fn calculate_coverage_stats(
    coverage_map: &HashMap<String, Vec<(usize, usize)>>,
    sequence_lengths: &HashMap<String, usize>,
) -> CoverageStats {
    let mut coverage_percentages = Vec::new();

    for (seq_name, intervals) in coverage_map {
        if let Some(&seq_length) = sequence_lengths.get(seq_name) {
            let covered_bases = calculate_covered_bases(intervals);
            let coverage_pct = (covered_bases as f64 / seq_length as f64) * 100.0;
            coverage_percentages.push(coverage_pct);
        }
    }

    // Add sequences with zero coverage
    for (seq_name, &_seq_length) in sequence_lengths {
        if !coverage_map.contains_key(seq_name) {
            coverage_percentages.push(0.0);
        }
    }

    if coverage_percentages.is_empty() {
        return CoverageStats {
            coverage_percentages: Vec::new(),
            mean: 0.0,
            variance: 0.0,
            skewness: 0.0,
            kurtosis: 0.0,
            total_sequences: 0,
        };
    }

    coverage_percentages.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let mean = coverage_percentages.iter().sum::<f64>() / coverage_percentages.len() as f64;
    let variance = calculate_variance(&coverage_percentages, mean);
    let skewness = calculate_skewness(&coverage_percentages, mean, variance);
    let kurtosis = calculate_kurtosis(&coverage_percentages, mean, variance);

    CoverageStats {
        coverage_percentages,
        mean,
        variance,
        skewness,
        kurtosis,
        total_sequences: sequence_lengths.len(),
    }
}

fn calculate_covered_bases(intervals: &[(usize, usize)]) -> usize {
    if intervals.is_empty() {
        return 0;
    }

    let mut sorted_intervals = intervals.to_vec();
    sorted_intervals.sort_by(|a, b| a.0.cmp(&b.0));

    let mut covered = 0;
    let mut current_end = 0;

    for &(start, end) in &sorted_intervals {
        if start >= current_end {
            covered += end - start;
            current_end = end;
        } else if end > current_end {
            covered += end - current_end;
            current_end = end;
        }
    }

    covered
}

fn calculate_variance(values: &[f64], mean: f64) -> f64 {
    if values.len() <= 1 {
        return 0.0;
    }

    let sum_squared_diff: f64 = values.iter().map(|x| (x - mean).powi(2)).sum();

    sum_squared_diff / (values.len() - 1) as f64
}

fn calculate_skewness(values: &[f64], mean: f64, variance: f64) -> f64 {
    if values.len() < 3 || variance == 0.0 {
        return 0.0;
    }

    let std_dev = variance.sqrt();
    let sum_cubed_diff: f64 = values.iter().map(|x| ((x - mean) / std_dev).powi(3)).sum();

    sum_cubed_diff / values.len() as f64
}

fn calculate_kurtosis(values: &[f64], mean: f64, variance: f64) -> f64 {
    if values.len() < 4 || variance == 0.0 {
        return 0.0;
    }

    let std_dev = variance.sqrt();
    let sum_fourth_diff: f64 = values.iter().map(|x| ((x - mean) / std_dev).powi(4)).sum();

    (sum_fourth_diff / values.len() as f64) - 3.0 // Excess kurtosis
}

fn print_coverage_stats(stats: &CoverageStats) {
    if stats.coverage_percentages.is_empty() {
        println!("  No data");
        return;
    }

    let median = if stats.coverage_percentages.len() % 2 == 0 {
        let mid = stats.coverage_percentages.len() / 2;
        (stats.coverage_percentages[mid - 1] + stats.coverage_percentages[mid]) / 2.0
    } else {
        stats.coverage_percentages[stats.coverage_percentages.len() / 2]
    };

    println!("  Total sequences: {}", stats.total_sequences);
    println!("  Min coverage: {:.2}%", stats.coverage_percentages[0]);
    println!(
        "  Max coverage: {:.2}%",
        stats.coverage_percentages[stats.coverage_percentages.len() - 1]
    );
    println!("  Mean coverage: {:.2}%", stats.mean);
    println!("  Median coverage: {median:.2}%");
    println!("  Std deviation: {:.2}%", stats.variance.sqrt());
    println!("  Skewness: {:.3}", stats.skewness);
    println!("  Kurtosis: {:.3}", stats.kurtosis);
}

fn validate_paf(
    query_fasta: &str,
    target_fasta: &str,
    paf_path: &str,
    error_mode: &str,
    num_threads: usize,
) -> Result<()> {
    // Read all PAF lines into memory with line numbers
    let paf_file = File::open(paf_path).context("Failed to open PAF file")?;
    let reader = BufReader::new(paf_file);

    println!("[pafcheck] Reading PAF file...");
    let lines: Result<Vec<_>> = reader
        .lines()
        .enumerate()
        .map(|(i, line)| line.context("Failed to read PAF line").map(|l| (i + 1, l)))
        .collect();
    let lines = lines?;

    if lines.is_empty() {
        println!("[pafcheck] PAF file is empty. No validation needed.");
        return Ok(());
    }

    println!(
        "[pafcheck] Processing {} PAF records using {} threads",
        lines.len(),
        num_threads
    );

    // Create progress bar
    let progress_bar = Arc::new(ProgressBar::new(lines.len() as u64));
    progress_bar.set_style(
        ProgressStyle::default_bar()
            .template("[pafcheck] {spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} records ({percent}%) {msg}")
            .unwrap()
            .progress_chars("#>-")
    );

    // Determine chunk size
    let chunk_size = lines.len().div_ceil(num_threads);

    // Process chunks in parallel using scoped threads
    let results: Result<Vec<ThreadResult>> = thread::scope(|s| {
        let handles: Vec<_> = lines
            .chunks(chunk_size)
            .map(|chunk| {
                let chunk = chunk.to_vec();
                let query_fasta = query_fasta.to_string();
                let target_fasta = target_fasta.to_string();
                let error_mode = error_mode.to_string();
                let progress = Arc::clone(&progress_bar);

                s.spawn(move || -> Result<ThreadResult> {
                    process_chunk(&query_fasta, &target_fasta, &error_mode, chunk, progress)
                })
            })
            .collect();

        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });

    progress_bar.finish_with_message("Validation complete");

    let results = results?;

    // Aggregate results
    let mut total_error_count = 0;
    let mut error_type_counts: HashMap<ErrorType, usize> = HashMap::new();
    let mut all_error_messages = Vec::new();

    for result in results {
        total_error_count += result.error_count;
        for (error_type, count) in result.error_type_counts {
            *error_type_counts.entry(error_type).or_insert(0) += count;
        }
        all_error_messages.extend(result.error_messages);
    }

    // Sort error messages by line number (extract line number from message)
    all_error_messages.sort_by(|a, b| {
        let line_a = extract_line_number(a).unwrap_or(0);
        let line_b = extract_line_number(b).unwrap_or(0);
        line_a.cmp(&line_b)
    });

    // Print all error messages in order
    for message in all_error_messages {
        println!("{message}");
    }

    // Print summary
    if total_error_count > 0 {
        println!("[pafcheck] PAF validation completed with errors:");
        for (error_type, count) in error_type_counts.iter() {
            println!("[pafcheck]   - {error_type:?}: {count} errors");
        }
        println!("[pafcheck] Total errors: {total_error_count}");
        anyhow::bail!("PAF validation failed with {} errors", total_error_count);
    } else {
        println!("[pafcheck] PAF validation completed successfully. No errors found.");
        Ok(())
    }
}

fn process_chunk(
    query_fasta: &str,
    target_fasta: &str,
    error_mode: &str,
    chunk: Vec<(usize, String)>,
    progress: Arc<ProgressBar>,
) -> Result<ThreadResult> {
    let mut fasta_reader = MultiFastaReader::new(query_fasta, target_fasta)
        .context("Failed to create FASTA readers")?;

    let mut error_count = 0;
    let mut error_type_counts: HashMap<ErrorType, usize> = HashMap::new();
    let mut error_messages = Vec::new();

    for (i, (line_number, line)) in chunk.iter().enumerate() {
        // Update progress
        if i % 100 == 0 {
            progress.inc(100.min(chunk.len() - i) as u64);
            progress.set_message(format!("Processing record {line_number}"));
        }

        let record = PafRecord::from_line(line).context(format!(
            "Failed to parse PAF record at line {line_number}"
        ))?;

        let mut output = Vec::new();
        if let Err(e) = validate_record(&record, &mut fasta_reader, error_mode, &mut output) {
            if let Some(validation_error) = e.downcast_ref::<ValidationError>() {
                for (error_type, error_info) in &validation_error.errors {
                    let count = error_info.count;
                    *error_type_counts.entry(error_type.clone()).or_insert(0) += count;
                    error_count += count;
                    error_messages.push(format!(
                        "[pafcheck] Error at line {}: {:?}: {}",
                        line_number, error_type, error_info.first_message
                    ));
                    if count > 1 {
                        error_messages.push(format!(
                            "[pafcheck] {error_type:?}: Total occurrences: {count}"
                        ));
                    }
                }
            } else {
                error_count += 1;
                error_messages.push(format!("[pafcheck] Error at line {line_number}: {e}"));
            }
        }
    }

    // Update progress for remaining items
    let remaining = chunk.len() % 100;
    if remaining > 0 {
        progress.inc(remaining as u64);
    }

    Ok(ThreadResult {
        error_count,
        error_type_counts,
        error_messages,
    })
}

fn extract_line_number(message: &str) -> Option<usize> {
    // Extract line number from error message format "[pafcheck] Error at line X:"
    if let Some(start) = message.find("line ") {
        let start = start + 5; // Skip "line "
        if let Some(end) = message[start..].find(':') {
            return message[start..start + end].parse().ok();
        }
    }
    None
}
