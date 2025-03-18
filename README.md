# pafcheck

pafcheck is a tool for validating PAF (Pairwise Alignment Format) files against their corresponding FASTA sequences. It ensures that the alignments described in the PAF file match the actual sequences in the FASTA files.

## Installation

To install pafcheck from GitHub, follow these steps:

```bash
git clone https://github.com/ekg/pafcheck.git
cd pafcheck
cargo build --release
```

The compiled binary will be available at `target/release/pafcheck`.

## Usage

```bash
pafcheck -q <query_fasta> -t <target_fasta> -p <paf_file> [-e <error_mode>]
```

- `-q, --query-fasta`: Path to the bgzip-compressed and tabix-indexed query FASTA file
- `-t, --target-fasta`: Path to the bgzip-compressed and tabix-indexed target FASTA file (optional, defaults to query FASTA if not provided)
- `-p, --paf`: Path to the PAF file to validate
- `-e, --error-mode`: Error handling mode: "omit" (default) or "report"

## Error Types Checked

pafcheck validates the following types of errors:

1. **Mismatch**: When the CIGAR string indicates a match, but the actual sequences don't match.
2. **CigarMismatch**: When the CIGAR string indicates a mismatch, but the actual sequences match.
3. **LengthMismatch**: When the length implied by the CIGAR string doesn't match the actual sequence length.

## Generating Input Files

To generate input files that can be validated with pafcheck, you can use tools like wfmash and minimap2. Here are example commands:

### Using wfmash:

```bash
wfmash tests/data/a.fa tests/data/b.fa > w.paf
```

### Using minimap2:

```bash
minimap2 -cx asm20 --eqx tests/data/a.fa tests/data/b.fa > m.paf
```

Note: The `--eqx` option is important for minimap2 as it generates CIGAR strings with '=' for matches and 'X' for mismatches, which is required for proper validation with pafcheck.

After generating these PAF files, you can validate them using pafcheck:

```bash
pafcheck -q tests/data/a.fa -t tests/data/b.fa -p w.paf
pafcheck -q tests/data/a.fa -t tests/data/b.fa -p m.paf
```

Make sure to bgzip-compress and index your FASTA files before running pafcheck:

```bash
bgzip tests/data/a.fa
bgzip tests/data/b.fa
samtools faidx tests/data/a.fa.gz
samtools faidx tests/data/b.fa.gz
```

## Development

To compile pafcheck with guix or get a development shell see [guix.scm](./guix.scm).

## Contributing

Contributions to pafcheck are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT license. Have fun.
