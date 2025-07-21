# simulate_experiment_countmat.R

## Purpose and Role
This file serves as the **primary entrypoint** for count matrix-based RNA-seq simulation in the polyester package. It provides a high-level interface for generating simulated RNA-seq reads based on a predefined read count matrix, enabling users to control exactly how many reads are generated from each transcript for each sample.

## Pipeline Flow
1. **Input Validation**: Validates input parameters (FASTA vs GTF+seqpath)
2. **Transcript Loading**: Loads transcript sequences from either FASTA or GTF+seqpath
3. **Name Processing**: Handles missing transcript names by adding dummy identifiers
4. **Matrix Validation**: Ensures read count matrix dimensions match transcript count
5. **Parameter Validation**: Validates extra parameters via `.check_extras()`
6. **Output Setup**: Creates output directory structure
7. **Read Generation**: Delegates to `sgseq()` for actual read generation

## Key Variables

### Input Parameters
- `fasta`: Path to FASTA file containing transcript sequences
- `gtf`: Path to GTF file or data frame with transcript structures
- `seqpath`: Path to folder with chromosome FASTA files
- `readmat`: Matrix (transcripts × samples) specifying read counts
- `outdir`: Output directory for simulated reads
- `paired`: Boolean for paired-end vs single-end reads
- `seed`: Optional seed for reproducibility
- `ncores`: Number of cores for parallel processing
- `...`: Additional parameters captured as `extras`

### Internal Variables
- `extras`: List capturing all additional parameters via `...`
- `transcripts`: DNAStringSet object containing transcript sequences
- `sysoutdir`: Sanitized output directory path for system commands

## Complex Structures

### `extras` Structure
The `extras` variable is a list that captures all additional parameters passed via `...`. It contains:
- **Error Models**: `error_model`, `error_rate`, `error_path`
- **Fragment Parameters**: `fraglen`, `fragsd`, `distr`, `custdens`
- **Bias Models**: `bias` (positional), `gcbias` (GC content), `frag_GC_bias`
- **Library Parameters**: `lib_sizes`, `seq_depth`, `strand_specific`
- **Output Options**: `gzip`, `fastq`, `shuffle`, `verbose`
- **Advanced**: `adapter_contamination`, `pcr_rate`, `exon_junction_coverage`

## Key Dependencies
- **`.check_extras()`**: Validates and sets defaults for extra parameters
- **`seq_gtf()`**: Processes GTF files to extract transcript sequences
- **`sgseq()`**: Core function for fragment generation and read writing
- **`readDNAStringSet()`**: Loads sequences from FASTA files

## Usage Notes
This function is designed for users who want precise control over read counts per transcript/sample, as opposed to the more flexible `simulate_experiment()` which uses statistical models to determine read counts.