# sgseq.R

## Purpose and Role
This file contains the **core sequencing engine** of the polyester package. The `sgseq()` function performs the actual read generation process, handling fragment generation, read extraction, error addition, and file writing. It operates as the final stage in the simulation pipeline, converting validated parameters and transcript sequences into actual FASTA/FASTQ files.

## Pipeline Flow
1. **Initialization**: Sets up coverage tracking if requested
2. **Parallel Processing**: Uses mclapply for parallel sample processing
3. **Fragment Generation**: Generates fragments from transcripts based on count matrix
4. **Read Extraction**: Extracts reads from fragments (paired or single-end)
5. **Error Addition**: Applies sequencing errors based on specified model
6. **File Writing**: Writes reads to FASTA/FASTQ files in chunks
7. **Coverage Reporting**: Optionally writes coverage information

## Key Variables

### Input Parameters
- `readmat`: Count matrix (transcripts × samples) specifying read counts
- `transcripts`: DNAStringSet containing transcript sequences
- `paired`: Boolean for paired-end vs single-end reads
- `outdir`: Output directory for generated files
- `extras`: Validated parameter list from `.check_extras()`
- `reportCoverage`: Boolean for coverage reporting (default: FALSE)
- `ncores`: Number of cores for parallel processing

### Internal Variables
- `tObj`: Replicated transcript sequences based on read counts
- `iterations`: Number of 1M-read chunks needed
- `tFrags`: Generated fragments from transcripts
- `reads`: Extracted reads from fragments
- `errReads`: Reads with sequencing errors applied
- `region_counts`: Exon junction coverage counts (if enabled)

## Complex Structures

### Chunk Processing System
- **Memory Management**: Processes reads in 1 million read chunks to balance speed and memory
- **Offset Tracking**: Maintains read numbering across chunks
- **Progress Reporting**: Provides verbose progress messages per iteration

### Error Model Integration
- **Uniform Errors**: Uses `add_error()` for uniform error distribution
- **Platform Errors**: Uses `add_platform_error()` for empirical/custom models
- **Model Selection**: Chooses error addition method based on `extras$error_model`

### Coverage Tracking
- **Template Coverage**: Tracks coverage across unique transcript templates
- **Coordinate Mapping**: Maps reads to genomic coordinates for coverage calculation
- **File Output**: Writes coverage matrices to `sample_coverages.rda`

### Exon Junction Coverage
- **Data Structure**: Uses data.table for efficient junction counting
- **ID Tracking**: Tracks covered exon junctions with unique identifiers
- **Output Format**: Writes junction coverage to TSV files

## Key Dependencies
- **`generate_fragments()`**: Creates fragments from transcript sequences
- **`get_reads()`**: Extracts reads from fragments
- **`add_error()`**: Adds uniform sequencing errors
- **`add_platform_error()`**: Adds platform-specific errors
- **`write_reads()`**: Writes reads to output files
- **`reverse_complement()`**: Handles strand-specific simulation
- **parallel::mclapply**: Enables parallel processing across samples

## Usage Notes
This function is **not intended for direct user use** - it's designed to be called internally by `simulate_experiment()`, `simulate_experiment_countmat()`, and `simulate_experiment_empirical()`. It provides the computational core that transforms simulation parameters into actual sequencing data.