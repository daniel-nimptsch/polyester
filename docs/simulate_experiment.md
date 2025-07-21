# simulate_experiment.R

## Purpose and Role
This file provides the **statistical simulation engine** for RNA-seq experiments in the polyester package. Unlike `simulate_experiment_countmat.R` which uses predefined counts, this function uses statistical models (negative binomial) to generate read counts based on fold changes and baseline expression levels, making it ideal for differential expression studies.

## Pipeline Flow
1. **Parameter Validation**: Validates all input parameters and extra arguments
2. **Transcript Loading**: Loads transcript sequences from FASTA or GTF+seqpath
3. **Fold Change Validation**: Ensures fold change matrix dimensions are correct
4. **Baseline Calculation**: Computes baseline read counts with fold changes applied
5. **Statistical Modeling**: Uses negative binomial distribution to generate read counts
6. **Bias Application**: Applies GC bias and other optional modifications
7. **Output Setup**: Creates output directory and writes simulation info
8. **Read Generation**: Delegates to `sgseq()` for actual read generation

## Key Variables

### Input Parameters
- `fasta`: Path to FASTA file with transcript sequences
- `gtf`: Path to GTF file with transcript annotations
- `seqpath`: Path to folder with chromosome FASTA files
- `outdir`: Output directory for simulated reads
- `num_reps`: Vector specifying replicates per group (default: c(10,10))
- `reads_per_transcript`: Baseline reads per transcript (default: 300)
- `size`: Negative binomial size parameter (default: basemeans/3)
- `fold_changes`: Matrix of multiplicative fold changes between groups
- `paired`: Boolean for paired-end reads (default: TRUE)
- `reportCoverage`: Boolean for coverage reporting (default: FALSE)
- `ncores`: Number of cores for parallel processing
- `...`: Additional parameters captured as `extras`

### Internal Variables
- `basemeans`: Matrix of baseline read counts after fold changes
- `group_ids`: Vector mapping samples to groups
- `numreadsList`: List of negative binomial samples per group
- `readmat`: Final count matrix after all adjustments
- `counts_matrix`: Unbiased counts matrix for output

## Complex Structures

### Statistical Modeling
- **Negative Binomial**: Uses NB distribution with mean = basemeans, size = basemeans/3
- **Fold Change Application**: Multiplies baseline counts by fold changes
- **Variance Structure**: Variance = mean + mean²/size

### Parameter Calculation
- **Baseline Means**: Can use length-based model via `meanmodel=TRUE`
- **Size Parameter**: Automatically calculated as basemeans/3 if not provided
- **Library Size Adjustment**: Applies library size factors after count generation
- **Sequencing Depth**: Adjusts counts to match target sequencing depth

### Fold Change Matrix
- **Dimensions**: Must match number of transcripts × number of groups
- **Values**: Multiplicative factors applied to baseline counts
- **DE Status**: Transcripts with fold changes ≠ 1 marked as differentially expressed

### Output Information
- **sim_tx_info.txt**: Contains transcript IDs, fold changes, and DE status
- **sim_rep_info.txt**: Contains sample IDs, group assignments, and library sizes
- **sim_counts_matrix.rda**: R data file with unbiased count matrix

## Key Dependencies
- **`.check_extras()`**: Validates extra parameters
- **`.check_fold_changes()`**: Validates fold change matrix dimensions
- **`.check_error_model()`**: Validates error model parameters
- **`.write_info()`**: Writes simulation information files
- **`seq_gtf()`**: Processes GTF files for transcript extraction
- **`add_gc_bias()`**: Applies GC bias to read counts
- **`NB()`**: Negative binomial sampling function
- **`sgseq()`**: Core read generation function

## Usage Notes
This function is ideal for **differential expression studies** where you want to simulate realistic RNA-seq data with known fold changes between groups. The statistical modeling approach provides more realistic variation compared to fixed count matrices.