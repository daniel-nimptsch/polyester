# seq_gtf.R

## Purpose and Role
This file provides the **transcript extraction functionality** for GTF-based RNA-seq simulation. It processes GTF files (containing transcript structure information) along with genomic DNA sequences to generate transcript sequences, serving as a bridge between genomic annotation and simulated RNA-seq data.

## Pipeline Flow
1. **Input Validation**: Validates GTF input (file path vs data frame)
2. **GTF Processing**: Reads and validates GTF file structure
3. **Chromosome Validation**: Ensures all chromosomes have corresponding sequence files
4. **Sequence Loading**: Loads DNA sequences for each chromosome
5. **Transcript Assembly**: Extracts and assembles transcript sequences from exons
6. **Strand Handling**: Reverse complements sequences for negative strand transcripts
7. **Output**: Returns DNAStringSet of transcript sequences

## Key Variables

### Input Parameters
- `gtf`: Path to GTF file or data frame with transcript annotations
- `seqs`: Path to folder with chromosome FASTA files or named DNAStringSet
- `feature`: Either 'transcript' (default) or 'exon' for output type
- `exononly`: Boolean to use only exon features (default: TRUE)
- `idfield`: Field name in GTF attributes for transcript IDs (default: "transcript_id")
- `attrsep`: Separator for GTF attributes (default: "; ")

### Internal Variables
- `gtf_dat`: Processed GTF data frame with standardized column names
- `chrs`: Unique chromosome names from GTF
- `seqlist`: List of DNAStringSets for each chromosome
- `full_list`: Combined DNAStringSet of all sequences

## Complex Structures

### GTF Data Frame Structure
The processed GTF data frame contains these standardized columns:
- `seqname`: Chromosome/contig name
- `source`: Annotation source
- `feature`: Feature type (exon, transcript, etc.)
- `start`: Genomic start position
- `end`: Genomic end position
- `score`: Quality score
- `strand`: Strand orientation (+ or -)
- `frame`: Reading frame
- `attributes`: Additional feature attributes

### Sequence Processing
- **File Discovery**: Automatically finds chromosome files with extensions: 'fa', 'fasta', 'fna', 'fas'
- **Validation**: Ensures all chromosomes in GTF have corresponding sequence files
- **Assembly**: Combines exon sequences to form complete transcript sequences
- **Strand Correction**: Reverse complements sequences on negative strand

## Key Dependencies
- **`readDNAStringSet()`**: Loads DNA sequences from FASTA files
- **`subseq()`**: Extracts subsequences from genomic coordinates
- **`reverseComplement()`**: Handles strand orientation
- **`getAttributeField()`**: Parses transcript IDs from GTF attributes

## Usage Notes
This function is essential for GTF-based simulation workflows, providing the transcript sequences needed for read generation. It handles both file-based and in-memory sequence inputs, making it flexible for different use cases.