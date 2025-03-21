# Important Note

This is a customized polyester repository used for the [ASimulatoR package](https://github.com/daisybio/ASimulatoR). Please do not install directly. This is only to be used by [ASimulatoR](https://github.com/daisybio/ASimulatoR).

# Introduction

Polyester is an R package designed to simulate RNA sequencing experiments with differential transcript expression. 

Given a set of annotated transcripts, Polyester will simulate the steps of an RNA-seq experiment (fragmentation, reverse-complementing, and sequencing) and produce files containing simulated RNA-seq reads. Simulated reads can be analyzed using your choice of downstream analysis tools. 

Polyester has a built-in wrapper function to simulate a case/control experiment with differential transcript expression and biological replicates. Users are able to set the levels of differential expression at transcripts of their choosing. This means they know which transcripts are differentially expressed in the simulated dataset, so accuracy of statistical methods for differential expression detection can be analyzed. 

Polyester offers several unique features:
* Built-in functionality to simulate differential expression at the transcript level
* Ability to explicitly set differential expression signal strength
* Simulation of small datasets, since large RNA-seq datasets can require lots of time and computing resources to analyze
* Generation of raw RNA-seq reads, as opposed to alignments or transcript-level abundance estimates
* Transparency/open-source code

# Installation

Start R and run:


```r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("polyester")
```

# Required Input
You'll need to provide transcript annotation from which reads should be simulated. There are several public data repositories where you can download this annotation. You can simulate reads from any organism for which annotation is available.

Annotation must be provided in one of two formats:

1. [FASTA](https://earray.chem.agilent.com/earray/helppages/fasta_format_files.htm): text file containing names and sequences of transcripts from which reads should be simulated. Known transcripts from human chromosome 22 (hg19 build) are available in `extdata/chr22.fa`. 
* [GTF format](http://www.ensembl.org/info/website/upload/gff.html) + FASTA sequence files. The GTF file should denote the transcript structures, and you'll need a FASTA file of the full DNA sequence for each chromosome in the GTF file. All the chromosome-specific FASTA files should be in the same directory. 

GTF files and DNA sequences for many widely-studied organisms can be downloaded [here](http://ccb.jhu.edu/software/tophat/igenomes.shtml) (sequences are in the `<organism>/<source>/<build>/Sequence/Chromosomes` subdirectory, e.g., `Homo_sapiens/UCSC/hg19/Sequence/Chromosomes`).

# Simulating reads

Simulating an RNA-seq experiment with Polyester generally requires just one function call. In the simplest case, you will use the `simulate_experiment()` function.

### `simulate_experiment` example

A FASTA file called `chr22.fa` is provided with `polyester`. This file contains sequences for 918 transcripts on chromosome 22, as annotated in hg19. For this very small example, we will only simulate from the first 20 of these transcripts.

We will set the first 2 transcripts to be overexpressed in group A and the next 2 transcripts to be overexpressed in group B, each at a fold change of 4. The way to do this in Polyester is to provide a "fold change matrix": for each transcript and each group, specify a fold change. Polyester will generate baseline read numbers (assuming no differential expression), and will then multiply those mean numbers by the fold change you specify for the replicates in that group. The fold change matrix for this simple 2-group experiment looks like this:


```r
fold_changes = matrix(c(4,4,rep(1,18),1,1,4,4,rep(1,16)), nrow=20)
head(fold_changes)
```

```
##      [,1] [,2]
## [1,]    4    1
## [2,]    4    1
## [3,]    1    4
## [4,]    1    4
## [5,]    1    1
## [6,]    1    1
```
The matrix has two columns, since there will be two groups (cases and controls) in this experiment.

The rest of the experiment can be simulated with code like the chunk below. 


```r
library(polyester)
library(Biostrings)

# FASTA annotation
fasta_file = system.file('extdata', 'chr22.fa', package='polyester')
fasta = readDNAStringSet(fasta_file)

# subset the FASTA file to first 20 transcripts
small_fasta = fasta[1:20]
writeXStringSet(small_fasta, 'chr22_small.fa')

# ~20x coverage ----> reads per transcript = transcriptlength/readlength * 20
# here all transcripts will have ~equal FPKM
readspertx = round(20 * width(small_fasta) / 100)

# simulation call:
simulate_experiment('chr22_small.fa', reads_per_transcript=readspertx, 
    num_reps=c(10,10), fold_changes=fold_changes, outdir='simulated_reads') 
```

The documentation for the `simulate_experiment` function explains these arguments in detail. Briefly, we will need baseline abundances for each transcript, the number of replicates per group (here we have 2 groups with 10 replicates each, but you can have as many groups as you like), the fold change matrix (fill this with 1's if you don't want any differential expression), and a directory where the fasta files containing sequencing reads and simulation information should be written.

The statistical model behind this function is a negative binomial distribution, which appropriately captures both biological and technical variability in expression measurements (Anders and Huber 2010; Robinson, McCarthy, and Smyth 2010). The important parameters for the negative binomial model are:

* `reads_per_transcript`: The baseline mean number of reads for each transcript.  
    - Fold changes will multiply these baseline mean numbers for the specified group. 
    - Long transcripts usually produce more reads in RNA-seq experiments than short ones, so you may want to specify `reads_per_transcript` as a function of transcript length
    - Default is 300 (regardless of transcript length).
* `size`: controls the per-transcript mean/variance relationship. In the negative binomial distribution, the mean/variance relationship is: ```mean = mean + (mean^2) / size```. You can specify the size for each transcript. By default, size is defined as 1/3 of the transcript's mean, which (in our experience) creates an idealized, low-variance situation.  Decrease the value of `size` to introduce more variance into your simulations.


### `simulate_experiment_countmat` example
If you're an experienced user requiring more flexibility, you can use the `simulate_experiment_countmat` function to directly specify the number of reads you'd like to simulate for each transcript and each replicate in the data set. This function takes a count matrix as an argument. Each row of this matrix represents a transcript, and each column represents a sample in the experiment. Entry `i,j` of the matrix specifies how many reads should be sampled from transcript `i` for sample `j`.

For example, let's say we want to simulate timecourse data. To do this, we will explicitly specify a number of reads to generate for each transcript (row), at each timepoint (column). We will again only simulate from 20 transcripts.


```r
# set up transcript-by-timepoint matrix:
num_timepoints = 12
countmat = matrix(readspertx, nrow=length(small_fasta), ncol=num_timepoints)

# add spikes in expression at certain timepoints to certain transcripts:
up_early = c(1,2) 
up_late = c(3,4)
countmat[up_early, 2] = 3*countmat[up_early, 2]
countmat[up_early, 3] = round(1.5*countmat[up_early, 3])
countmat[up_late, 10] = 6*countmat[up_late, 10]
countmat[up_late, 11] = round(1.2*countmat[up_late, 11])

# simulate reads:
simulate_experiment_countmat('chr22_small.fa', readmat=countmat, 
    outdir='timecourse_reads') 
```

In this scenario, we simulated 12 total timepoints. We also added differential expression: transcripts 1 and 2 are overexpressed (compared to baseline) at timepoints 2 and 3, with a fold change of 3 at timepoint 2 and a fold change of 1.5 at timepoint 3. Similarly, transcripts 3 and 4 are overexpressed at timepoints 10 and 11, with a fold change of 6 at timepoint 10 and a fold change of 1.2 at timepoint 11. 

### Parallel processing
Using the `parallel` package, you can now specify how many cores you wish to use to process multiple samples simultaneously. For any of the `simulate_experiment*`, you can specify the number of cores to use. Note that the internal function will split each sample across cores, so using multiple cores will only give a benefit if you are simulating more than one sample.

### Changing other parameters and adding bias 
The following parameters can be provided to `simulate_experiment` and `simulate_experiment_countmat`:

* `readlen`: Read length (default 100)
* `paired`: Whether the reads should be paired-end (default TRUE)
* `distr`: Distribution from which to draw the fragment lengths. Default is `'normal'`, mean=250 and sd=25. Other options are `'empirical'` and `'custom'`. `'empirical'` means fragment lengths are drawn from a length distribution Frazee et al estimated from real data, and `'custom'` requires you to provide a `logspline` density object from which you'd like to draw the fragment lengths.
* `fraglen`: The mean fragment length, if using a normal distribution (default 250)
* `fragsd`: Standard devation of fragment lengths, if using a normal distribution (default 25)
* `error_model`: How should sequencing errors be simulated? The default is that sequencing errors are uniformly distributed across samples, reads, and nucleotides, at an error rate of 0.5\%. Other options are `'illumina4'`, `'illumina5'`, and `'custom'`, where the Illumina error models were estimated from a real data set and ship with GemSIM (McElroy, Luciani, and Thomas 2012), and `'custom'` allows you to use GemSIM to estimate an error model from your data set. See `?add_platform_error` and/or the [GemSIM paper](http://www.biomedcentral.com/1471-2164/13/74) for details. Code Frazee et al used to modify GemSIM's Illumina error models for Polyester are available at the [GitHub repository](https://github.com/biomedbigdata/polyester/blob/master/ErrorModels.md). 
* `error_rate`: In the uniform error model, probability that the sequencer records the wrong nucleotide at any given base (default 0.005). 
* `bias`: Positional bias model to use when fragmenting transcripts. By default, all fragments from a transcript are equally likely (`'none'`). Other choices are `'rnaf'` and `'cdnaf'`, which mimic positional bias arising from different fragmentation protocols. See `?generate_fragments` and the manuscript (Frazee et al, 2014) for details.
* `gc_bias`: sample-specific GC bias models to be used to change expression values after read numbers are assigned. Frazee et al modeled transcript expression as a function of GC content for 7 biological replicates in a real data set, and shift expression values accordingly. See `?add_gc_bias` for details. Ignored in `simulate_experiment_countmat`.
* `frag_GC_bias`: A sample-specific GC content bias on the
  fragment level instead of the transcript level. See `?simulate_experiment`
* `strand_specific`: Whether the experiment should be strand-specific
  or not, default is `FALSE` so an unstranded experiment.
* `meanmodel`: Set to TRUE to estimate `reads_per_transcript` as a data-driven function of transcript length. Ignored in `simulate_experiment_countmat`. 
* `lib_sizes`: multiplicative library size factor for each replicate in the experiment. Ignored in `simulate_experiment_countmat`. 
* `shuffle`: Set to TRUE to shuffle each chunk of simulated reads before writing them to file. This is important for most downstream quantification tools, as they expect the reads to not be in order.
* `fastq`: Set to TRUE to write reads into fastq files instead of fasta files.
* `verbose`: Set to TRUE to print progress messages during the sequencing process
* `seq_depth`: Number of reads to be sequenced per sample. Can be a vector of length one or `sum(num_reps)`. Readcounts will be multiplied by a factor to be equal to `seq_depth`. If used with `lib_sizes`, `seq_depth` will be applied first.
* `pcr_rate`: Fraction of fragments that will be duplicated. Reads from these fragments will have PCR_DUP in the name.
* `pcr_lambda`: If `!is.null(pcr_rate)` lambda for the poisson distribution to draw the number of duplicates.
* `adapter_contamination`: If the fragment is smaller than the readlength, should we sequence into the `adapter_sequence`?
* `adapter_sequence`:If `adapter_contamination`: adapter sequence

For most of these parameters, you can see additional, precise documentation using `?simulate_experiment`. Also, [this review paper](http://genomebiology.com/2010/11/12/220) (Oshlack, Robinson, and Young, _Genome Biology_ 2010, open access) provides a good overview of the RNA sequencing process, and might be particularly useful for understanding where some of these simulation parameters come into play. If you'd like to explore or change specific steps in the sequencing process (fragmentation, reverse-complementing, error-adding), the internal functions called within `simulate_experiment` are available and individually documented in Polyester.

### Using real data to guide simulation

We also provide a function `simulate_experiment_empirical` that takes a real transcript expression matrix (in FPKM or RPKM units) and corresponding annotation to simulate an experiment with abundances and differential expression fold changes similar to those given in the expression matrix. This function is compatible with the [Ballgown package](http://www.bioconductor.org/packages/release/bioc/html/ballgown.html), or you can simply provide a transcript-by-replicate expression matrix generated with your favorite abundance estimation software. 

To create a count matrix that resembles a real dataset, use the `create_read_numbers` function. To run this example, you will need to install the Ballgown package from Bioconductor if you do not already have it:


```r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("ballgown")
```

```r
library(ballgown)
data(bg)
bg = subset(bg, "chr=='22'")

# load gtf file for annotation:
gtfpath = system.file('extdata', 'bg.gtf.gz', package='polyester')
gtf = subset(gffRead(gtfpath), seqname=='22')

# load chromosome sequence corresponding to gtf file (just for this example)
system('wget https://www.dropbox.com/s/04i6msi9vu2snif/chr22seq.rda')
load('chr22seq.rda')
names(chr22seq) = '22'

# simulate reads based on this experiment's FPKMs
simulate_experiment_empirical(bg, grouplabels=pData(bg)$group, gtf=gtf,
    seqpath=chr22seq, mean_rps=5000, outdir='empirical_reads', seed=1247)
```


```
## [1] "fold change direction: 1/0"
```

# Output
A call to `simulate_experiment`, `simulate_experiment_countmat`, or `simulate_experiment_empirical` will write FASTA files to the directory specified by the `outdir` argument. Reads in the FASTA file will be labeled with the transcript from which they were simulated, based on the identifiers provided in the initial input, as well as the 1-based start and end positions of the simulated reads. The positions for mate 1 and mate 2 are provided in the **same orientation**, as if the transcript were being scanned for read coverage, even though mate 2 is aligned in the **opposite** orientation. Make sure to parse these positions correctly based on how your read alignments are reported.

If `paired` is true, you'll get two FASTA files per biological replicate (left mates are designated by the suffix `_1.fasta`; right mates by `_2.fasta`). If single-end reads are generated (`paired=FALSE`) you'll get one FASTA file per replicate. 

Files will be named `sample_01` through `sample_N` where `N` is the total number of replicates (i.e., `sum(num_reps)` in `simulate_experiment`). Samples are grouped in order: e.g., if you provide `num_reps = c(3,6,2)`, samples 01 through 03 belong to group A, samples 04 through 09 belong to group B, and samples 10 and 11 belong to group C. 

In `simulate_experiment` and `simulate_experiment_emprical`, simulation information is written to disk by default. There will be 2 files: `sim_tx_info.txt`, giving transcript IDs and differential expression fold changes and statuses, and `sim_rep_info.txt`, giving sample IDs, library sizes, and group designations. These files are essential in downstream analyses. You will need to keep track of this information separately if you use `simulate_experiment_countmat.`

# Bug reports
Bug reports are very welcome as issues on our [GitHub repository](https://github.com/alyssafrazee/polyester). 

# Session Information


```r
sessionInfo()
```

```
## R version 3.1.1 (2014-07-10)
## Platform: x86_64-apple-darwin13.1.0 (64-bit)
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] ballgown_1.0.1      Biostrings_2.34.0   XVector_0.6.0      
##  [4] IRanges_2.0.0       S4Vectors_0.4.0     BiocGenerics_0.12.1
##  [7] polyester_1.2.0     knitr_1.8.2         colorout_1.0-3     
## [10] devtools_1.6.1     
## 
## loaded via a namespace (and not attached):
##  [1] annotate_1.44.0         AnnotationDbi_1.28.1   
##  [3] base64enc_0.1-2         BatchJobs_1.5          
##  [5] BBmisc_1.8              Biobase_2.26.0         
##  [7] BiocParallel_1.0.0      bitops_1.0-6           
##  [9] brew_1.0-6              checkmate_1.5.0        
## [11] codetools_0.2-9         DBI_0.3.1              
## [13] digest_0.6.4            evaluate_0.5.5         
## [15] fail_1.2                foreach_1.4.2          
## [17] formatR_1.0             genefilter_1.48.1      
## [19] GenomeInfoDb_1.2.3      GenomicAlignments_1.2.1
## [21] GenomicRanges_1.18.3    grid_3.1.1             
## [23] iterators_1.0.7         lattice_0.20-29        
## [25] limma_3.22.1            logspline_2.1.5        
## [27] markdown_0.7.4          Matrix_1.1-4           
## [29] mgcv_1.8-4              nlme_3.1-118           
## [31] RColorBrewer_1.1-2      RCurl_1.95-4.5         
## [33] Rsamtools_1.18.2        RSQLite_1.0.0          
## [35] rtracklayer_1.26.2      sendmailR_1.2-1        
## [37] splines_3.1.1           stringr_0.6.2          
## [39] survival_2.37-7         sva_3.12.0             
## [41] tools_3.1.1             XML_3.98-1.1           
## [43] xtable_1.7-4            zlibbioc_1.12.0
```

