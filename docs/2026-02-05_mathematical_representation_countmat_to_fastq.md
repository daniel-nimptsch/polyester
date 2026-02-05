# Mathematical Representation of Polyester Pipeline: Count Matrix to FASTQ Generation

**Date:** 2026-02-05  
**Pipeline Section:** Transcript count matrix → FASTQ generation  
**Focus:** Bioinformatic science details and mathematical models

## Overview of Pipeline Flow

The Polyester RNA-seq simulation pipeline transforms a predefined read count matrix into simulated FASTQ files through the following sequential steps:

```
Input: Read count matrix R[N×M]
       ↓
1. Fragment Generation: Generate fragments from transcripts
       ↓
2. Read Extraction: Extract reads from fragments (single/paired-end)
       ↓
3. Error Introduction: Add sequencing errors
       ↓
4. Output: Write FASTQ files
```

Where:
- $N$ = number of transcripts
- $M$ = number of samples/replicates
- $R_{ij}$ = number of reads to simulate from transcript $i$ in sample $j$

## 1. Mathematical Models

### 1.1 Count Matrix Model

The input is a user-defined matrix $R \in \mathbb{N}^{N \times M}$ where:

$$
R_{ij} \in \mathbb{N} \quad \text{for } i = 1,\ldots,N; \quad j = 1,\ldots,M
$$

This represents **exact control** over read counts per transcript per sample, unlike the statistical model in `simulate_experiment()` which uses negative binomial distributions.

### 1.2 Fragment Generation Model

For each transcript $i$ in sample $j$, generate $R_{ij}$ fragments. The fragment generation involves:

#### 1.2.1 Fragment Length Distribution

Three distributions available:

1. **Normal distribution** (default):
   $$
   L_{ij} \sim \mathcal{N}(\mu_{fl}, \sigma_{fl}^2)
   $$
   where $\mu_{fl} = 250$ bp, $\sigma_{fl} = 25$ bp by default.

2. **Empirical distribution** (from GEUVADIS data):
   $$
   L_{ij} \sim \text{Empirical}(\text{logspline density})
   $$

3. **Custom logspline distribution**:
   $$
   L_{ij} \sim \text{Custom}(f_{\text{cust}})
   $$
   where $f_{\text{cust}}$ is a user-provided logspline density.

#### 1.2.2 Fragment Starting Position Bias

Three positional bias models:

1. **Uniform (no bias)**:
   $$
   P(\text{start} = s) = \frac{1}{L_i - L_{ij} + 1}
   $$
   for $s = 1, \ldots, L_i - L_{ij} + 1$

2. **RNA fragmentation bias (`rnaf`)**:
   $$
   P(\text{start} = s) \propto f_{\text{rnaf}}\left(\frac{s}{L_i}\right)
   $$
   where $f_{\text{rnaf}}$ is the probability model from Li & Jiang (2012)

3. **cDNA fragmentation bias (`cdnaf`)**:
   $$
   P(\text{start} = s) \propto f_{\text{cdnaf}}\left(\frac{s}{L_i}\right)
   $$
   where $f_{\text{cdnaf}}$ is the probability model from Li & Jiang (2012)

#### 1.2.3 GC Content Bias

For fragment $f$ with GC content $g_f \in [0,1]$:

$$
P(\text{retain fragment } f) = \phi_j(g_f)
$$

where $\phi_j: [0,1] \rightarrow [0,1]$ is a sample-specific GC bias function. Built-in models $\phi_1, \ldots, \phi_7$ are loess fits from GEUVADIS data.

#### 1.2.4 PCR Amplification Bias

Optional PCR duplication with parameters:
- $p_{\text{PCR}}$ = rate of fragments undergoing PCR
- $\lambda_{\text{PCR}}$ = Poisson mean for duplication count

For each fragment $f$:
$$
\text{PCR}(f) = 
\begin{cases}
\text{Not duplicated} & \text{with probability } 1 - p_{\text{PCR}} \\
\text{Duplicated } D_f \text{ times} & \text{with probability } p_{\text{PCR}}
\end{cases}
$$
where $D_f \sim \text{Poisson}(\lambda_{\text{PCR}}) + 1$

### 1.3 Strand Orientation Model

For unstranded protocols (default), each fragment is reverse-complemented with probability 0.5:

$$
\text{Orientation}(f) = 
\begin{cases}
\text{Forward} & \text{with probability } 0.5 \\
\text{Reverse complement} & \text{with probability } 0.5
\end{cases}
$$

### 1.4 Read Extraction Model

#### 1.4.1 Single-end reads:
$$
\text{Read}(f) = \text{subseq}(f, 1, \min(R, |f|))
$$
where $R$ = read length (default 100 bp)

#### 1.4.2 Paired-end reads:
- Mate 1: $\text{subseq}(f, 1, \min(R, |f|))$
- Mate 2: $\text{reverse\_complement}(\text{subseq}(f, |f|-R+1, |f|))$

For fragments shorter than read length $R$:
- Single-end: Use entire fragment
- Paired-end: Both mates use entire fragment, mate 2 is reverse complement

### 1.5 Sequencing Error Models

#### 1.5.1 Uniform error model:
For each nucleotide position $k$ in read $r$:
$$
P(\text{error at } k) = \epsilon
$$
$$
P(\text{base } b' | \text{correct base } b) = \frac{1}{4} \quad \text{for } b' \neq b
$$
Default $\epsilon = 0.005$

#### 1.5.2 Empirical error models (from GemSim):
Position- and base-dependent error probabilities:
$$
P_{b,k}(\text{error}) = \epsilon_{b,k}
$$
$$
P_{b,k}(b' | b) = p_{b,k,b'} \quad \text{for } b' \in \{A,C,G,T,N\}
$$
where $\epsilon_{b,k}$ and $p_{b,k,b'}$ are estimated from real sequencing data.

#### 1.5.3 Custom error models:
User-provided error probability matrices.

## 2. Mathematical Formulation of Core Functions

### 2.1 `generate_fragments()` Function

For transcript $t$ of length $L_t$:

1. Draw fragment length: $L_f \sim \text{Distribution}(\theta)$
2. If $L_f < L_t$:
   - Draw start position: $s \sim \text{BiasModel}(L_t, L_f)$
   - Fragment: $f = t[s:(s+L_f-1)]$
3. Apply GC bias: retain with probability $\phi(g_f)$
4. Apply PCR bias: duplicate if selected
5. Apply strand orientation: reverse complement with probability 0.5

### 2.2 `get_reads()` Function

For fragment $f$ of length $L_f$:

If single-end:
$$
\text{read} = 
\begin{cases}
f[1:R] & \text{if } L_f \geq R \\
f & \text{if } L_f < R
\end{cases}
$$

If paired-end:
$$
\text{mate1} = 
\begin{cases}
f[1:R] & \text{if } L_f \geq R \\
f & \text{if } L_f < R
\end{cases}
$$
$$
\text{mate2} = 
\begin{cases}
\text{rc}(f[L_f-R+1:L_f]) & \text{if } L_f \geq R \\
\text{rc}(f) & \text{if } L_f < R
\end{cases}
$$

### 2.3 `add_error()` Function

For read $r$ of length $L_r$:

For each position $k = 1,\ldots,L_r$:
1. With probability $\epsilon_{\text{adj}} = \epsilon \times \frac{5}{4}$: replace base
2. New base drawn uniformly from $\{A,C,G,T,N\}$

### 2.4 `sgseq()` Function (Sequencing Core)

Algorithm for sample $j$:

1. **Input**: $R_{:,j}$ (column $j$ of count matrix)
2. **Replicate transcripts**: $T_j = \text{rep}(\text{transcripts}, R_{:,j})$
3. **Process in chunks** of 1M fragments:
   - Generate fragments: $F = \text{generate\_fragments}(T_j^{\text{chunk}})$
   - Extract reads: $R_{\text{reads}} = \text{get\_reads}(F)$
   - Add errors: $R_{\text{err}} = \text{add\_error}(R_{\text{reads}})$
   - Write to FASTQ

## 3. Key Parameters and Their Effects

### 3.1 Fragment Parameters
- $\mu_{fl}$, $\sigma_{fl}$: Affect insert size distribution
- Bias model: Affects coverage uniformity across transcripts

### 3.2 Error Parameters
- $\epsilon$: Affects base quality and mapping accuracy
- Error model type: Affects error patterns (uniform vs. empirical)

### 3.3 Bias Parameters
- GC bias functions $\phi_j$: Affect expression estimates for GC-rich/poor transcripts
- PCR parameters: Affect duplicate rate and coverage uniformity

## 4. Biological Interpretation

### 4.1 Transcript Representation
Each transcript $i$ is represented by its DNA sequence $S_i$. The count matrix $R_{ij}$ represents the **molecular abundance** of transcript $i$ in sample $j$.

### 4.2 Fragmentation Process
Models the physical/chemical fragmentation of RNA/cDNA during library preparation.

### 4.3 Sequencing Process
Models the Illumina sequencing-by-synthesis process, including:
- Read length constraints
- Paired-end geometry
- Sequencing errors

### 4.4 Technical Biases
Models known technical artifacts:
- GC bias: Differential amplification based on GC content
- Positional bias: Non-uniform fragmentation
- PCR bias: Duplication artifacts

## 5. Mathematical Summary

The complete transformation from count matrix to FASTQ can be expressed as:

$$
\text{FASTQ}_j = \mathcal{F}_j(R_{:,j}, S, \Theta)
$$

where:
- $\text{FASTQ}_j$ = FASTQ files for sample $j$
- $R_{:,j}$ = read counts for sample $j$
- $S = \{S_1, \ldots, S_N\}$ = transcript sequences
- $\Theta$ = simulation parameters (fragment, error, bias models)

The function $\mathcal{F}_j$ implements:
1. Fragment generation with biases
2. Read extraction with strand orientation
3. Error introduction
4. FASTQ formatting

This mathematical representation provides a complete description of the bioinformatic processes simulated by Polyester from count matrix input to FASTQ output, enabling precise control over RNA-seq simulation experiments for method evaluation and benchmarking.

## References

1. Frazee et al. (2014). Polyester: simulating RNA-seq datasets with differential transcript expression. *Bioinformatics*.
2. Li & Jiang (2012). Transcriptome assembly and isoform expression level estimation from biased RNA-Seq reads. *Bioinformatics* 28(22): 2914-2921.
3. McElroy et al. (2012). GemSIM: general, error-model based simulator of next-generation sequencing data. *BMC Genomics* 13(1): 74.
4. 't Hoen et al. (2013). Reproducibility of high-throughput mRNA and small RNA sequencing across laboratories. *Nature Biotechnology* 31(11): 1015-1022.