#' @name cdnaf
#' @title
#'   Model of positional bias that can arise when RNA-seq is performed
#'   using protocols relying on cDNA fragmentation.
#' @description
#'   This positional bias model was estimated in Li and Jiang (2012).
#'   With cDNA fragmentation, reads are more likely to have come from the 3'
#'   end of the transcript. The probabilities included in this dataset were
#'   estimated from Supplementary Figure S3 in Li and Jiang's manuscript.
#'   Data points from the figure were inferred and exported as CSV files
#'   using WebPlotDigitizer. The CSV files and the code used to process them
#'   and create the datasets are available in the Polyester GitHub
#'   repository (\url{https://github.com/alyssafrazee/polyester}).
#' @docType data
#' @format
#'   data frame with 100 rows and 2 columns. Column 1 is position along a
#'   transcript (in percent), while Column 2 is the probability of getting a
#'   fragment at that position. Column 2 sums to 1.
#' @references
#'   Li W and Jiang T (2012): Transcriptome assembly and isoform expression
#'   level estimation from biased RNA-Seq reads. Bioinformatics 28(22):
#'   2914-2921.
#'
#'   Rohatgi A (2014): WebPlotDigitizer: Version 3.4 of WebPlotDigitizer.
#'   ZENODO. 10.5281/zenodo.11835
NULL
