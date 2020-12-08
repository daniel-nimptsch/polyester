# function to make sure everything is compatible and assign sane defaults
# internal

.check_extras = function(extras, paired, total.n){

    if(!('distr' %in% names(extras))){
        extras$distr = 'normal'
    }else{
        extras$distr = match.arg(extras$distr,
            c('normal', 'empirical', 'custom'))
        if(extras$distr == 'custom' & !('custdens' %in% names(extras))){
            stop(.makepretty('to use custom fragment distribution, provide
                "custdens", a logspline object representing the distribution.'))
        }
    }

    # I don't love this--fraglen and fragsd aren't needed unless distr is normal.
    # but we store them anyway. should code better?
    if (!('fraglen' %in% names(extras))) {
        extras$fraglen = rep(250, total.n)
    } else {
      if (length(extras$fraglen) == 1) {
        extras$fraglen = rep(extras$fraglen, total.n)
      } else {
        stopifnot(length(extras$fraglen) == total.n)
      }
    }
    if (!('fragsd' %in% names(extras))) {
        extras$fragsd = rep(25, total.n)
    } else {
      if (length(extras$fragsd) == 1) {
        extras$fragsd = rep(extras$fragsd, total.n)
      } else {
        stopifnot(length(extras$fragsd) == total.n)
      }
    }

    if(!('readlen' %in% names(extras))){
        extras$readlen = 100
    }

    if(!('bias' %in% names(extras))){
        extras$bias = 'none'
    }else{
        extras$bias = match.arg(extras$bias, c('none', 'rnaf', 'cdnaf'))
    }

    if(!('error_model' %in% names(extras))){
        extras$error_model = 'uniform'
    }
    .check_error_model(extras, paired)

    if(!('error_rate' %in% names(extras))){
        extras$error_rate = 0.005
    }
    if(extras$error_model == 'custom'){
        extras$path = paste0(extras$model_path, '/', extras$model_prefix)
    }#this should work beause we already checked stuff.

    if(!('bias' %in% names(extras))){
        extras$bias = 'none'
    }else{
        extras$bias = match.arg(extras$bias, c('none', 'rnaf', 'cdnaf'))
    }

    if(!('lib_sizes' %in% names(extras))){
        extras$lib_sizes = rep(1, total.n)
    }else{
        stopifnot(is.numeric(extras$lib_sizes))
        stopifnot(length(extras$lib_sizes) == total.n)
    }

    if (!('frag_GC_bias' %in% names(extras))) {
      extras$frag_GC_bias <- 'none'
    } else {
      stopifnot(is.matrix(extras$frag_GC_bias))
      stopifnot(nrow(extras$frag_GC_bias) == 101)
      stopifnot(ncol(extras$frag_GC_bias) == total.n)
      stopifnot(all(extras$frag_GC_bias >= 0 & extras$frag_GC_bias <= 1))
    }

    if (!('strand_specific' %in% names(extras))) {
      extras$strand_specific <- FALSE
    }

    if (!('shuffle' %in% names(extras))) {
       extras$shuffle <- FALSE
    }
    if (!('fastq' %in% names(extras))) {
      extras$fastq <- FALSE
    }
    if (!('verbose' %in% names(extras))) {
      extras$verbose <- FALSE
    }
    if ('seq_depth' %in% names(extras)) {
      stopifnot(is.numeric(extras$seq_depth))
      if(length(extras$seq_depth) == 1)
        extras$seq_depth <- rep(extras$seq_depth, total.n)
      else 
        stopifnot(length(extras$seq_depth) == total.n)
    }
    # for alternative splicing simulator:
    if (!('exon_junction_coverage' %in% names(extras))){
      extras$exon_junction_coverage <- FALSE
    } else {
      if (!('exon_junction_table' %in% names(extras))) {
        stop("to use exon junction coverage please provide an exon junction table. This option should only be called inside of 'simulate_alternative_splicing'.")
      } else {
        if (!(data.table::is.data.table(extras$exon_junction_table))) {
          stop("please provide an exon junction table as data.table. This option should only be called inside of 'simulate_alternative_splicing'.")
        } else {
          if (any(!(c('transcript_id', 'type', 'tr_start', 'tr_end') %in% names(extras$exon_junction_table))))
            stop("please provide columns 'transcript_id', 'type', 'tr_start', 'tr_end' in exon junction table. This option should only be called inside of 'simulate_alternative_splicing'.")
          else
          extras$exon_junction_table = extras$exon_junction_table[, ASS_ID := as.character(.I)]
        }
      }
    }
    if ('adapter_contamination' %in% names(extras)) {
      stopifnot(is.logical(extras$adapter_contamination))
      if (extras$adapter_contamination) {
        if ('adapter_sequence' %in% names(extras)){
          stopifnot(is.character(extras$adapter_sequence))
        } else 
          extras$adapter_sequence <- 'CTGTCTCTTATACACATCT'
      }
    } else 
      extras$adapter_contamination <- FALSE
    if ('pcr_rate' %in% names(extras)) {
      stopifnot(extras$pcr_rate >= 0 & extras$pcr_rate <= 1)
      if ('pcr_lambda' %in% names(extras)) {
        stopifnot(is.numeric(extras$pcr_lambda))
      } else
        extras$pcr_lambda <- 1
    } else 
      extras$pcr_rate <- NULL
    return(extras)

}
