% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/swissArmyKnife.R
\name{get_cbs_per_chromosome}
\alias{get_cbs_per_chromosome}
\title{Single chromosome run circular binary segmentation (CBS) algorithm on read depth profile}
\usage{
get_cbs_per_chromosome(
  chrom_for_cbs,
  chromosome_names,
  signal,
  position,
  sample_id
)
}
\arguments{
\item{chrom_for_cbs}{The chromosome to run CBS algorithm on}

\item{chromosome_names}{A vector of chromosome strings, equivalent to the seqnames
column. Must contain at least the \code{chrom_for_cbs}}

\item{signal}{A vector of read depth signal to be used as input, equivalent to
the foreground or read ratio column. These values will be \code{log}ged before use in CBS}

\item{position}{A vector of the genomic position of the signal measurement, equivalent
to the start column}

\item{sample_id}{A string used to populate the ID column in the data.table}
}
\value{
data.table of CBS segmentation values per chromosome
}
\description{
Run CBS algorithm using DNAcopy on a per-chromosome basis on read depth profile. This
allows for parallel computation of the segmentation, rapidly reducing run time.
There should be no \code{NA}s in the columns or \code{0}s in the \code{signal} column.

Note, this function was designed to be run as part of the \code{get_dryclean_segmentation()}
workflow. Also, see \code{chrompar()} for executing in parallel.
}
\examples{
# For this example, there are no NAs or zeros in the `signal`
# Check for these and may need to filter them out
read_depth_demo_hg38

# Example of what the input looks like for
# chromosome_names
head(as.character(GenomicRanges::seqnames(read_depth_demo_hg38)))
# signal
head(as.double(GenomicRanges::values(read_depth_demo_hg38)[, "foreground"]))
# position
head(GenomicRanges::start(read_depth_demo_hg38))

}
\keyword{workflow}
