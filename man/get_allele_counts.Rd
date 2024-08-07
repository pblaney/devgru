% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/swissArmyKnife.R
\name{get_allele_counts}
\alias{get_allele_counts}
\title{Read in a BAM file, collect read counts for alleles at set of mutation loci using alleleCounter, and convert to data.table object}
\usage{
get_allele_counts(
  bam_file_path,
  mut_loci_obj,
  min_base_qual = 20,
  min_map_qual = 35,
  threads = 1
)
}
\arguments{
\item{bam_file_path}{Path to BAM file}

\item{mut_loci_obj}{Mutation loci file in data.table or GRanges format, required columns: chrom, pos, ref, alt}

\item{min_base_qual}{Minimum base quality required for a read to be counted, default: 20}

\item{min_map_qual}{Minimum mapping quality required for a read to be counted, default: 35}

\item{threads}{Number of threads to use for parallel execution, default: 1}
}
\value{
data.table object with the columns: #CHR, POS, Count_A, Count_C, Count_G, Count_T, Good_depth
}
\description{
Wrapper command to use alleleCounter to read in a BAM file and count all reads supporting each possible allele at set of mutation loci, then convert it to a data.table object with refactored seq details.
The process will be split across each chromosome and then merged for the final output.
Expects the first column to be chromosome name, and 3 other columns to exist in the mutation loci file: pos/start/end, ref/REF/Reference_Allele, alt/ALT/Tumor_Seq_Allele2.
}
