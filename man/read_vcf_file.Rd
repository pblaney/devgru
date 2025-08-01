% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/swissArmyKnife.R
\name{read_vcf_file}
\alias{read_vcf_file}
\title{Read in a VCF file and convert to GenomicRanges object}
\usage{
read_vcf_file(
  vcf_file_path,
  tumor_sample = NULL,
  normal_sample = NULL,
  caller = NULL,
  mut_type = NULL,
  seq_lengths = gUtils::hg_seqlengths()
)
}
\arguments{
\item{vcf_file_path}{Path to VCF file}

\item{tumor_sample}{Name of tumor sample as reported in VCF}

\item{normal_sample}{Name of normal sample as reported in VCF}

\item{caller}{Name of the caller that generated input VCF to be converted, supported: mutect, strelka, varscan, svaba, caveman, deepvariant}

\item{mut_type}{Type of mutations within VCF, supported: snv, indel, snp}

\item{seq_lengths}{Named vector object used as the template for new seq details, see \code{gUtils::hg_seqlengths()} for example}
}
\value{
GenomicRanges object with VCF FILTER/INFO/FORMAT columns, if present, and updated seqinfo, seqnames, seqlengths, seqlevels
}
\description{
Read in a VCF file and convert it to a GenomicRanges object with refactored seq details.
Currently supports conversion of somatic SNV/InDel VCFs from Mutect, Strelka, Varscan, SvABA, CaVEMan.
Also supports conversion of germline SNP/InDel VCF from DeepVariant.
The VCF file can be either zipped or unzipped.
}
