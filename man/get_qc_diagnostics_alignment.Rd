% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/swissArmyKnife.R
\name{get_qc_diagnostics_alignment}
\alias{get_qc_diagnostics_alignment}
\title{Generate diagnostic plots for alignment QC checks using Alfred summary files}
\usage{
get_qc_diagnostics_alignment(
  path_to_tumor_dir = NULL,
  path_to_normal_dir = NULL,
  dataset_title = "",
  seq_center = "",
  seq_protocol = "WGS",
  output_path = getwd()
)
}
\arguments{
\item{path_to_tumor_dir}{Path to directory of tumor sample Alfred summary files}

\item{path_to_normal_dir}{Path to directory of normal sample Alfred summary files}

\item{dataset_title}{Name of the dataset}

\item{seq_center}{Institute/site where sequencing was performed}

\item{seq_protocol}{Type of sequencing protocol for display purposes, default: WGS}

\item{output_path}{Location to save report PDF, default: \code{getwd()}}
}
\value{
Patchwork 'quilt'-like plot of ggplots
}
\description{
Single command to generate a comprehensive set of diagnostic plots that assist in
performing QC checks of tumor/normal read alignment of various sequencing protocol flavors
using Alfred summary files.
}
\examples{
# After downloading `*.qc.summary.txt` Alfred output file from the MGP1000 run,
# place tumor and normal sample files into respective subdirectories
# get_qc_diagnostics_alignment(
# path_to_tumor_dir = "data/sanger/examplePDexport_999/preprocessing/alfred/tumor/",
# path_to_normal_dir = "data/sanger/examplePDexport_999/preprocessing/alfred/normal/",
# dataset_title = "examplePDexport_999",
# seq_center = "Wellcome Sanger Institute",
# output_path = "data/sanger/examplePDexport_999/preprocessing/")

}
\keyword{workflow}
