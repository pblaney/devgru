% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/swissArmyKnife.R
\name{geom_deseq2_volcano}
\alias{geom_deseq2_volcano}
\title{Volcano plot for exploring results from DESeq2 differential gene expression analysis}
\usage{
geom_deseq2_volcano(
  deseq2_diff_expr_output,
  log2_fc_threshold = 2,
  p_adj_threshold = 0.05,
  n_positive_labels = 25,
  n_negative_labels = 25,
  label_vector = NULL,
  point_label_nudge_x = 0,
  point_label_nudge_y = 0,
  condition_label_x = 5
)
}
\arguments{
\item{deseq2_diff_expr_output}{DESeq2 object, works directly from \code{get_deseq2_diff_expr()} output}

\item{log2_fc_threshold}{Numeric cutoff for coloring points based on log fold change, default: 2}

\item{p_adj_threshold}{Numeric cutoff for coloring points based on adjust P-value significance, default: 0.05}

\item{n_positive_labels}{Number of top genes to label on positive side of log fold change, default: 25}

\item{n_negative_labels}{Number of top genes to label on negative side of log fold change, default: 25}

\item{label_vector}{Set of a priori genes to label, cannot be used in tandem with top N genes}

\item{point_label_nudge_x}{Distance on x-axis to adjust the gene labels, can be positive or negative value, default: 0}

\item{point_label_nudge_y}{Distance on y-axis to adjust the gene labels, can be positive or negative value, default: 0}

\item{condition_label_x}{Location on x-axis to put condition labels, default: 5}
}
\value{
ggplot object
}
\description{
Creates a ggplot2 style volcano plot for exploring the results from a DESeq2 DGE
analysis with parameters for easy styling
}
\examples{
# Run the DESeq2 DGE workflow
# dge_demo <- get_deseq2_diff_expr(
# counts_file_path = fs::path_package("extdata", "hmcl_counts.hg38.txt.gz", package = "devgru"),
# condition_file_path = fs::path_package("extdata", "hmcl_conditions.hg38.txt", package = "devgru"),
# gene_universe = gene_body_hg38,
# gene_symbol_column = "#GENE",
# ensembl_id_column = NULL,
# formula_string = "U266 ~ RPMI8226",
# min_transcripts = 10)

# Exploratory plot with top 25 genes labeled on both ends of spectrum
# geom_deseq2_volcano(deseq2_diff_expr_output = dge_demo)

# Targeted assessment with labels for genes of interest
# geom_deseq2_volcano(deseq2_diff_expr_output = dge_demo, label_vector = c("CCND1","MYC","MYCL"))

}
\keyword{workflow}
