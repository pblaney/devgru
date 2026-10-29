# Enrichment plot for exploring results from fgsea pathway enrichment analysis

Creates a ggplot2 style enrichment plot for exploring the results from a
fgsea GSEA analysis.

## Usage

``` r
geom_gsea_enrichment(pathways, pathway_of_interest, ranked_genes)
```

## Arguments

- pathways:

  GMT file with desired pathway gene sets to test for enrichment

- pathway_of_interest:

  Specific pathway to in GMT to explore

- ranked_genes:

  Named vector of genes and their ranks, works directly from
  [`get_fgsea_pathway_enrichment()`](https://pblaney.github.io/devgru/reference/get_fgsea_pathway_enrichment.md)
  output

## Value

ggplot object

## Examples

``` r
# Run the DESeq2 DGE workflow
# dge_demo <- get_deseq2_diff_expr(
# counts_file_path = fs::path_package("extdata", "hmcl_counts.hg38.txt.gz", package = "devgru"),
# condition_file_path = fs::path_package("extdata", "hmcl_conditions.hg38.txt", package = "devgru"),
# gene_universe = gene_body_hg38,
# gene_symbol_column = "#GENE",
# ensembl_id_column = NULL,
# formula_string = "U266 ~ RPMI8226",
# min_transcripts = 10)

# Get enrichment against B cell specific pathways from ImmuneSigDB
# gsea_demo <- get_fgsea_pathway_enrichment(
# deseq2_diff_expr_output = dge_demo,
# pathways = "c7.bcell.immunesigdb.v2025.1.Hs.symbols.gmt",
# cpus = 2)

# Plot the enrichment of specific pathway
# geom_gsea_enrichment(
# pathways = "~/morganLab/nsd2dTagSystem/data/c7.bcell.immunesigdb.v2025.1.Hs.symbols.gmt",
# pathway_of_interest = "GSE38304_MYC_NEG_VS_POS_GC_BCELL_UP",
# ranked_genes = gsea_demo$ranked_genes)
```
