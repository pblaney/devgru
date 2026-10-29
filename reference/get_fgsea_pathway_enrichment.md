# Run fgsea gene set enrichment analysis on DESeq2 differential gene expression output

Run fgsea pathway enrichment workflow from DESeq2 DGE analysis with
[`get_deseq2_diff_expr()`](https://pblaney.github.io/devgru/reference/get_deseq2_diff_expr.md)
using pathway gene set GMT file of interest.

## Usage

``` r
get_fgsea_pathway_enrichment(
  deseq2_diff_expr_output,
  pathways,
  cpus = 1,
  verbose = T
)
```

## Arguments

- deseq2_diff_expr_output:

  DESeq2 object, works directly from
  [`get_deseq2_diff_expr()`](https://pblaney.github.io/devgru/reference/get_deseq2_diff_expr.md)
  output

- pathways:

  GMT file with desired pathway gene sets to test for enrichment

- cpus:

  number of cpus for reading in data, used by
  [`fgsea::fgseaMultilevel()`](https://rdrr.io/pkg/fgsea/man/fgseaMultilevel.html),
  default: 1

- verbose:

  Output CLI during workflow, default: TRUE

## Value

list object with ranked genes vector and enrichment results data.table

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
# get_fgsea_pathway_enrichment(
# deseq2_diff_expr_output = dge_demo,
# pathways = "c7.bcell.immunesigdb.v2025.1.Hs.symbols.gmt",
# cpus = 2)
```
