# Package index

## Environment setup

- [`kit_loadout()`](https://pblaney.github.io/devgru/reference/kit_loadout.md)
  : Build the devgru environment by installing/loading packages

- [`chrompar()`](https://pblaney.github.io/devgru/reference/chrompar.md)
  :

  Scaffold function around `easypar` to run workflows in parallel by
  chromosome

## Readers

- [`read_bed_file()`](https://pblaney.github.io/devgru/reference/read_bed_file.md)
  : Read in a BED file and convert to GenomicRanges or data.table
- [`read_bedpe_file()`](https://pblaney.github.io/devgru/reference/read_bedpe_file.md)
  : Read BEDPE File with Flexible Output Format
- [`read_gtf_file()`](https://pblaney.github.io/devgru/reference/read_gtf_file.md)
  : Read in a GTF file and convert to GenomicRanges or data.table
- [`read_maf_file()`](https://pblaney.github.io/devgru/reference/read_maf_file.md)
  : Read MAF file and convert to GenomicRanges or data.table
- [`read_vcf_file()`](https://pblaney.github.io/devgru/reference/read_vcf_file.md)
  : Read in a VCF file and convert to GenomicRanges or data.table

## Conversions

- [`dt_to_gr()`](https://pblaney.github.io/devgru/reference/dt_to_gr.md)
  : Convert data.table to GenomicRanges Object
- [`dt_to_gi()`](https://pblaney.github.io/devgru/reference/dt_to_gi.md)
  : Convert BEDPE-like data.table to GInteractions
- [`gi_to_dt()`](https://pblaney.github.io/devgru/reference/gi_to_dt.md)
  : Convert GInteractions to BEDPE data.table format
- [`gi_to_grl()`](https://pblaney.github.io/devgru/reference/gi_to_grl.md)
  : Convert GInteractions to GRangesList
- [`grl_to_gi()`](https://pblaney.github.io/devgru/reference/grl_to_gi.md)
  : Convert GenomicRangesList to GInteractions
- [`gr_to_seg()`](https://pblaney.github.io/devgru/reference/gr_to_seg.md)
  : Convert segmentation GenomicRanges object to a data.table with
  standardized columns and format

## GenomicRanges utilities

- [`gr_flank()`](https://pblaney.github.io/devgru/reference/gr_flank.md)
  : Precise addition of flanking to GenomicRanges with boundary-aware
  options
- [`gr_sanitycheck()`](https://pblaney.github.io/devgru/reference/gr_sanitycheck.md)
  : Check if input is a GenomicRanges object, plus optional sanity check
  of column names
- [`dt_sanitycheck()`](https://pblaney.github.io/devgru/reference/dt_sanitycheck.md)
  : Check if input is a data.table object, plus optional sanity check of
  column names
- [`gi_refactor_seqs()`](https://pblaney.github.io/devgru/reference/gi_refactor_seqs.md)
  : Refactor seqinfo of GInteractions object for easy harmony
- [`gr_refactor_seqs()`](https://pblaney.github.io/devgru/reference/gr_refactor_seqs.md)
  : Refactor seqinfo, seqnames, seqlengths, seqlevels of GenomicRanges
  object for easy harmony
- [`grl_refactor_seqs()`](https://pblaney.github.io/devgru/reference/grl_refactor_seqs.md)
  : Refactor seqinfo of GRangesList object for easy harmony
- [`aggregate_these()`](https://pblaney.github.io/devgru/reference/aggregate_these.md)
  : Read and aggregate genomic or non-genomic data files with smart
  format detection

## Analysis

- [`get_cbs_per_chromosome()`](https://pblaney.github.io/devgru/reference/get_cbs_per_chromosome.md)
  : Single chromosome run circular binary segmentation (CBS) algorithm
  on read depth profile
- [`get_corrected_cnv_profile()`](https://pblaney.github.io/devgru/reference/get_corrected_cnv_profile.md)
  : Read in CNV profile data.table of various flavors and extract a
  corrected profile
- [`get_cov_gap_imputation()`](https://pblaney.github.io/devgru/reference/get_cov_gap_imputation.md)
  : Generate gap imputed coverage from fragCounter output
- [`get_deseq2_diff_expr()`](https://pblaney.github.io/devgru/reference/get_deseq2_diff_expr.md)
  : Run DESeq2 differential expression analysis on RNA-seq count data
- [`get_fgsea_pathway_enrichment()`](https://pblaney.github.io/devgru/reference/get_fgsea_pathway_enrichment.md)
  : Run fgsea gene set enrichment analysis on DESeq2 differential gene
  expression output
- [`get_fragcounter_segmentation()`](https://pblaney.github.io/devgru/reference/get_fragcounter_segmentation.md)
  : Generate a segmentation file from a fragCounter read depth profile
- [`get_igh_trx_subtype()`](https://pblaney.github.io/devgru/reference/get_igh_trx_subtype.md)
  : Detect Canonical IGH Translocations in Multiple Myeloma
- [`get_imputed_gaps_per_chromosome()`](https://pblaney.github.io/devgru/reference/get_imputed_gaps_per_chromosome.md)
  : Single chromosome run imputation of small, NA read gaps in
  fragCounter coverage
- [`get_jabba_qc_diagnostic()`](https://pblaney.github.io/devgru/reference/get_jabba_qc_diagnostic.md)
  : Run diagnostic workflow of after JaBbA run to help evaluate quality
  of graph
- [`get_maf_lite()`](https://pblaney.github.io/devgru/reference/get_maf_lite.md)
  : Convert SNV & InDel mutation table to MAF-lite format
- [`get_maf_nonsynonymous()`](https://pblaney.github.io/devgru/reference/get_maf_nonsynonymous.md)
  : Extract Non-Synonymous Mutations from MAF
- [`get_protein_coding_genes()`](https://pblaney.github.io/devgru/reference/get_protein_coding_genes.md)
  : Extract protein-coding genes from GTF file
- [`get_qc_diagnostics_alignment()`](https://pblaney.github.io/devgru/reference/get_qc_diagnostics_alignment.md)
  : Generate diagnostic plots for alignment QC checks using Alfred
  summary files
- [`get_qc_diagnostics_snvindel()`](https://pblaney.github.io/devgru/reference/get_qc_diagnostics_snvindel.md)
  : Generate diagnostic plots for SNV & InDel variant calling QC checks
- [`get_vaf()`](https://pblaney.github.io/devgru/reference/get_vaf.md) :
  Quick pull or explicitly calculate the VAF for mutation records of
  various flavors

## Plotting

- [`geom_deseq2_volcano()`](https://pblaney.github.io/devgru/reference/geom_deseq2_volcano.md)
  : Volcano plot for exploring results from DESeq2 differential gene
  expression analysis
- [`geom_gap_imputation()`](https://pblaney.github.io/devgru/reference/geom_gap_imputation.md)
  : Diagnostic plot used for reviewing the imputation of gaps in
  coverage profile
- [`geom_gsea_enrichment()`](https://pblaney.github.io/devgru/reference/geom_gsea_enrichment.md)
  : Enrichment plot for exploring results from fgsea pathway enrichment
  analysis
- [`geom_jabba_stats()`](https://pblaney.github.io/devgru/reference/geom_jabba_stats.md)
  : Diagnostic plot used for evaluating the performance of a JaBbA run
- [`copynumber_palettier()`](https://pblaney.github.io/devgru/reference/copynumber_palettier.md)
  : Create color palette for integer copy number states

## CLI helpers

- [`cli_pio_colorscheme()`](https://pblaney.github.io/devgru/reference/cli_pio_colorscheme.md)
  : Set the color scheme options for pio CLI text in devgru
- [`cli_stopwatch_end()`](https://pblaney.github.io/devgru/reference/cli_stopwatch_end.md)
  : Verbose CLI output to designate a workflow end
- [`cli_stopwatch_start()`](https://pblaney.github.io/devgru/reference/cli_stopwatch_start.md)
  : Verbose CLI output to designate a workflow start
- [`function_cli_intro()`](https://pblaney.github.io/devgru/reference/function_cli_intro.md)
  : Verbose CLI output for extended and rich description of steps along
  a complex workflow, typically used in a pipeline

## Datasets

- [`braf_dnase_demo_dt_hg19`](https://pblaney.github.io/devgru/reference/braf_dnase_demo_dt_hg19.md)
  : DNAaseI hypersensitivity sites for BRAF on hg19 as data.table
- [`braf_dnase_demo_dt_hg38`](https://pblaney.github.io/devgru/reference/braf_dnase_demo_dt_hg38.md)
  : DNAaseI hypersensitivity sites for BRAF on hg38 as data.table
- [`chromosome_arms_hg38`](https://pblaney.github.io/devgru/reference/chromosome_arms_hg38.md)
  : Genomic coordinates of each chromosome arm on hg38 as GenomicRanges
- [`example_dnase_hg38`](https://pblaney.github.io/devgru/reference/example_dnase_hg38.md)
  : DNAaseI hypersensitivity sites on hg38
- [`exclusion_regions_hg38`](https://pblaney.github.io/devgru/reference/exclusion_regions_hg38.md)
  : Collection of genomic regions to exclude from analysis on hg38 as
  GenomicRanges
- [`gene_body_hg38`](https://pblaney.github.io/devgru/reference/gene_body_hg38.md)
  : Gene body segments from Ensemble v108 on hg38 as GenomicRanges
- [`hmcl_conditions_demo_dt_hg38`](https://pblaney.github.io/devgru/reference/hmcl_conditions_demo_dt_hg38.md)
  : Condition labels for associated RNA-seq counts matrix as data.table
- [`hmcl_counts_demo_dt_hg38`](https://pblaney.github.io/devgru/reference/hmcl_counts_demo_dt_hg38.md)
  : RNA-seq counts matrix on hg38 as data.table
- [`kras_dnase_demo_gr_hg19`](https://pblaney.github.io/devgru/reference/kras_dnase_demo_gr_hg19.md)
  : DNAaseI hypersensitivity sites for KRAS on hg19 as GenomicRanges
- [`kras_dnase_demo_gr_hg38`](https://pblaney.github.io/devgru/reference/kras_dnase_demo_gr_hg38.md)
  : DNAaseI hypersensitivity sites for KRAS on hg38 as GenomicRanges
- [`myeloma_translocation_ig_loci_hg38`](https://pblaney.github.io/devgru/reference/myeloma_translocation_ig_loci_hg38.md)
  : Ig loci of recurrent translocations in multiple myeloma as
  GenomicRanges
- [`myeloma_translocation_partner_loci_hg38`](https://pblaney.github.io/devgru/reference/myeloma_translocation_partner_loci_hg38.md)
  : Partner loci of recurrent translocations in multiple myeloma as
  GenomicRanges
- [`read_depth_demo_hg38`](https://pblaney.github.io/devgru/reference/read_depth_demo_hg38.md)
  : dryclean-fragCounter read depth profile on hg38 as GenomicRanges
