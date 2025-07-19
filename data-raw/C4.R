# }}}}------->>> Code used to generate the internal data used for demos.

library(devgru)
library(gUtils)
library(dplyr)
library(data.table)
library(fs)

# library(BSgenome.Hsapiens.UCSC.hg19)
# library (rtracklayer)

#
# }}}}------->>> GenomicRanges Demos
#

# # From https://github.com/mskilab-org/gUtils/blob/master/data-raw/prep_UCSC.R
# # prepare a Seqinfo object
# sl <- GenomeInfoDb::seqlengths(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)[1:25]
# si <- GenomeInfoDb::Seqinfo(seqnames = names(sl), seqlengths = sl, genome = "hg19")
# GenomeInfoDb::isCircular(si) <- rep(FALSE, length(sl))
#
# # download data from UCSC Table Browser to include as example data in gUtils
# mySession <- rtracklayer::browserSession("UCSC")
# rtracklayer::genome(mySession) <- "hg19"
#
# # DNAaseI hypersensitivity sites
# track.name <- "wgEncodeUwDgf"
# table.name <- "wgEncodeUwDgfK562Hotspots"
# tbl.DNAase <- rtracklayer::getTable( ucscTableQuery(mySession, track = track.name, table = table.name))
# example_dnase <- with(tbl.DNAase, GenomicRanges::GRanges(gsub("chr(.*?)", "\\1", chrom), IRanges::IRanges(chromStart, chromEnd), seqinfo = si))
# mcols(example_dnase) <- tbl.DNAase[,-which(colnames(tbl.DNAase) %in% c("chrom", "chromStart","chromEnd","strand","name","qValue","bin","score"))]

# hg38
# Original data is hg19 so let's liftOver
liftchain_hg19_to_hg38 <- rtracklayer::import.chain("~/morganLab/mgp1000/development/nonCoding/hg19ToHg38.over.chain")
example_dnase_hg38 <- rtracklayer::liftOver(x = gUtils::gr.chr(gUtils::example_dnase), chain = liftchain_hg19_to_hg38)
# Unlist output to get GR
example_dnase_hg38 <- gUtils::grl.unlist(example_dnase_hg38)
# Clean-up footprint of example
example_dnase_hg38 <- example_dnase_hg38[,1:2]
# Export for usage
usethis::use_data(example_dnase_hg38, internal = FALSE, overwrite = TRUE)


# Use example DNase data from gUtils and find overlap with hg19 GR object of KRAS and 1000 bp flank
# hg19
kras_dnase_demo_gr_hg19 <- gUtils::gr.findoverlaps(query = gUtils::example_dnase,
                                                   subject = gUtils::gr.nochr(gUtils::parse.gr("chr12:25358180-25403863")) + 1000,
                                                   qcol = c("signalValue",'pValue'))
# Add gene name and biospecimen column
kras_dnase_demo_gr_hg19$biospecimen <- "K562"
kras_dnase_demo_gr_hg19$gene <- "KRAS"
# Clean-up footprint of example
kras_dnase_demo_gr_hg19 <- kras_dnase_demo_gr_hg19[,3:6]
# And complete the remaining seqinfo columns for genome and isCircular
GenomeInfoDb::genome(kras_dnase_demo_gr_hg19) <- "GRCh37"
GenomeInfoDb::isCircular(kras_dnase_demo_gr_hg19) <- rep(FALSE,25)
# Export for usage
usethis::use_data(kras_dnase_demo_gr_hg19, internal = FALSE, overwrite = TRUE)

# hg38
kras_dnase_demo_gr_hg38 <- gUtils::gr.findoverlaps(query = example_dnase_hg38,
                                                   subject = gUtils::parse.gr("chr12:25205246-25250929") + 1000,
                                                   qcol = c("signalValue",'pValue'))
# Add gene name column
kras_dnase_demo_gr_hg38$biospecimen <- "K562"
kras_dnase_demo_gr_hg38$gene <- "KRAS"
# Clean-up footprint of example and seqinfo
kras_dnase_demo_gr_hg38 <- kras_dnase_demo_gr_hg38[,3:6]
kras_dnase_demo_gr_hg38 <- gr_refactor_seqs(kras_dnase_demo_gr_hg38)
# Export for usage
usethis::use_data(kras_dnase_demo_gr_hg38, internal = FALSE, overwrite = TRUE)

#
# }}}}------->>> data.table Demos
#

# hg19
braf_dnase_demo_dt_hg19 <- gUtils::gr.findoverlaps(query = gUtils::example_dnase,
                                                   subject = gUtils::gr.nochr(gUtils::parse.gr("chr7:140419137-140624729")) + 1000,
                                                   qcol = c("signalValue",'pValue'),
                                                   return.type = "data.table") %>%
                                                dplyr::select(seqnames,start,end,strand,signalValue,pValue)
# Add gene name and biospecimen column
braf_dnase_demo_dt_hg19$biospecimen <- "K562"
braf_dnase_demo_dt_hg19$gene <- "KRAS"
# Export for usage
usethis::use_data(braf_dnase_demo_dt_hg19, internal = FALSE, overwrite = TRUE)

# hg38
braf_dnase_demo_dt_hg38 <- gUtils::gr.findoverlaps(query = example_dnase_hg38,
                                                   subject = gUtils::parse.gr("chr7:140719337-140924929") + 1000,
                                                   qcol = c("signalValue",'pValue'),
                                                   return.type = "data.table") %>%
                                                dplyr::select(seqnames,start,end,strand,signalValue,pValue)
# Add gene name and biospecimen column
braf_dnase_demo_dt_hg38$biospecimen <- "K562"
braf_dnase_demo_dt_hg38$gene <- "KRAS"
# Export for usage
usethis::use_data(braf_dnase_demo_dt_hg38, internal = FALSE, overwrite = TRUE)


#
#
# }}}}------->>> Exclusion list comprised of UCSC and ENCODE blacklists from excluderanges (https://github.com/dozmorovlab/excluderanges)
#
#

# The following code was used to generate the file exclusion_regions.bed stored in the data/ directory
# # AnnotationHub has some compatibility issues with newer versions of dbplyr
# dbplyr_v2.3.0 <- "https://github.com/tidyverse/dbplyr/archive/refs/tags/v2.3.1.tar.gz"
# install.packages(dbplyr_v2.3.0, repos = NULL, type = "source")
# library(dbplyr)
# library(AnnotationHub)
#
# ah <- AnnotationHub::AnnotationHub()
# query_data <- AnnotationHub::subset(ah, preparerclass == "excluderanges")
#
# # telomeres, extracted from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gaps.txt.gz
# hg38.UCSC.telomere <- query_data[["AH95938"]]
# telomeres <- gUtils::gr2dt(hg38.UCSC.telomere) %>%
#   dplyr::rename("start_0based" = start, "exclusion_segment" = type) %>%
#   dplyr::mutate("start" = dplyr::case_when(start_0based == 0 ~ 1,
#                                            start_0based != 0 ~ start_0based)) %>%
#   dplyr::select(seqnames,start,end,exclusion_segment)
#
# # centromeres, extracted from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gaps.txt.gz
# centromeres <- data.table::fread(file = "~/Downloads/centromeres.txt",
#                                  sep = "\t",
#                                  header = F) %>%
#   dplyr::select(V2,V3,V4) %>%
#   dplyr::rename("seqnames" = V2, "start" = V3, "end" = V4) %>%
#   dplyr::mutate("exclusion_segment" = "centromere")
#
# # acrocentric short arm chromosomes
# hg38.UCSC.short_arm <- query_data[["AH95937"]]
# acrocentric <- gUtils::gr2dt(hg38.UCSC.short_arm) %>%
#   dplyr::mutate("exclusion_segment" = "acrocentric_short_arm") %>%
#   dplyr::select(seqnames,start,end,exclusion_segment)
#
# # heterochromatin
# hg38.UCSC.heterochromatin <- query_data[["AH95935"]]
# heterochromatin <- gUtils::gr2dt(hg38.UCSC.heterochromatin) %>%
#   dplyr::mutate("exclusion_segment" = "heterochromatin") %>%
#   dplyr::select(seqnames,start,end,exclusion_segment)
#
# # DAC Exclusion List Regions https://www.encodeproject.org/annotations/ENCSR636HFF/
# hg38.Kundaje.GRCh38_unified_Excludable <- query_data[["AH95917"]]
# blacklist <- gUtils::gr2dt(hg38.Kundaje.GRCh38_unified_Excludable) %>%
#   dplyr::mutate("exclusion_segment" = "dac_blacklist") %>%
#   dplyr::select(seqnames,start,end,exclusion_segment)

# Don't rerun the code snippet but simply load the data here
exclusion_regions_hg38 <- read_bed_file(bed_file_path = fs::path_package("extdata", "exclusion_regions.hg38.bed", package = "devgru"),
                                        has_header = T)
# Export for usage
usethis::use_data(exclusion_regions_hg38, internal = FALSE, overwrite = TRUE)

#
#
# }}}}------->>> Chromosome arm coordinates comprised of genomic span from each start to per arm side of the centromere
#
#

# Don't rerun the code snippet but simply load the data here
chromosome_arms_hg38 <- read_bed_file(bed_file_path = fs::path_package("extdata", "chromosome_arms.hg38.bed", package = "devgru"),
                                 has_header = T)
# Export for usage
usethis::use_data(chromosome_arms_hg38, internal = FALSE, overwrite = TRUE)

#
#
# }}}}------->>> dryclean-fragCounter read depth profile demo
#
#

# # Read in dryclean-fragCounter read depth profile
# dryclean_fragcounter_cbs_demo <- readRDS(file = "~/morganLab/mgp1000/development/jAbBa/01-076_B_S10_001.dryclean.fragcounter.cov.rds")
# # Extract a subset of data from chr19
# read_depth_demo_chr19 <- gUtils::gr.findoverlaps(query = dryclean_fragcounter_cbs_demo,
#                                                  subject = gUtils::parse.gr("chr19:23000000-25000000"),
#                                                  qcol = c("background.log","foreground.log","input.read.counts",
#                                                           "median.chr","foreground","background","log.reads",
#                                                           "germline.status"),
#                                                  return.type = "data.table") %>%
#   dplyr::select(seqnames,start,end,background.log,foreground.log,
#                 input.read.counts,median.chr,foreground,background,
#                 log.reads,germline.status)
# # Extract a subset of data from chr22
# read_depth_demo_chr22 <- gUtils::gr.findoverlaps(query = dryclean_fragcounter_cbs_demo,
#                                                  subject = gUtils::parse.gr("chr21:23000000-25000000"),
#                                                  qcol = c("background.log","foreground.log","input.read.counts",
#                                                           "median.chr","foreground","background","log.reads",
#                                                           "germline.status"),
#                                                  return.type = "data.table") %>%
#   dplyr::select(seqnames,start,end,background.log,foreground.log,
#                 input.read.counts,median.chr,foreground,background,
#                 log.reads,germline.status)
# # Combine the chr19 and chr22 data for parallel demo as well
# read_depth_demo_hg38 <- gUtils::rrbind(read_depth_demo_chr19,
#                                        read_depth_demo_chr22)

# Don't rerun the code snippet but simply load the data here
read_depth_demo_hg38 <- read_bed_file(bed_file_path = fs::path_package("extdata", "read_depth_profile.hg38.bed.gz", package = "devgru"),
                                      has_header = T)
# Export for usage
usethis::use_data(read_depth_demo_hg38, internal = FALSE, overwrite = TRUE)

#
#
# }}}}------->>> Representative Transcripts for Genes
#
#

# # Load full gene annotation GTF
# full_genes_gtf <- devgru::read_gtf_file(gtf_file_path = paste0(base_data_path, "Homo_sapiens.GRCh38.108.gtf.gz"))
# full_genes_gtf_dt <- gr2dt(full_genes_gtf)
#
# # Filter to the same set of gene biotypes: protein-coding, long non-coding, and micro RNA
# # Must have at least 1 transcript record
# # Also filter out the records for selenocysteine, start codon, stop codon, and CDS as these regions are
# # covered within the exons
# genes_of_interest <- full_genes_gtf_dt %>%
#   dplyr::filter(gene_biotype %in% c("IG_C_gene","IG_D_gene","IG_J_gene","IG_V_gene",
#                                     "lncRNA","miRNA","protein_coding","TR_C_gene","TR_D_gene",
#                                     "TR_J_gene","TR_V_gene") &
#                   !is.na(gene_name) &
#                   !is.na(transcript_id) &
#                   !type %in% c("Selenocysteine","start_codon","stop_codon","CDS"))
#
# # The MANE Select and Plus Clinical transcripts will be used for protein coding genes
# # https://www.ncbi.nlm.nih.gov/refseq/MANE/
# mane_transcripts_full_summary <- data.table::fread(file = paste0(base_data_path, "MANE.GRCh38.v1.4.summary.txt.gz"))
#
# # Convert Ensembl string to match GTF and select single transcript for all genes
# # To construct the final set of transcripts, the first MANE Plus Clinical transcript will be selected if both are present (n=66)
# # For this version, the gene GNAS has 2 MANE Plus Clinical transcripts, in this case the first record will be kept (verified in UCSC browser this is most representative)
# mane_clinical_plus_transcripts <- mane_transcripts_full_summary %>%
#   dplyr::filter(MANE_status == "MANE Plus Clinical")
# mane_clinical_plus_transcripts <- mane_clinical_plus_transcripts[!duplicated(mane_clinical_plus_transcripts$symbol)]
#
# # Checked for duplicates, none present
# mane_select_transcripts <- mane_transcripts_full_summary %>%
#   dplyr::filter(!symbol %in% mane_clinical_plus_transcripts$symbol)
#
# mane_transcripts <- rbind(mane_clinical_plus_transcripts, mane_select_transcripts) %>%
#   dplyr::mutate("Ensembl_transcript" = str_remove(string = Ensembl_nuc, pattern = "\\..*"))
#
# # Grab the transcript name that corresponds to the gene annotation
# # This match is done using the width of the gene record and transcript record
# representative_transcript <- c()
# for(i in 1:n_distinct(genes_of_interest$gene_name)) {
#   # Isolate the available Ensembl transcript IDs for each gene
#   gene_by_transcript <- genes_of_interest %>%
#     dplyr::filter(gene_name == unique(genes_of_interest$gene_name)[i] & type == "transcript")
#
#   # Check for the MANE transcript
#   gene_by_mane_transcript <- mane_transcripts %>%
#     dplyr::filter(symbol == unique(genes_of_interest$gene_name)[i])
#
#   # Prioritize the MANE transcript if present
#   if(nrow(gene_by_mane_transcript) > 0) {
#     representative_transcript[i] <- gene_by_mane_transcript$Ensembl_transcript
#
#   } else {
#     # otherwise, choose the longest
#     representative_transcript[i] <- genes_of_interest %>%
#       dplyr::filter(gene_name == unique(genes_of_interest$gene_name)[i] & type == "transcript") %>%
#       dplyr::arrange(desc(width)) %>%
#       dplyr::select(transcript_id) %>%
#       dplyr::first() %>%
#       as.character()
#   }
#
#   if(i %% 1000 == 0) {
#     message(paste0(i, " representative transcripts selected ..."))
#   }
# }
#
# # Subset each GR obj to only include the representative transcripts
# per_representative_transcript_gr <- gr_refactor_seqs(dt2gr(genes_of_interest %>%
#                                                              dplyr::filter(transcript_id %in% representative_transcript &
#                                                                              type == "transcript")))
# per_representative_transcript_exon_gr <- gr_refactor_seqs(dt2gr(genes_of_interest %>%
#                                                                   dplyr::filter(transcript_id %in% representative_transcript &
#                                                                                   type == "exon")))
# per_representative_transcript_utr_gr <- gr_refactor_seqs(dt2gr(genes_of_interest %>%
#                                                                  dplyr::filter(transcript_id %in% representative_transcript &
#                                                                                  type %in% c("five_prime_utr","three_prime_utr"))))
#
# # Now build the CDS for each representative transcript per gene
# # subtract out 5' UTR and 3' UTR
# per_representative_transcript_cds <- data.table()
# for(i in 1:n_distinct(per_representative_transcript_utr_gr$gene_name)) {
#   per_representative_transcript_cds <- rrbind(per_representative_transcript_cds,
#                                               gr2dt(gr.setdiff(query = per_representative_transcript_exon_gr %Q% (gene_name == unique(per_representative_transcript_utr_gr$gene_name)[i]),
#                                                                subject = per_representative_transcript_utr_gr %Q% (gene_name == unique(per_representative_transcript_utr_gr$gene_name)[i]),
#                                                                ignore.strand = T)))
# }
#
# # Need to add in the set of genes that do not have any UTR to the CDS segments
# per_representative_transcript_cds_gr <- gr_refactor_seqs(dt2gr(rrbind(per_representative_transcript_cds,
#                                                                       gr2dt(per_representative_transcript_exon_gr %Q% (!gene_name %in% unique(per_representative_transcript_utr_gr$gene_name))))))
#
# # To find the intronic regions subtract the exonic regions from the transcript
# per_representative_transcript_intron <- data.table()
# for(i in 1:n_distinct(per_representative_transcript_exon_gr$gene_name)) {
#   per_representative_transcript_intron <- rrbind(per_representative_transcript_intron,
#                                                  gr2dt(gr.setdiff(query = per_representative_transcript_gr %Q% (gene_name == unique(per_representative_transcript_exon_gr$gene_name)[i]),
#                                                                   subject = per_representative_transcript_exon_gr %Q% (gene_name == unique(per_representative_transcript_exon_gr$gene_name)[i]),
#                                                                   ignore.strand = T)))
#   if(i %% 1000 == 0) {
#     message(paste0(i, " done ..."))
#   }
# }
#
# # Replace the type designation to intron then convert to GR obj
# per_representative_transcript_intron$type <- "intron"
# per_representative_transcript_intron_gr <- gr_refactor_seqs(dt2gr(per_representative_transcript_intron))
#
# # The final complete gene body set is the CDS + UTR + intron
# complete_gene_body_gr <- gr_refactor_seqs(grbind(per_representative_transcript_cds_gr,
#                                                  per_representative_transcript_intron_gr,
#                                                  per_representative_transcript_utr_gr))
# # Convert to DT to build the final GB DT
# complete_gene_body <- gr2dt(complete_gene_body_gr)
#
# # Final GDT / GGR
# # GRCh38 coordinates in GR obj format (seqnames, start, end, strand)
# # Genomic context of the gene body segment (type)
# # Prioritize gene annotations by gene name (gene_name,gene_id,gene_biotype,gene_source,transcript_id,transcript_name)
# # Some additional classifications for plotting (orientation,coding_v_noncoding,exonic_v_intronic)
# final_gene_dt <- complete_gene_body %>%
#   dplyr::select(seqnames,start,end,type,gene_name,gene_id,gene_biotype,gene_source,transcript_id,transcript_name,exon_number,tag) %>%
#   dplyr::left_join(y = gr2dt(per_representative_transcript_gr) %>% dplyr::select(transcript_id,strand),
#                    by = "transcript_id") %>%
#   dplyr::mutate("coding_v_noncoding" = case_when(type == "exon" ~ "coding",
#                                                  type != "exon" ~ "non_coding"),
#                 "exonic_v_intronic" = case_when(type %in% c("exon","five_prime_utr","three_prime_utr") ~ "exonic",
#                                                 type == "intron" ~ "intronic"),
#                 "orientation" = case_when(strand == "+" ~ 1,
#                                           strand == "-" ~ 0)) %>%
#   dplyr::rename("segment" = type) %>%
#   dplyr::select(seqnames,start,end,strand,segment,
#                 gene_name,gene_id,gene_biotype,
#                 gene_source,transcript_id,transcript_name,exon_number,tag,
#                 orientation,coding_v_noncoding,exonic_v_intronic)

# final_gene_gr <- gr_refactor_seqs(dt2gr(final_gene_dt))

# Don't rerun the code snippet but simply load the data here
gene_body_hg38 <- read_bed_file(bed_file_path = fs::path_package("extdata", "gene_body.Ensembl_v108.hg38.bed.gz", package = "devgru"),
                                has_header = T)
# Export for usage
usethis::use_data(gene_body_hg38, internal = FALSE, overwrite = TRUE)


#
#
# }}}}------->>> DESeq2 Differential Gene Expression Demo
#
#

# Raw RNA-seq data for the human myeloma cell lines (HMCLs) RPMI8226 and U266
# were obtained from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87585

# The data was processed using the following workflow
# https://igordot.github.io/sns/routes/rna-star.html
#   - Trim adapters and low quality bases (Trimmomatic).
#   - Align to the reference genome (STAR).
#   - Align to other species and common contaminants (fastq_screen).
#   - Generate normalized genome browser tracks.
#   - Determine the distribution of the bases within the transcripts and 5’/3’ biases (Picard).
#   - Determine if the library is stranded and the strand orientation.
#   - Generate genes-samples counts matrix (featureCounts)

hmcl_counts_demo_dt_hg38 <- data.table::fread(input = fs::path_package("extdata", "hmcl_counts.hg38.txt.gz", package = "devgru"))

hmcl_conditions_demo_dt_hg38 <- data.table::fread(input = fs::path_package("extdata", "hmcl_conditions.hg38.txt", package = "devgru"))

# Export for usage
usethis::use_data(hmcl_counts_demo_dt_hg38, internal = FALSE, overwrite = TRUE)
usethis::use_data(hmcl_conditions_demo_dt_hg38, internal = FALSE, overwrite = TRUE)




#
#
# }}}}------->>> My Color Pal
#
#

devgru_pal <- c("#183e70","#448bca","#b5b2d9","#e01183","#fdbe3a","#833894",
                "#149ba7","#eb1d25","#006748","#fce200","#191970","#EEA9B8")




