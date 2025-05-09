# }}}}------->>> Code used to generate the internal data used for demos.

library(devgru)
library(gUtils)
library(dplyr)
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
read_depth_demo_hg38 <- read_bed_file(bed_file_path = fs::path_package("extdata", "read_depth_profile.hg38.bed", package = "devgru"),
                                      has_header = T)
# Export for usage
usethis::use_data(read_depth_demo_hg38, internal = FALSE, overwrite = TRUE)
