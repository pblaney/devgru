# Gene body segments from Ensemble v108 on hg38 as GenomicRanges

This genomic ranges contains the gene body segments (UTRs, Exons,
Introns) for a single-transcript representation of Ensembl v108 genes on
hg38.

## Format

`GRanges`

## Details

                         |  gene body schematic  |

\[5` prime UTR]{Exon}<---Intron--->{Exon}<---Intron--->{Exon}[3` prime
UTR\]

The genes included had the following `gene_biotype`: protein_coding
IG_C_gene,IG_D_gene,IG_J_gene,IG_V_gene,
TR_C_gene,TR_D_gene,TR_J_gene,TR_V_gene lncRNA,miRNA

The transcript for protein-coding genes were determined using MANE
Select & MANE Plus Clinical v1.4 annotations

Built from

- Ensembl v108 GTF
  (https://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.gtf.gz)

- MANE v1.4 summary file
  (https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.4/MANE.GRCh38.v1.4.summary.txt.gz)
