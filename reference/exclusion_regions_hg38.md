# Collection of genomic regions to exclude from analysis on hg38 as GenomicRanges

A convenient single collection of genomic regions suggested to exclude
from analysis due to mappability issues, enrichment of duplicated
regions, telomeric, centromeric, acrocentric short-arms, and highly
heterochromatin/gene desert portions of the chromosomes.

## Format

`GRanges`

## Details

Built from UCSC and ENCODE blacklists. See excluderanges
(https://github.com/dozmorovlab/excluderanges).
