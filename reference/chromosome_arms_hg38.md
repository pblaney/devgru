# Genomic coordinates of each chromosome arm on hg38 as GenomicRanges

A convenient single collection of the genomic span of each chromosome
arm. The span is described as: \| ranges kept \| p-arm start
\|\<————-\>\| centromere q-arm centromere \|\<————-\>\| end

## Format

`GRanges`

## Details

Note, this data does not include a range across the centromeres, for
these see `exclusion_regions_hg38`.
