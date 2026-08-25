## Download Genome Reference

Use the UCSC hg38 **no-ALT analysis set** as the genome reference for this
pipeline — not a plain `hg38.fa.gz` download.

Why: a plain hg38 build carries ALT/patch/random contigs with no bwa `.alt`
companion file. Reads that are equally consistent with a primary chromosome
and its own ALT-haplotype contig (e.g. `chr21_KI270874v1_alt`) then get
MAPQ0 and are dropped by the default MethylDackel filter, causing large
false gaps in methylation coverage (one observed case: a 74kb dropout on
chr21). The analysis set is purpose-built for NGS alignment pipelines and,
in the no-ALT variant used here, removes ALT contigs entirely so this
ambiguity can't arise. It also hard-masks centromeric-array/WGS-duplicate
regions and the chrY PAR regions, and adds an EBV decoy sequence to absorb
contamination. See
[UCSC's analysis set docs](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/)
for details.

You may use the download script to fetch the reference and (optionally)
build + index the bisulfite hybrid reference (genome + methylation
spike-in) used by this pipeline:

```sh
mkdir -p <path_to_genomes>/hg38
bash download_hg38_analysis_set.sh <path_to_genomes>/hg38 <path_to_spikein.fa>
```

Omit the spike-in FASTA argument to just download the plain analysis set
without building/indexing a hybrid reference:

```sh
bash download_hg38_analysis_set.sh <path_to_genomes>/hg38
```
