# awk script to convert ucsc table dump to
# bed12 derived from https://gist.github.com/cbp44/e6c198d27fa67dcf6a749ed474fcff73#file-ncbi_refseq_tsv_to_exon_bed12-awk

BEGIN {
    IFS="\t"; OFS="\t";
}
{
    n=split($10, exon_starts, ",");
    split($11, exon_ends, ",");
    txStart=$5;

    exon_sizes_str=exon_ends[1]-exon_starts[1];
    exon_starts_str=exon_starts[1]-txStart;
    for (i=2; i<n; i++) {
        exon_sizes_str=exon_sizes_str","exon_ends[i]-exon_starts[i];
        exon_starts_str=exon_starts_str","(exon_starts[i]-txStart);
    };
    exon_sizes_str=exon_sizes_str",";
    exon_starts_str=exon_starts_str",";

    print $3,$5,$6,$2,0,$4,$7,$8,"33,33,33",$9,exon_sizes_str,exon_starts_str
}
