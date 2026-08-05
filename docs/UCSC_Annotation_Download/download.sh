# script to download and convert ucsc table dump to bed12 
# awk script derived from https://gist.github.com/cbp44/e6c198d27fa67dcf6a749ed474fcff73#file-ncbi_refseq_tsv_to_exon_bed12-awk

REF=${1:-hg38}


echo "downloading refGene"
wget -qO- http://hgdownload.cse.ucsc.edu/goldenpath/${REF}/database/refGene.txt.gz \
   | gunzip -c \
   | awk -f ucsc_refseq_tsv_to_exon_bed12.awk - \
   | gzip > refGene.${REF}.bed.gz

echo "downloading ncbiRefSeq"
wget -qO- http://hgdownload.cse.ucsc.edu/goldenpath/${REF}/database/ncbiRefSeq.txt.gz \
   | gunzip -c \
   | awk -f ucsc_refseq_tsv_to_exon_bed12.awk - \
   | gzip > ncbiRefSeq.${REF}.bed.gz

echo "downloading cpgIslandExt"
wget -qO- http://hgdownload.cse.ucsc.edu/goldenpath/${REF}/database/cpgIslandExt.txt.gz \
   | gunzip -c \
   | awk 'BEGIN{ IFS="\t"; OFS="\t"; }{ print substr($0, index($0, $2)); }' \
   | sort -k 1,1 -k2,2n - \
   | gzip > cpgIslandExt.${REF}.bed.gz

