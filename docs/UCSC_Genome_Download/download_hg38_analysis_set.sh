#!/usr/bin/env bash
# Downloads the UCSC hg38 no-ALT analysis set and, if a spike-in FASTA is
# given, builds the bwa-meth hybrid reference (genome + spike-in) and
# indexes it with bwameth.py.
#
# Use this instead of a plain hg38.fa.gz download: the plain build carries
# ALT/patch/random contigs with no bwa .alt companion file, which causes
# MAPQ0 coverage dropout when a read is equally consistent with a primary
# chromosome and its own ALT-haplotype contig. See README.md in this
# directory.
#
# Usage: ./download_hg38_analysis_set.sh <OUT_DIR> [SPIKEIN_FASTA]
#   OUT_DIR        directory to download/build the reference in
#   SPIKEIN_FASTA  optional; if given, appended to the analysis set to
#                  build the hybrid reference, which is then indexed with
#                  `bwameth.py index`

set -euo pipefail

OUT_DIR="${1:?Usage: $0 <OUT_DIR> [SPIKEIN_FASTA]}"
SPIKEIN_FASTA="${2:-}"

URL="http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz"

mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

echo "downloading hg38 no-ALT analysis set"
wget -q "${URL}"
gunzip -f hg38.analysisSet.fa.gz

if [ -n "${SPIKEIN_FASTA}" ]; then
    HYBRID_FASTA="hg38_analysisSet_spikein.fa"
    echo "building hybrid reference with ${SPIKEIN_FASTA} -> ${HYBRID_FASTA}"
    cat hg38.analysisSet.fa "${SPIKEIN_FASTA}" > "${HYBRID_FASTA}"

    echo "indexing hybrid reference (bwameth.py index)"
    bwameth.py index "${HYBRID_FASTA}"

    echo "done. reference: ${OUT_DIR}/${HYBRID_FASTA}"
else
    echo "done. reference: ${OUT_DIR}/hg38.analysisSet.fa (not indexed - no spike-in given)"
fi
