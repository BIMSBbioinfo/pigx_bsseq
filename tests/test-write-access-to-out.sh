#!/bin/sh

set -e
set -u

test_name=$(basename $0)
dir=$(mktemp --tmpdir -d "pigx.XXXX")
testdir="${dir}/${test_name}"
mkdir -p ${dir}/${test_name}
cd ${dir}/${test_name}

echo "TEST: $test_name"
echo "DIRECTORY: $dir"

# create test data
mkdir -p in
mkdir -p genome
touch genome/sample.fasta

cat <<EOF >settings.yaml
locations:
  input-dir: ${testdir}/in/
  output-dir: ${testdir}/out/
  genome-fasta: ${testdir}/genome/sample.fasta
general:
  assembly: hg19
EOF
cat <<EOF >samplesheet.csv
Read1,Read2,SampleID,Protocol,Treatment
EOF

# run the pipeline
export PIGX_UGLY=1
export PIGX_UNINSTALLED=1

# This should fail, because the genome directory does not contain
# bisulfite genome files, and we have no write access to the genome
# directory.
chmod -w genome

${builddir}/pigx-bsseq --target=genome-prep -s settings.yaml samplesheet.csv 2>&1 |
  grep "ERROR: reference genome has not been bisulfite-converted, and PiGx does not have permission to write to that directory."

chmod +w genome
rm -rf ${dir}
