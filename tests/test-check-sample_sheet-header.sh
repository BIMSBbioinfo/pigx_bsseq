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
EOF
cat <<EOF >sample_sheet.csv
foo,bar,baz
EOF

# run the pipeline
export PIGX_UGLY=1
export PIGX_UNINSTALLED=1

${builddir}/pigx-bsseq -s settings.yaml sample_sheet.csv 2>&1 |
  grep "First columns of the input table have to be"
rm -rf ${dir}
