# PiGx BSseq Pipeline.
#
# Copyright 2026 Alexander Blume <alexander.blume@mdc-berlin.de>
#
# This file is part of the PiGx BSseq Pipeline.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# ==========================================================================================
# GPU-accelerated bisulfite alignment using NVIDIA Parabricks fq2bam_meth.
#
# Replaces bwameth_align_trimmed + samblaster_markdup_sort in a single step.
# Output path is identical to the CPU bwameth track so all downstream rules
# (methylDackel, bigwig, diffmeth) are shared unchanged.
#
# Requires: Apptainer and a CUDA-capable GPU (compute capability >= 75, >= 16 GB VRAM).
# Container image is configured under tools.fq2bam_meth in settings.yaml.
# ==========================================================================================


rule bwameth_genome_preparation:
    input:
        ancient(GENOMEFILE)
    output:
        GENOMEFILE+".bwameth.c2t.sa",
        GENOMEFILE+".bwameth.c2t.amb",
        GENOMEFILE+".bwameth.c2t.ann",
        GENOMEFILE+".bwameth.c2t.pac",
        GENOMEFILE+".bwameth.c2t.bwt",
        GENOMEFILE+".bwameth.c2t"
    log:
        os.path.join(GENOMEPATH,'bwameth_genome_preparation.output_'+ASSEMBLY+'.log')
    message: "Converting {ASSEMBLY} Genome into Bisulfite analogue with bwa-meth"
    shell:
        nice("bwameth", ["index {input}"],"{log}")

rule bwameth_touch_index:
    input:
        GENOMEFILE+".bwameth.c2t.sa",
        GENOMEFILE+".bwameth.c2t.amb",
        GENOMEFILE+".bwameth.c2t.ann",
        GENOMEFILE+".bwameth.c2t.pac",
        GENOMEFILE+".bwameth.c2t.bwt"
    output:
        GENOMEFILE+".bwameth.c2t_was.touched"
    message: "Update timestamp for {ASSEMBLY} Genome Index"
    shell:
        "sleep 60; touch {input};touch {output}"


def bwameth_input(sample):
    files = list_files_TG(samplesheet(sample, 'files'), sample, '')
    return files

def fq2bam_meth_input_param(sample):
    files = list_files_TG(samplesheet(sample, 'files'), sample, '')
    if len(files) == 1:
        return f"--in-se-fq {' '.join(files)}"
    else:
        return f"--in-fq {' '.join(files)}"

rule fq2bam_meth_align:
    input:
        rules.bwameth_touch_index.output,
        index = rules.bwameth_genome_preparation.output,
        files = lambda wc: bwameth_input(wc.sample)
    output:
        bam   = DIR_sorted + "{sample}.bwameth.sorted.markdup.bam",
        index = DIR_sorted + "{sample}.bwameth.sorted.markdup.bam.bai"
    params:
        container = tool('fq2bam_meth'),
        outdir    = DIR_sorted,
        input_param = lambda wc: fq2bam_meth_input_param(wc.sample)
    resources:
        nvidia_gpu = config['execution']['rules']['fq2bam_meth_align']['gpus']
    threads:
        config['execution']['rules']['fq2bam_meth_align']['threads']
    log:
        DIR_sorted + "{sample}_fq2bam_meth.log"
    message: fmt("Aligning bisulfite reads with GPU fq2bam_meth for sample {wildcards.sample}.")
    shell:
        """
        apptainer run --nv \
          -B {GENOMEPATH} \
          -B {OUTDIR} \
          --pwd {params.outdir} \
          {params.container} \
          pbrun fq2bam_meth \
          --ref {GENOMEFILE} \
          {params.input_param} \
          --out-bam {output.bam} \
          --low-memory \
          --num-gpus {resources.nvidia_gpu} \
          > {log} 2>&1
        """
