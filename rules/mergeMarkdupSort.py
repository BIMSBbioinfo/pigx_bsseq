# PiGx BSseq Pipeline.
#
# Copyright 2022 Alexander Blume <alexander.blume@mdc-berlin.de>
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
# Merging replicates, marking duplicates and sorting bam files in one go without writing
# intermediate files

# Example:
#   samtools merge -o - SEsample.bwameth.bam SEsample_v2.bwameth.bam SEsample_v2copy.bwameth.bam |
#   samtools fixmate -m - - |
#   samtools sort - |
#   samtools markdup -f SEsample.merged.samtools_markdup.sorted.bwameth.samtools_markdup_metrics.txt - -  |
#   samtools sort -o SEsample.merged.samtools_markdup.sorted.bwameth.bam


def runFixmate(sample, log="", samples_dict=config["SAMPLES"]):
    samples = getSamplesPerMergeRep(mergeRep=sample, samples_dict=samples_dict)
    if not samples:
        samples = [sample]
    file_nums = [len(samples_dict[sample]["files"]) for sample in samples]
    cmd = ""
    if list(set(file_nums)) != [1]:
        ## FIXME: the logging is quite hacky, maybe sometime we can pass it via lambda
        log = DIR_sorted + f"{sample}.bwameth.sorted.markdup.log"
        cmd = f" ".join(
            [
                tool("samtools"),
                "collate -u -O",
                "2>> {log}",
                "|",
                tool("samtools"),
                "fixmate -m -u - - ",
                f"2>> {log}",
                "|",
            ]
        )
    return cmd


def mergeMarkdupSort_input(sample, samples_dict=config["SAMPLES"]):
    samples = getSamplesPerMergeRep(mergeRep=sample, samples_dict=samples_dict)
    if not samples:
        samples = [sample]
    files = [DIR_mapped + f"{sample}.bwameth.bam" for sample in samples]
    return files


rule mergeMarkdupSort_bam:
    input:
        lambda wc: mergeMarkdupSort_input(wc.sample),
    output:
        bam=DIR_sorted + "{sample}.bwameth.sorted.markdup.bam",
        bai=DIR_sorted + "{sample}.bwameth.sorted.markdup.bam.bai",
        metric=DIR_sorted + "{sample}.bwameth.sorted.markdup.metrics.txt",
    params:
        runFixmate=lambda wc: runFixmate(sample=wc.sample),
        numReps=lambda wc: len(mergeMarkdupSort_input(wc.sample)),
    resources:
        threads=config["execution"]["rules"]["mergeMarkdupSort_bam"]["threads"],
        memory=config["execution"]["rules"]["mergeMarkdupSort_bam"]["memory"],
    log:
        DIR_sorted + "{sample}.bwameth.sorted.markdup.log",
    message:
        fmt(
            "Merging {params.numReps} replicates, marking duplicates and sorting bam files for {wildcards.sample}"
        )
    shell:
        nice(
            "samtools",
            [
                "merge -o -u - {input}",
                "2>> {log}",
                " | ",
                "{params.runFixmate}",
                tool("samtools"),
                "sort -u",
                "-m {resources.memory}",
                "-@ {resources.threads}",
                "-",
                "2>> {log}",
                "|",
                tool("samtools"),
                "markdup",
                "-f {output.metric}",
                "-@ {resources.threads}",
                "- {output.bam}",
                "2>> {log};",
                tool("samtools"),
                "index {output.bam}",
                "2>> {log}",
            ],
            "{log}",
        )
