# PiGx BSseq Pipeline.
#
# Copyright 2019 Alexander Blume <alexander.blume@mdc-berlin.de>
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

# This file was originally taken and modified from
# https://github.com/katwre/makeWGBSnake/blob/master/Rules/Meth_preprocessing_methyldackel_rules.py

# WGBS pipeline
#
# Copyright © 2018 Katarzyna Wreczycka katarzyna.wreczycka@mdc-berlin.de
# This pipeline is heavily based on the PiGx BSseq pipeline github.com/BIMSBbioinfo/pigx_bsseq

# ==========================================================================================
# Extract methylation counts with methylDackel


def getContextArg(context):
    ## dependent on context we may return a combination of the following args
    ## ""                         Ouptut only CpG context (returned by default)
    ## "--noCpG --CHG"            Output only CHG context methylation metrics
    ## "--noCpG --CHH"            Output only CHH context methylation metrics
    contextRequest = str(context).upper()
    contextArg = ""
    if contextRequest != "CPG":
        contextArg = f"--noCpG --{contextRequest}"

    return contextArg


def protocol(sample):
    samples = getSamplesPerMergeRep(sample, config["SAMPLES"])
    protocol = [config["SAMPLES"][sample]['Protocol'].upper() for sample in samples]
    return str(set(protocol)).upper()


def keepDups(protocol):
    keepDupsFlag = str(config["general"]["methylation-calling"]["keep-duplicates"]).lower()
    keepDups = ""
    if keepDupsFlag == 'auto':
        if protocol == "RRBS":
            keepDups = "--keepDupes"
    elif keepDupsFlag in {'true', 'yes'}:
        keepDups = "--keepDupes"

    return keepDups


rule methyldackel_extract_methylKit_by_context:
    input:
        bamfile=DIR_sorted + "{sample}.bwameth.sorted.markdup.bam",
        genome=GENOMEFILE,
    output:
        contextCallFile=DIR_methcall + "methylDackel/" + "{sample}_methyldackel_{context}.methylKit.gz",
    wildcard_constraints:
        sample=".+(?<!deduped)",
    params:
        threads=config['execution']['rules']['methyldackel_extract']['threads'],
        prefix=DIR_methcall + "methylDackel/" + "{sample}_methyldackel",
        context=lambda wc: getContextArg(wc.context),
        protocol=lambda wc: protocol(wc.sample),
        keepDups=lambda wc: keepDups(protocol(wc.sample)),
        minqual=int(config['general']['methylation-calling']['minimum-quality']),
    log:
        DIR_methcall + "methylDackel/" + "{sample}.methyldackel_{context}_calls.log",
    message:
        fmt(
            "Calling methylation for {wildcards.context} context in sample {wildcards.sample} according to protocol {params.protocol}."
        )
    shell:
        nice(
            "methyldackel",
            [
                "extract",
                "{input.genome}",
                "{input.bamfile}",
                "-o {params.prefix}",
                "-@ {params.threads}",
                "{params.keepDups}",
                "--methylKit",
                "{params.context}",
                "-q {params.minqual}",
                ">> {log} 2>&1",
                ";\n",
                "gzip",
                "{params.prefix}_{wildcards.context}.methylKit",
            ],
            ("{log}"),
        )


rule methyldackel_extract_methylKit_deduped:
    input:
        bamfile=DIR_sorted + "{sample}.bwameth.sorted.markdup.bam",
        genome=GENOMEFILE,
    output:
        contextCallFile=DIR_methcall + "methylDackel/" + "{sample}.deduped_methyldackel_{context}.methylKit.gz",
    params:
        threads=config['execution']['rules']['methyldackel_extract']['threads'],
        prefix=DIR_methcall + "methylDackel/" + "{sample}.deduped_methyldackel",
        context=lambda wc: getContextArg(wc.context),
        protocol=lambda wc: protocol(wc.sample),
        keepDups=lambda wc: keepDups(protocol(wc.sample)),
        minqual=int(config['general']['methylation-calling']['minimum-quality']),
    log:
        DIR_methcall + "methylDackel/" + "{sample}.deduped.methyldackel_{context}_calls.log",
    message:
        fmt(
            "Calling methylation for {wildcards.context} context in sample {wildcards.sample} according to protocol {params.protocol}."
        )
    shell:
        nice(
            "methyldackel",
            [
                "extract",
                "{input.genome}",
                "{input.bamfile}",
                "-o {params.prefix}",
                "-@ {params.threads}",
                "{params.keepDups}",
                "--methylKit",
                "{params.context}",
                "-q {params.minqual}",
                ">> {log} 2>&1",
                ";\n",
                "gzip",
                "{params.prefix}_{wildcards.context}.methylKit",
            ],
            ("{log}"),
        )


# ==========================================================================================
# Extract methylation bias with methylDackel


def removeMapperSuffix(prefix):
    mapper_suffix = [".bwameth.sorted.markdup"]  # bwameth
    mapper_suffix += [
        "_1_val_1_bt2.sorted",
        "_se_bt2.sorted",
        "_1_val_1_bt2.sorted.deduped",
        "_se_bt2.sorted.deduped",
    ]  # bismark
    for suffix in mapper_suffix:
        if prefix.endswith(suffix):
            return prefix[: -len(suffix)]


rule methyldackel_mbias:
    input:
        bamfile=DIR_sorted + "{prefix}.bam",
        genome=ancient(GENOMEFILE),
    output:
        txt=DIR_methcall + "methylDackel/" + "{prefix}_methylDackel_mbias_{context}.txt",
        p1=DIR_methcall + "methylDackel/" + "{prefix}_methylDackel_mbias_{context}_OB.svg",
        p2=DIR_methcall + "methylDackel/" + "{prefix}_methylDackel_mbias_{context}_OT.svg",
    params:
        threads=config['execution']['rules']['methyldackel_extract']['threads'],
        prefix=DIR_methcall + "methylDackel/" + "{prefix}_methylDackel_mbias_{context}",
        sample=lambda wc: removeMapperSuffix(wc.prefix),
        context=lambda wc: getContextArg(wc.context),
        protocol=protocol("{params.sample}"),
        keepDups=keepDups("{params.protocol}"),
        minqual=int(config['general']['methylation-calling']['minimum-quality']),
    log:
        DIR_methcall + "methylDackel/" + "{prefix}_methylDackel_mbias_{context}.log",
    message:
        fmt(
            "Calculate methylation bias for context {wildcards.context} using MethylDackel for sample  with prefix {wildcards.prefix}"
        )
    shell:
        nice(
            "methyldackel",
            [
                "mbias",
                "{input.genome}",
                "{input.bamfile}",
                "{params.prefix}",
                "{params.keepDups}",
                "{params.context}",
                "-@ {params.threads}",
                "--txt",
                "-q {params.minqual}",
                "> {output.txt}",
                "2> {log}",
            ],
        )


# ==========================================================================================
# Convert to tabix files with methylKit
#


rule tabix_methyldackelfile:
    input:
        DIR_methcall + "methylDackel/" + "{prefix}_methyldackel_{context}.methylKit.gz",
    output:
        DIR_methcall + "methylDackel/" + "tabix_{context}/{prefix}_{context}.txt.bgz",
        DIR_methcall + "methylDackel/" + "tabix_{context}/{prefix}_{context}.txt.bgz.tbi",
    params:
        sampleid="{prefix}_{context}",
        assembly=ASSEMBLY,
        context="{context}",
        dbdir=DIR_methcall + "methylDackel/" + "/tabix_{context}/",
        mincov=int(config['general']['methylation-calling']['minimum-coverage']),
    log:
        DIR_methcall + "methylDackel/" + "tabix_{context}/{prefix}_{context}.makeTabix.log",
    message:
        fmt("Create Tabix file from MethylDackel file for sample {wildcards.prefix} and context {params.context}")
    shell:
        nice(
            'Rscript',
            [
                "{DIR_scripts}/makeTabix.R",
                "--location={input}",
                "--sample.id={params.sampleid}",
                "--assembly={params.assembly}",
                "--context={params.context}",
                "--mincov={params.mincov}",
                "--dbdir={params.dbdir}",
                "--logFile={log}",
            ],
            "{log}",
        )
