# PiGx BSseq Pipeline.
#
# Copyright © 2017 Bren Osberg <brendan.osberg@mdc-berlin.de>
# Copyright © 2017, 2018 Alexander Gosdschan <alexander.gosdschan@mdc-berlin.de>
# Copyright © 2017 Katarzyna Wreczycka <katwre@gmail.com>
# Copyright © 2017, 2018 Ricardo Wurmus <ricardo.wurmus@mdc-berlin.de>
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

import os
import sys
from datetime import datetime
from glob import glob

# ==============================================================================
#
#                                       HELPER FUNCTIONS
#
# put here helper functions that are used within more than one rule   
# ==============================================================================


# --------------------------------------
# general purpose
# --------------------------------------

def fmt(message):
    """Format the MESSAGE string."""
    return "----------  " + message + "  ----------"

def tool(name):
    return config['tools'][name]['executable']

def toolArgs(name):
    if 'args' in config['tools'][name]:
        return config['tools'][name]['args']
    else:
        return ""

# sample sheet accessor functions
def samplesheet(name, item=None):
    """Access the SAMPLES dict from config."""
    if item:
        return config["SAMPLES"][name][item]
    else:
        config["SAMPLES"][name]

# print datetime next to a message
def log_time(text):
    time = '$(date +"[%Y-%m-%d %T]")'
    log = '{} {}'.format(time, text)
    return(fmt(log))

# Generate a command line string that can be passed to snakemake's
# "shell".  The string is prefixed with an invocation of "nice".
def nice(cmd, args, log=None, fallback=None):
    executable = tool(cmd)
    line = ["nice -" + str(config['execution']['nice']),
            executable] + [toolArgs(cmd)] + args
    if log:
        line.insert(0, "echo {} > {};".format(log_time("Starting Now"),log))
        line.append(">> {} 2>&1".format(log))
    if fallback:
        line.append(" || {} ".format(fallback))
    if log:
        line.append("; echo {} >> {};".format(log_time("Done"),log)) 
    return " ".join(line)


# abandone current execution with a helpful error message:
def bail(msg):
    """Print the error message to stderr and exit."""
    print(msg, file=sys.stderr)
    exit(1)

def TrueOrFalse(value):
    value = repr(value)
    if value:
        answer = value.lower() in ["true","yes"]
    else:
        answer = False
    return(answer)



# --------------------------------------
# sample related      
# --------------------------------------
        
def dedupe_tag(protocol):
    if protocol.upper() == "WGBS":
        return ".deduped"
    elif protocol.upper() == "RRBS":
        return ""
    else:
        raise Exception("=== ERROR: unexpected protocol ===")


def get_fastq_name(full_name):
    # single end
    find_se_inx = full_name.find('_se_bt2')
    # paired-end
    find_pe_inx = full_name.find('_1_val_1_bt2')

    if(find_se_inx >= 0):
        output = full_name[:find_se_inx]
    elif(find_pe_inx >= 0):
        output = full_name[:find_pe_inx]
    else:
        bail("Unable to infer sample fastq name; cannot find trimming string in filename. \nHave the files been trimmed for adapter sequence and quality?")

    return(output)


# --------------------------------------
# context related      
# --------------------------------------
def destrand(context):
    return config["general"]["export-bigwig"]["context"][context.lower()]["destrand"]


def formatContext(context):
    contexts = {'cpg': 'CpG', 'chg': 'CHG', 'chh': 'CHH'}
    return contexts.get(context)


# --------------------------------------
# merging samples
# --------------------------------------
def getMergeRepPerSample(sample, samples_dict):
    """
    Extract the values of 'MergeReplicates' from dictionary 'samples_dict'
    for a given 'sample' (SampleID).

    Returns given 'sample' if column 'MergeReplicates' is empty or not found
    """
    mergeReps = (
        samples_dict[sample]["MergeReplicates"]
        if ("MergeReplicates" in samples_dict[sample])
        and (samples_dict[sample]["MergeReplicates"])
        else sample
    )
    # samples_dict[sample]["MergeReplicates"]
    # if not mergeReps: mergeReps = sample
    return mergeReps


def getSamplesPerMergeRep(mergeRep, samples_dict):
    """
    Extract the values of 'SampleID' from dictionary 'samples_dict'
    for a given 'mergeRep' (MergeReplicates).

    Returns given 'mergeRep' if column 'MergeReplicates' is empty or not found
    """
    samples = [
        sample
        for sample in samples_dict
        if ("MergeReplicates" in samples_dict[sample])
        and (samples_dict[sample]["MergeReplicates"] == mergeRep)
    ]
    return samples


# --------------------------------------
# generate dynamic output files per rule     
# --------------------------------------

def files_for_sample(proc):
    return [expand(proc(config['SAMPLES'][sample]['files'],
                        config['SAMPLES'][sample]['SampleID'],
                        config['SAMPLES'][sample]['Protocol']))
            for sample in config['SAMPLES']]


def list_files_rawQC(files, sampleID, protocol):
    PATH = DIR_rawqc
    infix = ["_fastqc"] if len(files) == 1 else ["_1_fastqc", "_2_fastqc"]
    return [os.path.join(PATH, f"{sampleID}{inf}.zip") for inf in infix]


def list_files_TG(files, sampleID, protocol):
    PATH = DIR_trimmed
    infix = ["_trimmed"] if len(files) == 1 else ["1_val_1", "_2_val_2"]
    return [os.path.join(PATH, f"{sampleID}{inf}.fq.gz") for inf in infix]


def list_files_posttrim_QC(files, sampleID, protocol):
    PATH = DIR_posttrim_QC
    infix = ["_trimmed"] if len(files) == 1 else ["1_val_1", "_2_val_2"]
    return [os.path.join(PATH, f"{sampleID}{inf}_fastqc.zip") for inf in infix]


def list_files_bismark(files, sampleID, protocol):
    PATH = DIR_mapped
    infix = "_trimmed" if len(files) == 1 else "_1_val_1"
    suffix = (
        ["_SE_report.txt", ".bam"] if len(files) == 1 else ["_PE_report.txt", "_pe.bam"]
    )
    return [os.path.join(PATH, f"{sampleID}{infix}_bismark_bt2{suf}") for suf in suffix]


def list_files_bwameth(files, sampleID, protocol):
    PATH = DIR_mapped
    return [os.path.join(PATH, f"{sampleID}.bwameth.bam")]


def list_files_dedupe(files, sampleID, protocol):
    PATH = DIR_sorted
    infix = "_se" if len(files) == 1 else "_1_val_1"
    return [os.path.join(PATH, f"{sampleID}{infix}.sorted{dedupe_tag(protocol)}.bam")]


def list_files_markdup(files, sampleID, protocol):
    PATH = DIR_sorted
    sampleID = getMergeRepPerSample(sample=sampleID, samples_dict=config["SAMPLES"])
    return [
        PATH + sampleID + ".bwameth.sorted.markdup.bam",
        PATH + sampleID + ".bwameth.sorted.picard_MarkDuplicates.metrics.txt",
    ]


def list_files_bwamethMappingStats(files, sampleID, protocol):
    PATH = DIR_sorted
    sampleID = getMergeRepPerSample(sample=sampleID, samples_dict=config["SAMPLES"])
    return [
        PATH + sampleID + ".bwameth.sorted.markdup.idxstats.txt",
        PATH + sampleID + ".bwameth.sorted.markdup.stats.txt",
        PATH + sampleID + ".bwameth.sorted.markdup.flagstat.txt",
    ]


def bam_processing(files, sampleID, protocol):
    PATH = DIR_methcall + "methylDackel/"
    # ---- change based on single or paired end
    mapper_suffix = "_se_bt2.sorted" if len(files) == 1 else "_1_val_1_bt2.sorted"
    prefix = PATH + sampleID + mapper_suffix + dedupe_tag(protocol)
    return [
        # contexts are cpg, chg, chh
        f"{prefix}_{context}.txt.bgz{ext}"
        for ext in ["", ".tbi"]
        for context in METH_CONTEXTS
    ]


def list_files_methyldackel_extract(files, sampleID, protocol):
    PATH = DIR_methcall + "methylDackel/"
    sampleID = getMergeRepPerSample(sample=sampleID, samples_dict=config["SAMPLES"])
    prefix = PATH + sampleID + dedupe_tag(protocol)
    return [
        # contexts are CpG, CHG, CHH
        f"{prefix}_methyldackel_{formatContext(context)}.methylKit.gz"
        for context in METH_CONTEXTS
    ]


def list_files_methyldackel_mbias_bismark(files, sampleID, protocol):
    PATH = DIR_methcall + "methylDackel/"
    # ---- change based on single or paired end
    mapper_suffix = "_se_bt2.sorted" if len(files) == 1 else "_1_val_1_bt2.sorted"
    prefix = PATH + sampleID + mapper_suffix + dedupe_tag(protocol)
    return [
        f"{prefix}_methylDackel_mbias_{context}{ext}"
        for ext in [".txt", "_OB.svg", "_OT.svg"]
        for context in METH_CONTEXTS
    ]


def list_files_methyldackel_mbias_bwameth(files, sampleID, protocol):
    PATH = DIR_methcall + "methylDackel/"
    sampleID = getMergeRepPerSample(sample=sampleID, samples_dict=config["SAMPLES"])
    prefix = PATH + sampleID + ".bwameth.sorted.markdup"
    return [
        # NOTE: this should be rewritten with snakemakes expand()
            # NOTE: this should be rewritten with snakemakes expand()  
        # NOTE: this should be rewritten with snakemakes expand()
        f"{prefix}_methylDackel_mbias_{context}{ext}"
        for ext in [".txt", "_OB.svg", "_OT.svg"]
        for context in METH_CONTEXTS
    ]


# FIXME: contexts should be generate output based on settings file
def list_files_maketabix_methyldackel(files, sampleID, protocol):
    PATH = DIR_methcall + "methylDackel/"
    # ---- change based on single or paired end
    prefix = sampleID + dedupe_tag(protocol)
    return [
        # contexts are CpG, CHG, CHH
        os.path.join(
            PATH, f"tabix_{formatContext(context)}", f"{prefix}_{context}.txt.bgz{ext}"
        )
        for ext in ["", ".tbi"]
        for context in METH_CONTEXTS
    ]


def bigwig_exporting_bismark(files, sampleID, protocol):
    PATH = DIR_bigwig
    DESTRAND = "_destranded" if destrand("cpg") else ""
    # ---- change based on single or paired end
    mapper_suffix = "_se_bt2.sorted" if len(files) == 1 else "_1_val_1_bt2.sorted"
    prefix = PATH + sampleID + mapper_suffix + dedupe_tag(protocol)
    return [f"{prefix}.{context}{DESTRAND}_methylKit.bw" for context in METH_CONTEXTS]


# FIXME: contexts should be generate output based on settings file
def bigwig_exporting_bwameth(files, sampleID, protocol):
    PATH = DIR_bigwig
    sampleID = getMergeRepPerSample(sample=sampleID, samples_dict=config["SAMPLES"])
    DESTRAND = "_destranded" if destrand("cpg") else ""
    prefix = PATH + sampleID + dedupe_tag(protocol)
    return [
        f"{prefix}.{formatContext(context)}{DESTRAND}_methylDackel.bw"
        for context in METH_CONTEXTS
    ]


def methSeg_bismark(files, sampleID, protocol):
    PATH = DIR_seg
    # ---- change based on single or paired end
    mapper_suffix = "_se_bt2.sorted" if len(files) == 1 else "_1_val_1_bt2.sorted"
    prefix = PATH + sampleID + mapper_suffix + dedupe_tag(protocol)
    return [
        f"{prefix}_{context}_methylKit.meth_segments.bed" for context in METH_CONTEXTS
    ]


def methSeg_bwameth(files, sampleID, protocol):
    PATH = DIR_seg
    sampleID = getMergeRepPerSample(sample=sampleID, samples_dict=config["SAMPLES"])
    prefix = PATH + sampleID + dedupe_tag(protocol)
    return [
        f"{prefix}_{formatContext(context)}_methylDackel.meth_segments.bed"
        for context in METH_CONTEXTS
    ]


def list_final_reports_bwameth(files, sampleID, protocol):
    SUFFIX = "CpG_methylDackel_{}".format(ASSEMBLY)
    PATH = os.path.join(DIR_final, "sample_reports/")
    sampleID = getMergeRepPerSample(sample=sampleID, samples_dict=config["SAMPLES"])
    prefix = PATH + sampleID + dedupe_tag(protocol)
    return [
        f"{prefix}_{formatContext(context)}_methylDackel_{ASSEMBLY}_final.html"
        for context in METH_CONTEXTS
    ]


def list_final_reports_bismark(files, sampleID, protocol):
    PATH = os.path.join(DIR_final, "sample_reports/")
    mapper_suffix = "_se_bt2.sorted" if len(files) == 1 else "_1_val_1_bt2.sorted"
    prefix = PATH + sampleID + mapper_suffix + dedupe_tag(protocol)
    return [
        f"{prefix}_{context}_methylKit_{ASSEMBLY}_final.html"
        for context in METH_CONTEXTS
    ]


# --------------------------------------
# diffmeth related
# --------------------------------------
def get_sampleids_from_treatment(treatment):
    """Get SampleIDs from treatment string."""
    sample_ids = list(config["SAMPLES"].keys())
    sample_treatments = [samplesheet(s, "Treatment") for s in sample_ids]
    sampleids_list = [
        sample_ids[i] for i, x in enumerate(sample_treatments) if x == treatment
    ]
    return sampleids_list


def get_sampleids_from_analysis(analysis):
    """Get SampleIDs for each Analysis group."""
    sampleids_list = []
    for group in config["DManalyses"][analysis]:
        for treatment in config["DManalyses"][analysis][group].split(","):
            sampleids_list += get_sampleids_from_treatment(treatment)

    return sampleids_list


def makeDiffMethPath(path, suffix, treatment):
    return path + str(treatment).replace('vs', '_') + dedupe_tag(config["SAMPLES"][get_sampleids_from_treatment(treatment[0])[0]]['Protocol']) + '_' + suffix

def diffmeth_input_function(treatments):
    treatments = treatments.replace(".deduped", "")
    sampleids = get_sampleids_from_treatment(treatments)

    inputfiles = []
    for sampleid in sampleids:
        files = config["SAMPLES"][sampleid]["fastq_name"]
        protocol = config["SAMPLES"][sampleid]["Protocol"]
        infix = "_se" if len(files) == 1 else "_1_val_1"
        inputfile = [
            os.path.join(
                DIR_methcall,
                f"{sampleid}{infix}_bt2.sorted{dedupe_tag(protocol)}_methylRaw.RDS",
            )
        ]
        inputfiles.append(inputfile)

    inputfiles = list(sum(inputfiles, []))

    return inputfiles


def files_for_treatment(proc):
    if "DManalyses" in config.keys():
        treatment_groups = config["DManalyses"]
        files = []
        if treatment_groups:
            files = [
                expand(proc(comparison))
                for comparison in treatment_groups
                if comparison
            ]
        return files


# FIXME: create files dependend on context
def list_files_unite_bismark(treatment):
    PATH = DIR_diffmeth + treatment + "/"  
    DESTRAND = "_destranded" if destrand("cpg") else ""
    return [ PATH + "methylBase_"+treatment+"_cpg"+DESTRAND+"_methylKit.txt.bgz" ] 
     
def list_files_unite_bwameth(treatment):
    PATH = DIR_diffmeth + treatment + "/"  
    DESTRAND = "_destranded" if destrand("cpg") else ""
    return [ PATH + "methylBase_"+treatment+"_CpG"+DESTRAND+"_methylDackel.txt.bgz" ] 


def list_files_diffmeth_bismark(treatment):
    PATH = DIR_diffmeth + treatment + "/"  
    DESTRAND = "_destranded" if destrand("cpg") else ""
    return [ 
            PATH + "methylDiff_"+treatment+"_cpg"+DESTRAND+"_methylKit_full.txt.bgz",
            PATH + "methylDiff_"+treatment+"_cpg"+DESTRAND+"_methylKit_results.tsv"
            ]

    
def list_files_diffmeth_bwameth(treatment):
    PATH = DIR_diffmeth + treatment + "/"  
    DESTRAND = "_destranded" if destrand("cpg") else ""
    return [ 
            PATH + "methylDiff_"+treatment+"_CpG"+DESTRAND+"_methylDackel_full.txt.bgz",
            PATH + "methylDiff_"+treatment+"_CpG"+DESTRAND+"_methylDackel_results.tsv"
            ]

def list_files_diffmeth_report_bwameth(treatment):
    PATH = DIR_final 
    DESTRAND = "_destranded" if destrand("cpg") else ""
    return [ PATH + treatment + "/" + treatment + "_CpG"+DESTRAND+"_methylDackel"+ ".diffmeth-report.html"]

def list_files_diffmeth_report_bismark(treatment):
    PATH = DIR_final 
    DESTRAND = "_destranded" if destrand("cpg") else ""
    return [ PATH + treatment + "/" + treatment + "_cpg"+DESTRAND+"_methylKit" +".diffmeth-report.html"]
