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


# print a warning
def warn(msg):
    """Print the warning message to stderr."""
    print(f"WARNING: {msg}", file=sys.stderr)

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
    return config['general']['export-bigwig']['context'][context.lower()]['destrand']


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
    if len(files) == 1:
        return [PATH+sampleID+"_fastqc.zip"]  # ---- single end
    elif len(files) == 2:
        # ---- paired end
        return [PATH+sampleID+"_1_fastqc.zip", PATH+sampleID+"_2_fastqc.zip"]
    else:
        raise Exception(
            "=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")


def list_files_TG(files, sampleID, protocol):
    PATH = DIR_trimmed
    if len(files) == 1:
        return [PATH+sampleID+"_trimmed.fq.gz"]  # ---- single end
    elif len(files) == 2:
        # ---- paired end
        return [PATH+sampleID+"_1_val_1.fq.gz", PATH+sampleID+"_2_val_2.fq.gz"]
    else:
        raise Exception(
            "=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")


def list_files_posttrim_QC(files, sampleID, protocol):
    PATH = DIR_posttrim_QC
    if len(files) == 1:
        return [PATH+sampleID+"_trimmed_fastqc.zip"]  # ---- single end
    elif len(files) == 2:
        # ---- paired end
        return [PATH+sampleID+"_1_val_1_fastqc.zip", PATH+sampleID+"_2_val_2_fastqc.zip"]
    else:
        raise Exception(
            "=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")


def list_files_bismark(files, sampleID, protocol):
    PATH = DIR_mapped
    if len(files) == 1:
        return [PATH+sampleID+"_trimmed_bismark_bt2_SE_report.txt",
                PATH+sampleID+"_trimmed_bismark_bt2.bam"]  # ---- single end
    elif len(files) == 2:
        return [PATH+sampleID+"_1_val_1_bismark_bt2_PE_report.txt",
                PATH+sampleID+"_1_val_1_bismark_bt2_pe.bam"]  # ---- paired end
    else:
        raise Exception(
            "=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")


def list_files_bwameth(files, sampleID, protocol):
    PATH = DIR_mapped
    if len(files) == 1:
        return [PATH+sampleID+".bwameth.bam"]  # ---- single end
    elif len(files) == 2:
        return [PATH+sampleID+".bwameth.bam"]  # ---- paired end
    else:
        raise Exception(
            "=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")




def list_files_dedupe(files, sampleID, protocol):
    PATH = DIR_sorted
    if len(files) == 1:
        # ---- single end
        return [PATH+sampleID+"_se_bt2.sorted" + dedupe_tag(protocol) + ".bam"]
    elif len(files) == 2:
        # ---- paired end
        return [PATH+sampleID+"_1_val_1_bt2.sorted" + dedupe_tag(protocol) + ".bam"]
    else:
        raise Exception(
            "=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")


def list_files_markdup(files, sampleID, protocol):
    PATH = DIR_sorted
    if len(files) == 1:
        return [PATH+sampleID+".bwameth.sorted.markdup.bam"]  # ---- single end
    elif len(files) == 2:
        return [PATH+sampleID+".bwameth.sorted.markdup.bam"]  # ---- paired end
    else:
        raise Exception(
            "=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")


def list_files_bwamethMappingStats(files, sampleID, protocol):
    PATH = DIR_sorted
    if len(files) == 1:
        return [PATH+sampleID+".bwameth.sorted.markdup.idxstats.txt",
                PATH+sampleID+".bwameth.sorted.markdup.stats.txt",
                PATH+sampleID+".bwameth.sorted.markdup.flagstat.txt"]  # ---- single end
    elif len(files) == 2:
        return [PATH+sampleID+".bwameth.sorted.markdup.idxstats.txt",
                PATH+sampleID+".bwameth.sorted.markdup.stats.txt",
                PATH+sampleID+".bwameth.sorted.markdup.flagstat.txt"]  # ---- paired end
    else:
        raise Exception(
            "=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")



# FIXME: contexts should be generate output based on settings file
def bam_processing(files, sampleID, protocol):
    PATH = DIR_methcall+ "methylKit/"
    if len(files) == 1:
        # ---- single end
        return [PATH+"tabix_cpg/"+sampleID+"_se_bt2.sorted" + dedupe_tag(protocol) + "_cpg.txt.bgz",
                PATH+"tabix_cpg/"+sampleID+"_se_bt2.sorted" + dedupe_tag(protocol) + "_cpg.txt.bgz.tbi",
                # PATH+"tabix_chg/"+sampleID+"_se_bt2.sorted" + dedupe_tag(protocol) + "_chg.txt.bgz",
                # PATH+"tabix_chg/"+sampleID+"_se_bt2.sorted" + dedupe_tag(protocol) + "_chg.txt.bgz.tbi",
                # PATH+"tabix_chh/"+sampleID+"_se_bt2.sorted" + dedupe_tag(protocol) + "_chh.txt.bgz",
                # PATH+"tabix_chh/"+sampleID+"_se_bt2.sorted" + dedupe_tag(protocol) + "_chh.txt.bgz.tbi"
                ]
    elif len(files) == 2:
        # ---- paired end
        return [PATH+"tabix_cpg/"+sampleID+"_1_val_1_bt2.sorted" + dedupe_tag(protocol) + "_cpg.txt.bgz",
                PATH+"tabix_cpg/"+sampleID+"_1_val_1_bt2.sorted" + dedupe_tag(protocol) + "_cpg.txt.bgz.tbi",
                # PATH+"tabix_chg/"+sampleID+"_1_val_1_bt2.sorted" + dedupe_tag(protocol) + "_chg.txt.bgz",
                # PATH+"tabix_chg/"+sampleID+"_1_val_1_bt2.sorted" + dedupe_tag(protocol) + "_chg.txt.bgz.tbi",
                # PATH+"tabix_chh/"+sampleID+"_1_val_1_bt2.sorted" + dedupe_tag(protocol) + "_chh.txt.bgz",
                # PATH+"tabix_chh/"+sampleID+"_1_val_1_bt2.sorted" + dedupe_tag(protocol) + "_chh.txt.bgz.tbi"
                ]
    else:
        raise Exception(
            "=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")
                
# FIXME: contexts should be generate output based on settings file
def list_files_methyldackel_extract(files, sampleID, protocol):
    PATH = DIR_methcall + "methylDackel/"
    return [PATH+sampleID+dedupe_tag(protocol)+"_methyldackel_CpG.methylKit",
            PATH+sampleID+dedupe_tag(protocol)+"_methyldackel_CHG.methylKit",
            PATH+sampleID+dedupe_tag(protocol)+"_methyldackel_CHH.methylKit",
            # PATH+sampleID+"_mbias_methyldackel.txt",
            # PATH+sampleID+"_mbias_OB.svg",
            # PATH+sampleID+"_mbias_OT.svg",
            # PATH+sampleID+"_methyldackel.cytosine_report.txt"
            ]

# FIXME: contexts should be generate output based on settings file
def list_files_maketabix_methyldackel(files, sampleID, protocol):
    PATH = DIR_methcall + "methylDackel/"
    return [
            PATH+"tabix_CpG/"+sampleID+dedupe_tag(protocol)+"_CpG.txt.bgz",
            PATH+"tabix_CpG/"+sampleID+dedupe_tag(protocol)+"_CpG.txt.bgz.tbi",
            PATH+"tabix_CpG/"+sampleID+dedupe_tag(protocol)+"_CpG.txt.bgz",
            PATH+"tabix_CpG/"+sampleID+dedupe_tag(protocol)+"_CpG.txt.bgz.tbi",
            PATH+"tabix_CpG/"+sampleID+dedupe_tag(protocol)+"_CpG.txt.bgz",
            PATH+"tabix_CpG/"+sampleID+dedupe_tag(protocol)+"_CpG.txt.bgz.tbi"
            ]

# FIXME: contexts should be generate output based on settings file
def bigwig_exporting_bismark(files, sampleID, protocol):
    PATH = DIR_bigwig
    DESTRAND = "_destranded" if destrand("cpg") else ""
    if len(files) == 1:
        # ---- single end
        return PATH+sampleID+"_se_bt2.sorted" + dedupe_tag(protocol) + ".cpg" + DESTRAND + "_methylKit"+ ".bw"
    elif len(files) == 2:
        # ---- paired end
        return [PATH+sampleID+"_1_val_1_bt2.sorted" + dedupe_tag(protocol) + ".cpg" + DESTRAND +  "_methylKit" + ".bw"]
    else:
        raise Exception(
            "=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")

# FIXME: contexts should be generate output based on settings file
def bigwig_exporting_bwameth(files, sampleID, protocol):
    PATH = DIR_bigwig
    DESTRAND = "_destranded" if destrand("cpg") else ""
    if len(files) == 1:
        # ---- single end
        return [PATH+sampleID+dedupe_tag(protocol) + ".CpG" + DESTRAND + "_methylDackel" + ".bw"]
    elif len(files) == 2:
        # ---- paired end
        return [PATH+sampleID+dedupe_tag(protocol) + ".CpG" + DESTRAND + "_methylDackel" + ".bw"]
    else:
        raise Exception(
            "=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")

def methSeg_bismark(files, sampleID, protocol):
    PATH = DIR_seg
    if len(files) == 1:
        # ---- single end
        return PATH+sampleID+"_se_bt2.sorted" + dedupe_tag(protocol) + "_cpg" + "_methylKit" +".meth_segments.bed"
    elif len(files) == 2:
        # ---- paired end
        return [PATH+sampleID+"_1_val_1_bt2.sorted" + dedupe_tag(protocol) + "_cpg" + "_methylKit" +".meth_segments.bed"]
    else:
        raise Exception(
            "=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")

def methSeg_bwameth(files, sampleID, protocol):
    PATH = DIR_seg
    if len(files) == 1:
        # ---- single end
        return []
    elif len(files) == 2:
        # ---- paired end
        return [PATH+sampleID+dedupe_tag(protocol) + "_CpG" + "_methylDackel" +".meth_segments.bed"]
    else:
        raise Exception(
            "=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")

def list_final_reports_bwameth(files, sampleID, protocol):
    SUFFIX = "CpG_methylDackel_{}".format(ASSEMBLY) 
    PATH = os.path.join(DIR_final,"sample_reports/" )
    if len(files) == 1:
        # ---- single end
        # return [PATH+sampleID+dedupe_tag(protocol) + ".CpG" + DESTRAND + "_methylDackel" + ".bw"]
        return [PATH+sampleID+dedupe_tag(protocol) + "_"+ SUFFIX+ "_final.html"]
    elif len(files) == 2:
        # ---- paired end
        return [PATH+sampleID+dedupe_tag(protocol) + "_"+ SUFFIX+ "_final.html"]
    else:
        raise Exception(
            "=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")

def list_final_reports_bismark(files, sampleID, protocol):
    SUFFIX = "cpg_methylKit_{}".format(ASSEMBLY) 
    PATH = os.path.join(DIR_final,"sample_reports/" )
    if len(files) == 1:
        # ---- single end
        return PATH+sampleID+"_se_bt2.sorted" + dedupe_tag(protocol) + "_"+ SUFFIX+ "_final.html"
    elif len(files) == 2:
        # ---- paired end
        return [PATH+sampleID+"_1_val_1_bt2.sorted" + dedupe_tag(protocol) + "_"+SUFFIX+"_final.html"]
    else:
        raise Exception(
            "=== ERROR: file list is neither 1 nor 2 in length. STOP! ===")


# --------------------------------------
# diffmeth related
# --------------------------------------
def get_sampleids_from_treatment(treatment):
    """Get SampleIDs from treatment string."""
    sample_ids = list(config["SAMPLES"].keys())
    sample_treatments = [samplesheet(s,"Treatment") for s in sample_ids]

    sampleids_list = [sample_ids[i]
                 for i, x in enumerate(sample_treatments) if x == treatment]

    return(sampleids_list)
    

def get_sampleids_from_analysis(analysis):
    """Get SampleIDs for each Analysis group."""
    sampleids_list = []
    for group in config['DManalyses'][analysis]:
        for treatment in config['DManalyses'][analysis][group].split(","):
            sampleids_list += get_sampleids_from_treatment(treatment)
            
    return(sampleids_list)
    



def makeDiffMethPath(path, suffix, treatment):
    return path + str(treatment).replace('vs', '_') + dedupe_tag(config["SAMPLES"][get_sampleids_from_treatment(treatment[0])[0]]['Protocol']) + '_' + suffix

def diffmeth_input_function(treatments):
    treatments = treatments.replace(".deduped", "")
    sampleids = get_sampleids_from_treatment(treatments)

    inputfiles = []
    for sampleid in sampleids:
        fqname = config["SAMPLES"][sampleid]['fastq_name']
        protocol = config["SAMPLES"][sampleid]['Protocol']
        if len(fqname) == 1:
            inputfile = [os.path.join(
                DIR_methcall, sampleid+"_se_bt2.sorted" + dedupe_tag(protocol) + "_methylRaw.RDS")]
        elif len(fqname) == 2:
            inputfile = [os.path.join(
                DIR_methcall, sampleid+"_1_val_1_bt2.sorted" + dedupe_tag(protocol) + "_methylRaw.RDS")]
        inputfiles.append(inputfile)

    inputfiles = list(sum(inputfiles, []))
    return(inputfiles)

# FIXME: create named treatment groups in settings file
def files_for_treatment(proc):
    if "DManalyses" in config.keys():
        treatment_groups = config['DManalyses']
        if treatment_groups:
            return [expand(proc(comparison)) for comparison in treatment_groups if comparison]
        else:
            return []

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
