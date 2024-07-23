# PiGx BSseq Pipeline.
#
# Copyright © 2022 Alexander blume <alexander.blume@mdc-berlin.de>
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
import csv
from glob import glob
from func_defs import bail, warn

# ==============================================================================
#
#                                       VALIDATION FUNCTIONS
#
# functions used to verify the integrity of sample sheet and settings file   
# ==============================================================================


# --------------------------------------
# validation of sample sheet
# --------------------------------------

def get_filenames(mylist):
    return list(map(lambda x: splitext_fqgz(x)[0], mylist))

def splitext_fqgz(string):
    def fq_suffix(filename):
        return any(filename.endswith(ext) for ext in [".fq", ".fastq", ".fasta"])

    def is_zipped(filename):
        return any(filename.endswith(ext) for ext in [".gz", ".bz2"])

    if is_zipped(string):
        string, zipext = os.path.splitext(string)
    else:
        zipext = ""
    if fq_suffix(string):
        base, ext = os.path.splitext(string)
        return (base, ext + zipext)
    else:
        bail("Input files are not fastq files!")

def parse_samples(lines):
    """
    Parse csv table with information about samples, eg:

    Read1,Read2,SampleID,ReadType,Treatment
    sampleB.pe1.fq.gz,sampleB.pe2.fq.gz,sampleB,WGBS,B,,
    pe1.single.fq.gz,,sampleB1,WGBS,B,,

    It returns a dictionary required for the config file.
    """
    sreader = csv.reader(lines, delimiter=',')
    all_rows = [row for row in sreader if row]

    header = list(map(lambda x: x.strip(), all_rows[0] ))
    rows   = all_rows[1:]
    minimal_header = ['Read1', 'Read2', 'SampleID', 'Protocol', 'Treatment']

    if header[:5] != minimal_header:
        raise Exception("First columns of the input table have to be " +
                        ",".join(minimal_header) + ".")

    sample_ids = [x[2] for x in rows]
    if len(set(sample_ids)) != len(sample_ids):
        raise Exception("Column 'SampleID' has non-unique values.")

    # Create a dictionary with all params, keys are samples ids
    outputdict = {}
    for row in rows:
        if ( len(row) != 5 ):
            bail("Invalid row format in Samplesheet. Each row should have five columns.")
        row = list(map(lambda x: x.strip(), row))
        files = list(filter(None, row[0:2]))
        if not files:
            raise Exception("Each sample has to have an entry in at least one of the columns 'Read1' or 'Read2'.")

        sampleid_dict = {}
        for idx in range(len(header[2:])):
            try:
                sampleid_dict[header[2:][idx]] = row[2:][idx]
            except IndexError:
                raise Exception("Number of columns in row " + idx + " doesn't match number of elements in header.")

        sampleid_dict['files']      = files
        sampleid_dict['fastq_name'] = get_filenames(files)
        outputdict[row[2]] = sampleid_dict
    return { 'SAMPLES': outputdict }

def update_config_layout(config):
    # check for new annotation layout
    if not 'differential-methylation' in config['general']:
        return config

    differential_methylation = config['general']['differential-methylation']
    if not 'annotation' in differential_methylation:
        return config

    warn("The specification of annotation files has changed in version 0.1.9 \n"+
        "For now updating configuration file with the new default settings layout.\n"+
        "In the future, please retrieve the new default settings layout with 'pigx-bsseq --init settings'.")
    ## move cpgisland and refgene files to location dict
    if 'cpgIsland-bedfile' in differential_methylation['annotation']:
        config['locations']['cpgIsland-bedfile'] = differential_methylation['annotation']['cpgIsland-bedfile']
    if 'refGenes-bedfile' in differential_methylation['annotation']:
        config['locations']['refGenes-bedfile'] = differential_methylation['annotation']['refGenes-bedfile']
    ## remove annotation section
    differential_methylation.pop('annotation')

    config['general']['differential-methylation'] = differential_methylation

    return config



# --------------------------------------
# validation of Config 
# --------------------------------------
# check for common input/configuration errors:
def validate_config(config):
    # Check that all locations exist
    required_locations = ['input-dir', 'output-dir', 'genome-fasta']
    for loc in required_locations:
        if ( (not loc == 'output-dir') and (not (os.path.isdir(config['locations'][loc]) or os.path.isfile(config['locations'][loc])))):
            bail("ERROR: The following necessary directory/file does not exist: {} ({})".format(
                config['locations'][loc], loc))
    optional_locations = ['cpgIsland-bedfile', 'refGenes-bedfile']
    for loc in optional_locations:
        if ( config['locations'][loc] and ( not os.path.isfile(config['locations'][loc]))):
            bail("ERROR: The following optional file does not exist: {} ({})".format(
                config['locations'][loc], loc))
    


    # Load parameters specific to samples
    with open(config['locations']['sample-sheet'], 'r') as f:
            lines = f.read().splitlines()
    sample_params = parse_samples(lines)
    config.update(sample_params)

    # Check that all of the requested differential methylation
    # treatment values are found in the sample sheet.
    treatments = set([config["SAMPLES"][sample]["Treatment"]
                      for sample in config["SAMPLES"]])
    if 'DManalyses' in config:
        if config['DManalyses']:
            for analysis in config['DManalyses']:
                    for group in config['DManalyses'][analysis]['treatment_sample_groups'].split(",") + config['DManalyses'][analysis]['control_sample_groups'].split(","):
                        group = group.strip() #remove any leading/trailing whitespaces in the sample group names
                        if not any(treat == group for treat in treatments):
                            bail("ERROR: Invalid treatment group '{}' in analysis '{}'".format(
                            group, analysis))
        else:
            bail("ERROR: The config file contains empty 'DManalyses' section, please consider removing or commenting out this section.\n")
            
                        
    if 'treatment-groups' in config['general']['differential-methylation']:
        bail("ERROR: The specification of treatment groups and differential analysis has changed.\n"+
        "Please retrieve the new default settings layout with 'pigx-bsseq --init settings'.\n")

    # Check for a any Assembly string
    if not config['general']['assembly']:
            bail("ERROR: Please set a genome assembly string in the settings file at general::assembly.")

    # Check for a any Assembly string
    if not (config['general']['use_bwameth'] or config['general']['use_bismark']):
            bail("ERROR: Please enable one or both bisulfite aligners at general::use_bwameth/use_bismark.")
    

    # Check if we have permission to write to the reference-genome directory ourselves
    # if not, then check if the ref genome has already been converted
    genome_dir = os.path.dirname(config['locations']['genome-fasta'])
    if (not os.access(genome_dir, os.W_OK) and
            not os.path.isdir(os.path.join(genome_dir, 'Bisulfite_Genome'))):
        bail("ERROR: reference genome has not been bisulfite-converted, and PiGx does not have permission to write to that directory. Please either (a) provide Bisulfite_Genome conversion directory yourself, or (b) enable write permission in {} so that PiGx can do so on its own.".format(
            genome_dir))

    # Check for a genome fasta file
    fasta = glob(os.path.join(genome_dir, '*.fasta'))
    fa    = glob(os.path.join(genome_dir, '*.fa'))
    if not len(fasta) + len(fa) == 1 :
        bail("ERROR: Missing (or ambiguous) reference genome: The number of files ending in either '.fasta' or '.fa' in the following genome directory does not equal one: {}".format(genome_dir))
