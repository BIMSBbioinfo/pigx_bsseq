# PiGx BSseq Pipeline.
#
# Copyright © 2018 Alexander Gosdschan <alexander.gosdschan@mdc-berlin.de>
# Copyright © 2019 Alexander Blume <alexander.Blume@mdc-berlin.de>
# Copyright © 2017, 2018 Bren Osberg <Brendan.Osberg@mdc-berlin.de>
# Copyright © 2017, 2018 Katarzyna Wreczycka <katwre@gmail.com>
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


# This file was inspired by
# https://github.com/katwre/makeWGBSnake/blob/master/Scripts/MethDiff.R
# taken from

# WGBS pipeline
#
# Copyright © 2018 Katarzyna Wreczycka katarzyna.wreczycka@mdc-berlin.de
# This pipeline is heavily based on the PiGx BSseq pipeline github.com/BIMSBbioinfo/pigx_bsseq


# methDiff.R - takes a methylKit methylBaseDB tabix file and performs
# 	differential methylation testing
# ---last updated Oct. 2019 by A. Blume

## Collect arguments
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if (length(args) < 1) {
  args <- c("--help")
}

## Help section
if ("--help" %in% args) {
  cat("
      Calculate Differential Methylation
      
      Arguments:
      --inputfile 	merged methylation call files generated by unite_meth_calls 
      --sampleids 	sample names used in this analysis 
      --treatments 	comma separated list of treatment values in same order as sample_ids
      --assembly 	assembly to be used in methlyDiff object
      --context		methylation context
      --destranded 	wether strands are merged or not
      --methylDiff_results_suffix 		suffix of tabix file containing testing results
      --treatment_group  comma separated list of treatment groups
      --control_group   comma separated list of control groups (put same as 
        treatment group again if you want to calculate methylation difference 
        as the difference of max(x) - min(x) where x is vector of 
        mean methylation per group per region)
      --resultsFile 	file containing testing results
      --outdir output directory
      --cores number of cores to use for calculateDiffMeth
      --logFile file to print the logs to
      --help              - print this text
      
      Example:
      ./test.R --arg1=1 --arg2='output.txt' --arg3=TRUE \n\n")

  q(save = "no")
}


parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
parseArgsList <- function(x) strsplit(as.character(x), " ")

## Parse arguments (we expect the form --arg=value)
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1
argsL <- sapply(argsL, parseArgsList)


## catch output and messages into log file
out <- file(argsL$logFile, open = "at")
sink(out, type = "output")
sink(out, type = "message")

message("========= Given Arguments ==========")
print(argsL)
message("====================================")

# Prepare Variables -----------------------------------------------------------

# load libraries
suppressPackageStartupMessages(library("methylKit"))
data.table::setDTthreads(8)

# load variables
input <- argsL$inputfile
sampleids <- strsplit(argsL$sampleids, ",", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]]
assembly <- argsL$assembly
destranded <- ifelse(tolower(argsL$destranded) %in% c("true", "yes"), TRUE, FALSE)
outdir <- argsL$outdir
resultsFile <- argsL$resultsFile

# split all treatment values, could be numeric or not
treatmentsStr <- strsplit(argsL$treatments, ",", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]]

treatment_group <- strsplit(argsL$treatment_group, ",", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]]
control_group <- strsplit(argsL$control_group, ",", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]]

if (!setequal(control_group, treatment_group)) {
  message("Remapping Treatments Descriptions into treatment and control groups.")
  # reorganize samples into treatment and control groups
  treatments <- ifelse(treatmentsStr %in% treatment_group, 1, 0)

  message("Control group (0): ", paste0(control_group, collapse = ", "))
  message("Treatment group (1): ", paste0(treatment_group, collapse = ", "))
} else {
  message("Remapping Treatment Description to number.")
  # convert into named numeric vector
  treatments <- as.numeric(as.factor(treatmentsStr))
  names(treatments) <- treatmentsStr

  message(paste(sort(unique(treatmentsStr)),
    sort(unique(treatments)),
    sep = ": ", collapse = "\n"
  ))
}


# convert variables and perform checks
context <- switch(tolower(argsL$context),
  cpg = "CpG",
  chh = "CHH",
  chg = "CHG",
  NULL
)

if (is.null(context)) {
  stop(
    "The given context <", argsL$context,
    "> is not among the supported ones ('CpG','CHG','CHH')"
  )
}

cores <- as.numeric(argsL$cores)

methylDiff_results_suffix <- argsL$methylDiff_results_suffix

# Run Functions -----------------------------------------------------------

# Read input files
message("Reading united samples ...")
methylBaseDB <- methylKit:::readMethylBaseDB(
dbpath = input,
dbtype = "tabix",
sample.ids = sampleids,
assembly = assembly,
context = context,
resolution = "base",
treatment = treatments,
destranded = destranded
)

# remove output file if already exists
methylDiff_outFile <- gsub("results.tsv",
						 sprintf("%s.txt.bgz",
								 methylDiff_results_suffix),
						 resultsFile)

if( file.exists(methylDiff_outFile)) {
unlink(c(methylDiff_outFile, paste0(methylDiff_outFile, ".tbi")))
}

# Find differentially methylated cytosines
message("Calculating differential methylation ...")
methylDiffDB <- methylKit::calculateDiffMeth(
.Object = methylBaseDB,
overdispersion = "MN",
test = "Chisq",
mc.cores = cores,
save.db = TRUE,
dbdir = outdir,
suffix = methylDiff_results_suffix
)


message("Exporting results to tsv...")
methylKit:::.write.table.noSci(methylKit::getData(methylDiffDB),
sep = "\t",
row.names = FALSE,
quote = FALSE,
file = resultsFile
)
