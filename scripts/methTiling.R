#!/usr/bin/env Rscript

# PiGx BSseq Pipeline.
#
# Copyright © 2026 Alexander Gosdschan <alexander.gosdschan@mdc-berlin.de>
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

# methTiling.R - tiles methylKit tabix methylation calls

## Collect arguments
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if (length(args) < 1) {
  args <- c("--help")
}

## Help section
if ("--help" %in% args) {
  cat(
    "
      Tile methylation calls using methylKit

      Arguments:
      --tabix location of input tabix file
      --outTabix name of tiled tabix output file
      --sampleId unique name of input sample
      --assembly assembly used to map the reads
      --context context of methylation
      --winSize tile window size in base pairs (default: 500)
      --stepSize tile step size in base pairs (default: 250)
      --covBases minimum covered bases per tile (default: 1)
      --cores number of processing cores
      --logFile file to print the logs to
      --help              - print this text

      Example:
      ./methTiling.R --tabix='input.txt.bgz' --outTabix='tiles.txt.bgz' --sampleId='sample' --assembly='hg19' --context='CpG' --cores=1 --logFile='tiles.log' \n\n"
  )

  q(save = "no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")

if (any(!grepl("^--[^=]+=.+$", args))) {
  stop("Arguments must use the form --arg=value")
}

argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1

requiredArgs <- c(
  "tabix",
  "outTabix",
  "sampleId",
  "assembly",
  "context",
  "cores",
  "logFile"
)
missingArgs <- setdiff(requiredArgs, names(argsL))
if (length(missingArgs) > 0) {
  stop("Missing required arguments: ", paste(missingArgs, collapse = ", "))
}

## catch output and messages into log file
out <- file(argsL$logFile, open = "at")
sink(out, type = "output")
sink(out, type = "message")

message("========= Given Arguments ==========")
print(argsL)
message("====================================")


# Run Functions -----------------------------------------------------------

### Methylation Tiling

## load methylKit
suppressPackageStartupMessages(expr = {
  library("methylKit")
})

tabixfile <- argsL$tabix
outTabix <- argsL$outTabix
sample.id <- argsL$sampleId
assembly <- argsL$assembly
# Coverage filtering is already applied by the methylation-calling rule.
# Tiling therefore only needs its own cov.bases threshold.
win.size <- as.integer(ifelse(is.null(argsL$winSize), 500, argsL$winSize))
step.size <- as.integer(ifelse(is.null(argsL$stepSize), 250, argsL$stepSize))
cov.bases <- as.integer(ifelse(is.null(argsL$covBases), 1, argsL$covBases))
cores <- as.integer(argsL$cores)

if (!file.exists(tabixfile)) {
  stop("Input tabix file does not exist: ", tabixfile)
}
if (
  any(is.na(c(win.size, step.size, cov.bases, cores))) ||
    any(c(win.size, step.size, cov.bases, cores) < 1)
) {
  stop(
    "--winSize, --stepSize, --covBases, and --cores must be positive integers"
  )
}

## format context string for methylKit's accepted context values
contextStr <- switch(
  tolower(argsL$context),
  "cpg" = "CpG",
  "chg" = "CHG",
  "chh" = "CHH",
  NULL
)
if (is.null(contextStr)) {
  stop(
    "The given context <",
    argsL$context,
    "> is not among the supported ones ('CpG','CHG','CHH')"
  )
}

## Read methylation calls from the existing, coverage-filtered tabix file
message("Reading tabix methylation calls from: ", tabixfile)
methylRawDB <- methRead(
  location = tabixfile,
  sample.id = sample.id,
  assembly = assembly,
  treatment = 0,
  context = contextStr,
  dbtype = "tabix"
)

## Tile methylation calls and save methylKit's generated tabix database
tileSuffix <- sprintf(
  "tiles_win%dbp_step%dbp_covBase%d",
  win.size,
  step.size,
  cov.bases
)
# methylKit derives its own filename.  Keep that intermediate database out of
# the public sample directory, then copy the data and index to Snakemake's
# parameterized output paths below.
workdir <- file.path(dirname(outTabix), ".methTiling", sample.id)
dir.create(workdir, recursive = TRUE, showWarnings = FALSE)

message("Tiling methylation calls")
methylTilesDB <- tileMethylCounts(
  methylRawDB,
  win.size = win.size,
  step.size = step.size,
  cov.bases = cov.bases,
  save.db = TRUE,
  mc.cores = cores,
  suffix = tileSuffix,
  dbdir = workdir
)

## Copy methylKit output to the declared Snakemake output paths.
## Both files are required because downstream methylKit reads the tabix index.
sourceTabix <- getDBPath(methylTilesDB)
sourceIndex <- paste0(sourceTabix, ".tbi")
outIndex <- paste0(outTabix, ".tbi")
if (!file.exists(sourceTabix) || !file.exists(sourceIndex)) {
  stop("methylKit did not produce the expected tiled tabix data and index")
}

dir.create(dirname(outTabix), recursive = TRUE, showWarnings = FALSE)
if (
  !file.copy(sourceTabix, outTabix, overwrite = TRUE) ||
    !file.copy(sourceIndex, outIndex, overwrite = TRUE)
) {
  stop("Could not copy tiled tabix output to declared destination: ", outTabix)
}
if (!all(file.exists(c(outTabix, outIndex)))) {
  stop("Expected tiled tabix output and index were not created: ", outTabix)
}

message("Tiled tabix saved to: \n\t", outTabix)
