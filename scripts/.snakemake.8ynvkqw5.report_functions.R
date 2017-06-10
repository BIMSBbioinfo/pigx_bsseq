
######## Snakemake header ########
library(methods)
Snakemake <- setClass(
    "Snakemake",
    slots = c(
        input = "list",
        output = "list",
        params = "list",
        wildcards = "list",
        threads = "numeric",
        log = "list",
        resources = "list",
        config = "list",
        rule = "character"
    )
)
snakemake <- Snakemake(
    input = list('/home/bosberg/bs/pigx_sandbox/report_templates/methCall.report.Rmd', '06_sorted/output2K_1_val_1_bt2.deduped.sorted.bam', "template" = '/home/bosberg/bs/pigx_sandbox/report_templates/methCall.report.Rmd', "bamfile" = '06_sorted/output2K_1_val_1_bt2.deduped.sorted.bam'),
    output = list('methylation_calls/output2K_1_val_1_bt2.deduped.sorted_meth_calls.nb.html', 'methylation_calls/output2K_1_val_1_bt2.deduped.sorted_methylRaw.RData', 'methylation_calls/output2K_1_val_1_bt2.deduped.sorted_CpG.txt', "report" = 'methylation_calls/output2K_1_val_1_bt2.deduped.sorted_meth_calls.nb.html', "callFile" = 'methylation_calls/output2K_1_val_1_bt2.deduped.sorted_CpG.txt', "rdatafile" = 'methylation_calls/output2K_1_val_1_bt2.deduped.sorted_methylRaw.RData'),
    params = list('ce10', 10, '/home/bosberg/bs/pigx_sandbox/test_dataset/out/06_sorted/output2K_1_val_1_bt2.deduped.sorted.bam', '/home/bosberg/bs/pigx_sandbox/test_dataset/out/methylation_calls/output2K_1_val_1_bt2.deduped.sorted_methylRaw.RData', 0, "assembly" = 'ce10', "minqual" = 10, "inBam" = '/home/bosberg/bs/pigx_sandbox/test_dataset/out/06_sorted/output2K_1_val_1_bt2.deduped.sorted.bam', "rdata" = '/home/bosberg/bs/pigx_sandbox/test_dataset/out/methylation_calls/output2K_1_val_1_bt2.deduped.sorted_methylRaw.RData', "mincov" = 0),
    wildcards = list('output2K_1_val_1_bt2.deduped', "prefix" = 'output2K_1_val_1_bt2.deduped'),
    threads = 1,
    log = list('/home/bosberg/bs/pigx_sandbox/test_dataset/out/output2K_1_val_1_bt2.deduped.sorted_meth_calls.log'),
    resources = list(),
    config = list("trim_galore_args" = '', "CHROM_INFO" = '/home/bosberg/bs/ref_genome/ce10/chromInfo.txt', "GENOMEPATH" = '/home/bosberg/bs/ref_genome/ce10/', "bam_methCall_args_mincov" = '0', "PATHIN" = '/home/bosberg/bs/pigx_sandbox/test_dataset/in/', "RCODE" = '2K_', "bam_methCall_args_minqual" = '10', "fastqc_args" = '', "INEXT" = '.fq.gz', "NON_DIR_FLAG" = '', "NUMTHREADS" = '2', "SAMPLES" = list("sampleA" = list("ReadType" = 'WGBS', "SampleID" = 'sampleA', "Treatment" = 'A', "fastq_ext" = c('fq.gz', 'fq.gz'), "fastq_name" = c('output2K_1', 'output2K_2'), "files" = c('output2K_1.fq.gz', 'output2K_2.fq.gz')), "sampleB" = list("ReadType" = 'WGBS', "SampleID" = 'sampleB', "Treatment" = 'B', "fastq_ext" = c('fq.gz'), "fastq_name" = c('output2K_1'), "files" = c('output2K_1.fq.gz'))), "bismark_args" = ' -N 1 -L 2 ', "LOG" = '/home/bosberg/bs/pigx_sandbox/test_dataset/out/', "NICE" = '19', "PATHOUT" = '/home/bosberg/bs/pigx_sandbox/test_dataset/out/', "PROGS" = list("BISMARK" = 'bismark', "BISMARK2REPORT" = 'bismark2report', "BISMARK_GENOME_PREPARATION" = 'bismark_genome_preparation', "BISMARK_METHYLATION_EXTRACTOR" = 'bismark_methylation_extractor', "BOWTIE2" = 'bowtie2', "CUTADAPT" = 'cutadapt', "DEDUPLICATE_BISMARK" = 'deduplicate_bismark', "FASTQC" = 'fastqc', "SAMTOOLS" = 'samtools', "TRIMGALORE" = 'trim_galore'), "GTOOLBOX" = '/home/bosberg/.guix-profile/bin/', "GENOME_VERSION" = 'ce10'),
    rule = 'bam_methCall'
)
######## Original script #########
## Wrapper function to run a specific Rmd script
## which does the computation and generates an HTML report with code included


render2Report <- function(reportFile,
                          outFile,
                          outDir,
                          report.params)
{
  
  #print(getwd())
  
  
  ## write stdout to log file
  # sink(snakemake@log[[1]])
  
 
  ## the logo is stored in the template directory
  pathToLogo <- paste0(normalizePath(dirname(reportFile)),"/pigx_bsseq_logo.html")
  
  ## we set the knitr root dir to be the base directory,
  ## such that all paths are relative from there
  rootDir <- dirname(dirname(reportFile))

  interDir <- paste0(outDir,"/inter")

  
  rmarkdown::render(
    input = reportFile,
    output_file = outFile,
    output_dir = outDir,
    # intermediates_dir = paste0(outDir,"/tmp"),
    intermediates_dir = interDir,
    knit_root_dir =  outDir,#rootDir
    output_format = rmarkdown::html_notebook(
      toc = TRUE,
      toc_float = TRUE,
      theme = 'lumen',
      number_sections = FALSE,
      code_folding = "hide",
      self_contained = TRUE,
      includes = list(in_header = pathToLogo)
    ),
    params = report.params,
    quiet = FALSE,
    clean = TRUE,
    envir = new.env()
  )
  #unlink(paste0(outDir,"/tmp"),recursive = TRUE)
  #unlink(list.files(path = outDir,pattern = "knit|utf8|nb_files"),recursive = TRUE)
  
  
  
  
}



## catch output and messages into log file
out <- file(snakemake@log[[1]], open = "wt")
sink(out,type = "output")
sink(out, type = "message")


## debugging
# save.image(file = "snakemakeObj.RData")

## check for filepaths to be normalized
# snakeParams <- snakemake@params[nchar(names(snakemake@params)) > 0]
# snakeParams <- lapply(snakeParams, function(x) {
#   if(class(x) == "character") return(normalizePath(x))
#   else return(x) })

cat(paste(  
  Sys.time(),"\n\n",
  "Rendering report:",basename(snakemake@output[["report"]]),"\n",
  "from template:",basename(snakemake@input[["template"]]),"\n",
  "into directory:",normalizePath(dirname(snakemake@output[["report"]])),"\n"
  ))


render2Report(reportFile = normalizePath(snakemake@input[["template"]]),
              outFile = basename(snakemake@output[["report"]]),
              outDir = normalizePath(dirname(snakemake@output[["report"]])),
              #report.params = snakeParams )
              report.params = snakemake@params[nchar(names(snakemake@params)) > 0] )



#load("snakemakeObj.RData")
