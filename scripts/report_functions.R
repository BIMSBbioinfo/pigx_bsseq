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
  
  rmarkdown::render(
    input = reportFile,
    output_dir = outDir,
    # intermediates_dir = paste0(outDir,"/tmp"),
    # intermediates_dir = outDir,
    output_file = outFile,
    # knit_root_dir = rootDir,
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
  unlink(paste0(outDir,"/tmp"),recursive = TRUE)
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



render2Report(reportFile = normalizePath(snakemake@input[["template"]]),
              outFile = basename(snakemake@output[["report"]]),
              outDir = normalizePath(dirname(snakemake@output[["report"]])),
              #report.params = snakeParams )
              report.params = snakemake@params[nchar(names(snakemake@params)) > 0] )

#load("snakemakeObj.RData")