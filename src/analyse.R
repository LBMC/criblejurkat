require(quantmod)
require(flowStats)
require(ggcyto)
require(MASS)
require(gridExtra)
require(gplots)
source("../../src/func/functions.R")

source("https://bioconductor.org/biocLite.R")
biocLite("flowCore")

setwd("~/projects/CribleJurkat")
require(flowCore)
data_dir <- "data/2018_01_08_20171212/"

get_files <- function(path, regexp) {
  file_list <- base::list.files(
    path = path, full.names = TRUE, recursive = TRUE
  )
  file_list <- file_list[grepl(regexp, file_list, perl = T)]
  return(as.vector(unlist(file_list)))
}

if (base::file.info(data_dir)$isdir) {
  fcs_files <- get_files(data_dir, ".fcs")
  x <- read.FCS(fcs_files[1], transformation = FALSE)
}
summary(x)

read.flowSet(fcs_files)
