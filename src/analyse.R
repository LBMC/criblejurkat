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
data_dir <- "data/examples/"

get_files <- function(path, regexp) {
  file_list <- base::list.files(
    path = path, full.names = TRUE, recursive = TRUE
  )
  file_list <- file_list[grepl(regexp, file_list, perl = T)]
  return(as.vector(unlist(file_list)))
}


load_annotation <- function(data_path) {
  annotation_path <- paste0(data_path, "annotation.csv")
  annotation <- read.table(annotation_path, h = T, sep = ";", stringsAsFactors = F)
  annotation$dapi <- as.factor(annotation$dapi)
  annotation$drug.time <- factor(paste(annotation$drug,".",  annotation$time, "UT", sep = ""))
  return(annotation)
}
annotation <- load_annotation(data_dir)

if (base::file.info(data_dir)$isdir) {
  fcs_files <- get_files(data_dir, ".fcs")
  x <- read.flowSet(fcs_files)
}
str(x)
pData(x) <- cbind(pData(x), annotation)
str(x)
