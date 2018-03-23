load_infos <- function(annotation) {
  info <- t(unname(sapply(unique(annotation$cell.line), function(x) {
    annotation.temp <- annotation[which(annotation$cell.line == x), ]
    n <- dim(annotation.temp)[1]
    tab.time <- table(annotation.temp$time)
    c(x, n, table(annotation.temp$dapi)[names(table(annotation.temp$dapi))[1]], table(annotation.temp$drug)[names(table(annotation.temp$drug))[1]], paste(names(tab.time), " UT (", tab.time, ")", sep = "", collapse = " / "))
  })))
  colnames(info) <- c("cell.line", "nb samples", paste("DAPI ", names(table(annotation$dapi))[1], sep = ""), names(table(annotation$drug))[1],"drug time UT (nb samples)")
}

load_crap <- function(annotation) {
  column.names <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
  row.names <- c("A", "B", "C", "D", "E", "F", "G", "H")
  mat.value.template <- matrix(NA, nrow = 8, ncol = 12, dimnames = list(row.names, column.names))
  column.template <- matrix(1:96, nrow = 1, dimnames = list(NULL, unlist(sapply(column.names, function(x) paste(row.names, x, sep = ""), simplify = F))))
  row.template <- matrix(1:96, nrow = 1, dimnames = list(NULL, unlist(sapply(row.names, function(x) paste(x, column.names, sep = ""), simplify = F))))
  levels.fact <- levels(annotation$drug.time)
  n.cl <- length(unique(annotation$cell.line))
}

#' read annotation.csv files in a folder
#'
#' @param data_path path to the csv folder
#' @return a data.frame
#' @examples
#' \dontrun{
#' load_annotation("data/examples/")
#' }
#' @export load_annotation
load_annotation <- function(data_path) {
  annotation_path <- paste0(data_path, "annotation.csv")
  annotation <- read.table(annotation_path, h = T, sep = ";", stringsAsFactors = F)
  annotation$dapi <- as.factor(annotation$dapi)
  annotation$drug.time <- factor(paste(annotation$drug,".",  annotation$time, "UT", sep = ""))
  return(annotation)
}

load_channel <- function(data_path) {
  channel_path <- paste0(data_path, "channels.csv")
  channels <- read.table(channel_path, h = T, sep = ";", stringsAsFactors = F)
}

#' get list of file from a folder
#'
#' @param data_path path to the data folder
#' @return list of file names
#' @examples
#' \dontrun{
#' get_files("data/examples/")
#' }
get_files <- function(path, regexp) {
  file_list <- base::list.files(
    path = path, full.names = TRUE, recursive = TRUE
  )
  file_list <- file_list[grepl(regexp, file_list, perl = T)]
  file_list <- as.vector(unlist(file_list))
  if (length(file_list) == 0) {
    stop(print0(
      "error: ", data_path, " contain no ", regexp, " files."
    ))
  }
  return(file_list)
}

#' read fcs files in a folder
#'
#' @param data_path path to the fcs folder
#' @return a object of class flowSet
#' @examples
#' \dontrun{
#' load_data("data/examples/")
#' }
#' @importFrom flowCore read.flowSet
#' @export load_data
load_data <- function(data_path) {
  if (base::file.info(data_path)$isdir) {
    fcs_files <- get_files(data_path, ".fcs")
    fcs_data <- read.flowSet(
      fcs_files,
      transformation = F,
      alter.names = T
    )
    pData(fcs_data) <- cbind(
      pData(phenoData(fcs_data)),
      load_annotation(data_path)
    )
  } else {
    stop(print0(
      "error: ", data_path, " is not a directory"
    ))
  }
  return(fcs_data)
}
