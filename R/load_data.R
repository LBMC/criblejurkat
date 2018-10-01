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
#' @importFrom utils read.table
#' @examples
#' \dontrun{
#' load_annotation("data/examples/")
#' }
#' @export load_annotation
load_annotation <- function(data_path) {
  annotation_path <- paste0(data_path, "annotation.csv")
  annotation <- utils::read.table(annotation_path, h = T, sep = ";", stringsAsFactors = F)
  annotation$dapi <- as.factor(annotation$dapi)
  annotation$drug.time <- factor(paste(annotation$drug,".",  annotation$time, "UT", sep = ""))
  return(annotation)
}

#' @importFrom utils read.table
load_channel <- function(data_path) {
  channel_path <- paste0(data_path, "channels.csv")
  channels <- utils::read.table(channel_path, h = T, sep = ";", stringsAsFactors = F)
}

#' get list of file from a folder
#'
#' @param path path to the data folder
#' @param regexp regexp that match the file names
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
    stop(paste0(
      "error: ", path, " contain no ", regexp, " files."
    ))
  }
  return(file_list)
}

#' read fcs files in a folder
#'
#' @param data_path path to the fcs folder
#' @return an object of class flowSet
#' @examples
#' \dontrun{
#' load_data("data/examples/")
#' }
#' @importFrom flowCore read.flowSet
#' @importFrom Biobase AnnotatedDataFrame
#' @importFrom utils read.csv
#' @importFrom methods as
#' @export load_data
load_data <- function(data_path) {
  if (base::file.info(data_path)$isdir) {
    fcs_files <- get_files(data_path, ".fcs")
  } else {
    stop(paste0(
      "error: ", data_path, " is not a directory"
    ))
  }
  if (base::length(fcs_files) < 1) {
    stop(paste0(
      "error: ", data_path, " doesn't contain .fcs files"
    ))
  }
  if (base::file.exists(paste0(data_path, "/annotation.csv"))) {
    annotation <- utils::read.csv(
      paste0(data_path, "/annotation.csv"), sep = ";", header = TRUE
    )
  } else {
    annotation <- tryCatch({
      annotation_files <- get_files(data_path, ".csv")
      print(paste0(
        "warning: no annotation.csv file found. Loading ",
        annotation_files[1]
      ))
      utils::read.csv(
        annotation_files[1], sep = ";", header = TRUE
      )
    }, error = function(e){
      stop(paste0(
        "error: load_data(). no annotation.csv file found in ",
        data_path
      ))
    })
  }
  annotation <- methods::as(annotation, "AnnotatedDataFrame")
  rownames(annotation) <- fcs_files
  fcs_data <- flowCore::read.flowSet(
    transformation = F,
    alter.names = T,
    phenoData = annotation,
    truncate_max_range = TRUE
  )
  return(fcs_data)
}

#' get project name
#'
#' @param fcs_data an object of class flowSet
#' @return the name of the project
#' @importFrom flowCore pData
#' @examples
#' \dontrun{
#' project_name(fcs_data)
#' }
#' @export project_name
project_name <- function(fcs_data) {
  if (is.data.frame(fcs_data)) {
    well <- fcs_data$file_path
  } else {
    well <- rownames(flowCore::pData(fcs_data))
  }
  gsub("data/(.+)/.*fcs", "\\1", well, perl=T)[1]
}

#' create outdir
#'
#' @param fcs_data an object of class flowSet
#' @param folder a directory name
#' @return the name of the project
#' @examples
#' \dontrun{
#' mk_outdir(fcs_data, "gating")
#' }
#' @export project_name
mk_outdir <- function(fcs_data, folder) {
  outdir <- paste0(
    "results/",
    project_name(fcs_data),
    "/", folder, "/"
  )
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  return(outdir)
}

#' convert data from anova_rlm to well position in chemical data base
#'
#' @param data data.frame outputed by anova_rlm()
#' @param col (default: c("01", "12")) new position of the "None" columns
#' @return data.frame with the None well in col
#' @examples
#' \dontrun{
#' data_chem <- move_none_well(data)
#' }
#' @export project_name
move_none_well <- function(data, col = c("01", "12")) {
  data$column <- as.factor(data$column)
  new_column_well <- levels(data$column)
  column_well_none <- levels( as.factor( as.vector(
    data$column[data$drug_status %in% "None"]
  )))
  for (i in 1:min(length(col), length(column_well_none))) {
    new_column_well <- move_column(new_column_well, column_well_none[i], col[i])
  }
  levels(data$column) <- new_column_well
  data <- data[order(as.numeric(as.factor(data$line)),
                     as.numeric(data$column)), ]
  data$code.well <- paste0(data$line, data$column)
  return(data)
}

column_name <- function(x) {
  if (is.numeric(x)) {
    x <- as.character(x)
    x[as.numeric(x) < 10] <- paste0("0", x[as.numeric(x) < 10])
  } else {
    x <- as.numeric(x)
  }
  return(x)
}

move_column <- function(x, old_col, new_col) {
    x <- as.numeric(x)
    new_col <- as.numeric(new_col)
    old_col <- as.numeric(old_col)
    old_col_pos <- which(x == old_col)
    if (old_col < new_col) {
      to_decrease <- old_col < x & x <= new_col
      x[to_decrease] <- x[to_decrease] - 1
    }
    if (old_col > new_col) {
      to_increase <- old_col > x & x >= new_col
      x[to_increase] <- x[to_increase] + 1
    }
    x[old_col_pos] <- new_col
    return(column_name(x))
}


