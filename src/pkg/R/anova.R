#' convert annotated flowSet to data.frame
#'
#' @param fcs_data an object of class flowSet
#' @param channels (default = c("Y1.A", "B1.A)) the column ratio will be the
#' ratio of Y1.A / B1.A
#' @param line_lenght (default = 12) length of the lines in the plate
#' @param column_length (default = 8) length of the columns in the plate
#' @return a data.frame
#' @examples
#' \dontrun{
#' data <- flowset2dataframe(fsc_data)
#' }
#' @importFrom flowCore sampleNames pData phenoData
#' @importFrom flowWorkspace GatingSet
#' @export flowset2dataframe
flowset2dataframe <- function(fsc_data, channels = c("Y1.A", "B1.A"),
                              line_length = 12,
                              column_length = 8) {
  data <- fcs_data
  data <- apply(matrix(flowCore::sampleNames(data), ncol = 1), 1,
    FUN = function(x, fset, infos){
      data <- data.frame(exprs(fset[[x]]))
      infos <- infos[rownames(infos) %in% x, ]
      data <- base::data.frame(
        x,
        1:nrow(data),
        data,
        infos,
        row.names = NULL
      )
      return(data)
    },
    fset = data,
    infos = flowCore::pData(flowCore::phenoData(fcs_data))
  )
  data <- data.frame(do.call(rbind, data))
  names(data)[1:2] <- c("well", "step")
  data$name <- as.factor(data$name)
  data$drug <- as.vector(data$drug)
  data <- compute_ratio(data, channels)
  data <- parse_drug(data)
  data <- compute_line_column(data, line_lenght, column_length)
  return(data)
}

compute_ratio <- function(data, channels) {
  b_ratio_1 <- colnames(data) %in% channels[1]
  b_ratio_2 <- colnames(data) %in% channels[2]
  data$ratio <- as.vector(data[, b_ratio_1] / data[, b_ratio_2])
  return(data)
}

parse_drug <- function(data) {
  data$drug_status <- data$drug
  b_drug <- !(data$drug %in% "None")
  data$drug[b_drug] <- paste0(as.vector(data$drug[b_drug]), "_",
                              as.vector(data$code.well[b_drug]))
  data$drug <- as.factor(data$drug)
  data$drug <- relevel(data$drug, "None")
  return(data)
}

compute_line_column <- function(data, line_lenght = 12, column_length = 8) {
  well_number <- as.numeric(as.factor(data$code.well))
  data$line <- floor(( well_number - 1 ) / 12)
  data$column <- floor(well_number / 12)
  return(data)
}

batch_effect <- function(data) {
  well_number <- as.numric(as.factor(data$code.well))
  b_drug <- !(data$drug %in% "None")
}
