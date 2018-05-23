#' convert annotated flowSet to data.frame
#'
#' @param fcs_data an object of class flowSet
#' @param channels (default = c("Y1.A", "B1.A)) the column ratio will be the
#' ratio of Y1.A / B1.A
#' @return a data.frame
#' @examples
#' \dontrun{
#' data <- flowset2dataframe(fsc_data)
#' }
#' @importFrom flowCore sampleNames pData phenoData
#' @importFrom flowWorkspace GatingSet
#' @export flowset2dataframe
flowset2dataframe <- function(fsc_data, channels = c("Y1.A", "B1.A")) {
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
  data <- compute_line_column(data)
  data <- batch_effect(data)
  return(data)
}

compute_ratio <- function(data, channels) {
  b_ratio_1 <- colnames(data) %in% channels[1]
  b_ratio_2 <- colnames(data) %in% channels[2]
  data$ratio <- as.vector(data[, b_ratio_1] / data[, b_ratio_2])
  return(data)
}

parse_drug <- function(data) {
  data$drug_status <- as.factor(data$drug)
  b_drug <- !(data$drug %in% "None")
  data$drug[b_drug] <- paste0(as.vector(data$drug[b_drug]), "_",
                              as.vector(data$code.well[b_drug]))
  data$drug <- as.factor(data$drug)
  data$drug <- relevel(data$drug, "None")
  return(data)
}

compute_line_column <- function(data) {
  well <- as.vector(data$code.well)
  data$line <- gsub("([A-Z])[0-9]{2}", "\\1", well)
  data$column <- gsub("[A-Z]([0-9]{2})", "\\1", well)
  return(data)
}

batch_effect <- function(data) {
  b_drug <- data$drug %in% "None"
  well_number <- as.numeric(as.factor(data$code.well[!b_drug]))
  drug_number <- as.numeric(as.factor(data$code.well[b_drug]))
  well_number <- as.numeric(levels(as.factor(well_number)))
  drug_number <- as.numeric(levels(as.factor(drug_number)))
  none_dist <- matrix(
    data = rep(
      NA, length(well_number) * length(drug_number)
    ),
    ncol = length(drug_number)
  )
  j <- 1
  for (i in drug_number) {
    none_dist[, j] <- abs(well_number - i)
    j <- j + 1
  }
  none_closest <- apply(none_dist, 1, FUN = function(x){
    which(x %in% min(x))[1]
  })
  data$batch <- drug_number[none_closest]
  return(data)
}
