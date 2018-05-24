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
  data <- parse_drug(data)
  data <- compute_line_column(data)
  data <- batch_effect(data)
  data <- compute_ratio(data, channels)
  data <- power_trans(data)
  return(data)
}

#' @importFrom MASS boxcox
#' @importFrom scales boxcox_trans
power_trans <- function(data, formula = "ratio ~ drug + batch",
                        sample_size =  nrow(data)/100) {
  s_data <- data[sample(1:nrow(data), sample_size), ]
  model <- MASS::boxcox(as.formula(formula), data = s_data,
                        lambda = seq(-2, 10, 1/10))
  lambda <- model$x[model$y == max(model$y)]
  power_tr <- scales::boxcox_trans(lambda)
  power_tr <- power_trans$transform
  variable_name <- gsub("(.*) ~.*", "\\1", formula)
  data[[paste0(variable_name, "_norm")]] <-power_tr(data[[variable_name]])
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
  well_number <- as.numeric(as.factor(data$code.well))
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
  none_closest <- drug_number[none_closest]
  names(none_closest) <- 1:length(none_closest)
  data$batch <- as.factor(none_closest[as.numeric(as.factor(data$code.well))])
  return(data)
}

anova_rlm <- function(data, formula = "ratio ~ drug + batch") {
  variable_name <- gsub("(.*) ~.*", "\\1", formula)
  model <- MASS::rlm(as.formula(formula),
                     data = data,
                     psi = psi.huber,
                     k = quantile(data[[variable_name]], 0.90))
  outdir <- mk_outdir(fcs_data, "test")
  pdf(
    paste0(
      outdir,
      "anova_rlm.pdf"
    ),
    width = 80.3, height = 110.7
  )
  par(mfrow=c(4,2))
  plot(model)
  plot(1:nrow(data), residuals(model))
  plot(data$drug, residuals(model))
  plot(data$batch, residuals(model))
  dev.off()
  summodel <- summary(model)
  model_anova <- data.frame(summodel$coefficients)
  model_anova$p.value =  2*pt(
    abs(model_anova$t.value),
    summodel$df[2],
    lower.tail=FALSE
  )
  model_anova$signif <- model_anova$p.value < 0.0001
  write.csv(model_anova, file = paste0(outdir, "anova_rlm.csv"))
  data$signif <- NA
  data$coef <- NA
  data$coef_std <- NA
  data$pval <- NA
  data$tval <- NA
  for (drug in levels(data$drug)) {
    if (is.na(drug %in% "None")) {
      data$signif[data$drug %in% drug] <- model_anova$signif[grepl(drug, rownames(model_anova))]
      data$coef[data$drug %in% drug] <- model_anova$Value[grepl(drug, rownames(model_anova))]
      data$coef_std[data$drug %in% drug] <- model_anova[grepl(drug, rownames(model_anova)), 2]
      data$tval[data$drug %in% drug] <- model_anova$t.value[grepl(drug, rownames(model_anova))]
      data$pval[data$drug %in% drug] <- model_anova$p.value[grepl(drug, rownames(model_anova))]
    }
  }
  return(data)
}
