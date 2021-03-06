#' convert annotated flowSet to data.frame
#'
#' @param fcs_data an object of class flowSet
#' @param channels (default = c("Y1.A", "B1.A)) the column ratio will be the
#' ratio of Y1.A / B1.A
#' @param norm (default: TRUE) normalization step (power_trans) for the ratio
#' @return a data.frame
#' @examples
#' \dontrun{
#' data <- flowset2dataframe(fsc_data)
#' }
#' @importFrom flowCore sampleNames pData phenoData exprs
#' @importFrom flowWorkspace GatingSet
#' @export flowset2dataframe
flowset2dataframe <- function(fcs_data, channels = c("Y1.A", "B1.A"),
                              norm = TRUE) {
  data <- fcs_data
  data <- apply(matrix(flowCore::sampleNames(data), ncol = 1), 1,
    FUN = function(x, fset, infos){
      data <- data.frame(flowCore::exprs(fset[[x]]))
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
  names(data)[1:2] <- c("file_path", "step")
  data$name <- as.factor(data$name)
  data$drug <- as.vector(data$drug)
  data <- parse_drug(data)
  data <- compute_line_column(data)
  data <- batch_effect(data)
  data <- compute_ratio(data, channels)
  if (norm) {
    data <- power_trans(data)
  }
  return(data)
}

#' @importFrom MASS boxcox
#' @importFrom scales boxcox_trans
power_trans <- function(data, formula = "ratio ~ batch",
                        sample_size =  nrow(data)/100) {
  s_data <- data[sample(1:nrow(data), sample_size), ]
  model <- MASS::boxcox(stats::as.formula(formula), data = s_data,
                        lambda = seq(-2, 10, 1/10),
                        plotit = FALSE)
  lambda <- model$x[model$y == max(model$y)]
  power_tr <- scales::boxcox_trans(lambda)
  power_tr <- power_tr$transform
  variable_name <- gsub("(.*) ~.*", "\\1", formula)
  if (min(data[[variable_name]]) < 0) {
    data[[variable_name]] <- data[[variable_name]] +
      abs(min(data[[variable_name]]))
  }
  data[[paste0(variable_name, "_norm")]] <-power_tr(data[[variable_name]])
  return(data)
}

compute_ratio <- function(data, channels) {
  b_ratio_1 <- colnames(data) %in% channels[1]
  b_ratio_2 <- colnames(data) %in% channels[2]
  data$ratio <- as.vector(data[, b_ratio_1] / data[, b_ratio_2])
  return(data)
}

#' @importFrom stats relevel
parse_drug <- function(data) {
  data$drug_status <- as.factor(data$drug)
  b_drug <- !(data$drug %in% "None")
  data$drug[b_drug] <- paste0(as.vector(data$drug[b_drug]), "_",
                              as.vector(data$code.well[b_drug]))
  data$drug <- as.factor(data$drug)
  data$drug <- stats::relevel(data$drug, "None")
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

#' build lm model between drugs accounting for batch effect
#'
#' @param data a data.frame
#' @param formula (default: "ratio ~ drug + batch") the formula of the model
#' @param lower (default: TRUE) tests if "y" (the ratio) is lower than in controls
#' @param chunk (default: 20000) size of data chunk to do the computation on
#' @param outdir set the outdir directory to a path (default extracted from fcs)
#' @return a data.frame
#' @examples
#' \dontrun{
#' data <- anova_lm(data)
#' }
#' @import biglm
#' @importFrom grDevices dev.off pdf
#' @importFrom stats as.formula quantile
#' @export anova_lm
anova_lm <- function(data, formula = "ratio ~ drug + batch",
                     lower = TRUE,
                     outdir,
                     chunk = nrow(data)) {
  if (file.exists(paste0(outdir, "anova_lm.Rdata"))) {
    message("model file found. skipping computation")
    load(paste0(outdir, "anova_lm.Rdata"))
    model_anova <- compute_pval(model, lower = lower)
    data <- export_lm_results(data, model_anova)
  } else {
    data$drug <- as.factor(data$drug)
    data$drug <- relevel(data$drug, "None")
    if (nrow(data) > chunk) {
      model <- split_lm(data = data, formula = formula, chunk = chunk)
    } else {
      model <- biglm::biglm(stats::as.formula(formula),
                        data = data)
    }
    model_anova <- compute_pval(model, lower = lower)
    if (missing(outdir)) {
      outdir <- mk_outdir(data, "/test")
    }
    data <- export_lm_results(data, model_anova)
    save(data, model, file = paste0(outdir, "anova_lm.Rdata"))
  }
  export_drug_table(data, model_anova, outdir, model = "lm")
  return(data)
}

#' split lm model by chunck with biglm
#'
#' @param data a data.frame
#' @param formula (default: "ratio ~ drug + batch") the formula of the model
#' @param chunk (default: 200000) chunk size
#' @return a data.frame
#' @import biglm
#' @importFrom stats as.formula quantile
#' @export split_lm
split_lm <- function(data, formula = "ratio ~ drug + batch", chunk = nrow(data)) {
  print(paste0("spliting lm model by chunk of ", chunk, " lines"))
  data_index <- seq_len(nrow(data))
  data_index <- sample(data_index, nrow(data))
  index_start <- 1
  data$drug <- as.factor(data$drug)
  data$drug <- relevel(data$drug, "None")
  for (index_stop in seq(from = chunk, to = nrow(data), by = chunk)) {
    row_select <- data_index[index_start:index_stop]
    if (index_start == 1) {
      model <- biglm::biglm(stats::as.formula(formula),
                        data = data[row_select, ])
    } else {
      update(model, data[row_select, ])
    }
    index_start <- index_stop + 1
  }
  return(model)
}

#' build rlm model between drugs accounting for batch effect
#'
#' @param data a data.frame
#' @param formula (default: "ratio ~ drug + batch") the formula of the model
#' @param lower (default: TRUE) tests if "y" (the ratio) is lower than in controls
#' @param outdir set the outdir directory to a path (default extracted from fcs)
#' @return a data.frame
#' @examples
#' \dontrun{
#' data <- anova_rlm(data)
#' }
#' @importFrom MASS rlm psi.huber
#' @importFrom grDevices dev.off pdf
#' @importFrom stats as.formula quantile
#' @export anova_rlm
anova_rlm <- function(data, formula = "ratio ~ drug + batch", lower = TRUE,
                      outdir) {
  if (file.exists(paste0(outdir, "anova_rlm.Rdata"))) {
    message("model file found. skipping computation")
    load(paste0(outdir, "anova_rlm.Rdata"))
    model_anova <- compute_pval(model, lower = lower)
    data <- export_rlm_results(data, model_anova)
  } else {
    variable_name <- gsub("(.*) ~.*", "\\1", formula)
    model <- MASS::rlm(stats::as.formula(formula),
                      data = data,
                      psi = MASS::psi.huber,
                      k = stats::quantile(data[[variable_name]], 0.90))
    model_anova <- compute_pval(model, lower = lower)
    if (missing(outdir)) {
      outdir <- mk_outdir(data, "/test")
    }
    data <- export_rlm_results(data, model_anova)
    save(data, model, file = paste0(outdir, "anova_rlm.Rdata"))
  }
  export_drug_table(data, model_anova, outdir, model = "rlm")
  return(data)
}

#' @import biglm
#' @importFrom stats pt
compute_pval <- function(model, lower = TRUE) {
  model_anova <- data.frame(coef = coef(model),
                            se = sqrt(diag(vcov(model))),
                            t.value = coef(model) / sqrt(diag(vcov(model))),
                            row.names = names(coef(model))
  )
  model_df <- model$df
  if ("n" %in% names(model)) {
    model_df <- model$n - length(model$names)
  }
  model_anova$pval =  stats::pt(
    model_anova$t.value,
    model_df,
    lower.tail=TRUE
  )
  return(model_anova)
}

#' @importFrom utils write.csv
export_rlm_results <- function(data, model_anova) {
  data$coef <- NA
  data$coef_std <- NA
  data$pval <- NA
  data$tval <- NA
  data$drug <- as.factor(data$drug)
  for (drug in levels(data$drug)) {
    if (!(drug %in% "None")) {
      data$coef[data$drug %in% drug] <- model_anova$Value[grepl(drug,
        rownames(model_anova))]
      data$coef_std[data$drug %in% drug] <- model_anova[grepl(drug,
        rownames(model_anova)), 2]
      data$tval[data$drug %in% drug] <- model_anova$t.value[grepl(drug,
        rownames(model_anova))]
      data$pval[data$drug %in% drug] <- model_anova$pval[grepl(drug,
        rownames(model_anova))]
    }
  }
  return(data)
}

#' @importFrom utils write.csv
export_lm_results <- function(data, model_anova) {
  data$coef <- NA
  data$coef_std <- NA
  data$pval <- NA
  data$tval <- NA
  data$drug <- as.factor(data$drug)
  for (drug in levels(data$drug)) {
    if (!(drug %in% "None")) {
      data_select <- data$drug %in% drug
      model_select <- grepl(drug, rownames(model_anova))
      if (any(model_select)) {
        data$coef[data_select] <- model_anova$coef[model_select]
        data$coef_std[data_select] <- model_anova$se[model_select]
        data$tval[data_select] <- model_anova$t.value[model_select]
        data$pval[data_select] <- model_anova$pval[model_select]
      } else {
        message(paste0(drug ," not tested"))
      }
    }
  }
  return(data)
}

export_drug_table <- function(data, model_anova, outdir,
                              channels = c("Y1.A", "B1.A"),
                              model = "rlm") {
  drug_table <- model_anova[grepl("drug", rownames(model_anova)), ]
  drug_name <- rownames(model_anova)[grepl("drug", rownames(model_anova))]
  for (channel in channels) {
    drug_mean <- c(by(data[, which(colnames(data) %in% channel)],
                                  data$drug, mean))
    drug_test <- paste0("drug",names(drug_mean)) %in% rownames(drug_table)
    drug_table[[channel]] <- drug_mean[drug_test]
  }
  utils::write.csv(model_anova, file = paste0(outdir, "anova_", model,
                                              ".csv"))
  drug_table <- drug_table[order(drug_table$pval), ]
  utils::write.csv(drug_table, file = paste0(outdir, "anova_", model,
                                             "_drug.csv"))
}
