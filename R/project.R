#' perform all the analysis for a given data set
#'
#' @param data_path path to the data set
#' @param meta (default: FALSE) is the set_analysis part of a sets analysis
#' @return TRUE if everythings ran correclty
#' @examples
#' \dontrun{
#' set_analysis("data/2018_01_08_20171212/")
#' }
#' @export set_analysis
set_analysis <- function(data_path = "data/", meta = F, output_data = F) {
  fcs_raw <- load_data(data_path)
  data_raw <- flowset2dataframe(fcs_raw, norm = F)
  plot_well(data_raw, sufix = "_raw")
  plot_line(data_raw, sufix = "_raw")
  plot_column(data_raw, sufix = "_raw")
  fcs_raw <- remove_negatives(fcs_raw, data_raw)
  rm(data_raw)

  fcs_nonDebris <- rm_debris(fcs_raw)
  fcs_nonSinglets <- rm_nonsinglets(fcs_nonDebris)
  fcs_data <- rm_nonfluo(fcs_nonSinglets)

  data <- flowset2dataframe(fcs_data, norm = T)
  rm(fcs_data)
  if (output_data) {
    write.csv(data, file = paste0(
      "results/",
      gsub("data/(.+)", "\\1/", data_path, perl = T)[1],
      "data.csv"
    ))
  }
  if (meta) {
    return(data)
  }
  outdir_rlm <- paste0("results/",
                       gsub("data/(.+)", "\\1/", data_path, perl=T)[1])
  if (!dir.exists(outdir_rlm)) {
    dir.create(outdir_rlm, recursive = TRUE)
  }
  data <- anova_lm(data, outdir = outdir_rlm)
  data <- scaled_pval(data)
  plot_well(data)
  plot_line(data)
  plot_column(data)
  data_chem <- move_none_well(data)
  plot_well(data_chem, sufix = "_chem")
}

#' perform all the analysis for a set of data sets
#'
#' @param data_path path to folder containing the data sets
#' @param rlm_model (default: TRUE) should rlm model be use or lm ?
#' @param chunk (default: 20000) size of data chunk to do the computation on for lm
#' @return TRUE if everythings ran correclty
#' @examples
#' \dontrun{
#' analysis("data/set_test")
#' }
#' @export analysis
analysis <- function(data_path = "data/", rlm_model = FALSE, output_data = F) {
  if (base::file.info(data_path)$isdir) {
    set_folders <- list.dirs(data_path, full.names = F)[-1]
  } else {
    stop(paste0(
      "error: ", data_path, " is not a directory"
    ))
  }
  if (base::length(set_folders) < 1) {
    stop(paste0(
      "error: ", data_path, " doesn't contains any folders"
    ))
  }
  outdir_rlm <- paste0("results/",
                       gsub("data/(.+)", "\\1/", data_path, perl=T)[1])
  if (!dir.exists(outdir_rlm)) {
    dir.create(outdir_rlm, recursive = TRUE)
  }
  min_sets_factors <- c()
  for (folder in set_folders) {
    message(paste0("gating for ", data_path, "/", folder))
    if (!file.exists(paste0(outdir_rlm, "/", folder, ".Rdata"))) {
      set_data <- set_analysis(paste0(data_path, "/", folder),
                               meta = T,
                               output_data = output_data)
      set_data$set <- folder
      set_data$drug <- ifelse(set_data$drug %in% "None",
                              "None",
                              paste0(set_data$drug, "_", folder)
      )
      save(set_data,
           file = paste0(outdir_rlm, "/", folder, ".Rdata"))
    } else {
      message(paste0("gating file for ",
                     folder,
                     " found. skipping computation"))
      load(paste0(outdir_rlm, "/", folder, "/", folder, ".Rdata"))
    }
    if (length(min_sets_factors) == 0) {
      min_sets_factors <- colnames(set_data)
    }
    min_sets_factors <- intersect(min_sets_factors,
                                  colnames(set_data))
    rm(set_data)
  }
  data <- setNames(data.frame(matrix(ncol = length(min_sets_factors),
                                     nrow = 0)),
                   min_sets_factors)
  for (folder in set_folders) {
    load(paste0(outdir_rlm, "/", folder, "/", folder, ".Rdata"))
    data <- rbind(data, set_data[, min_sets_factors])
    rm(set_data)
  }
  data$set <- as.factor(data$set)
  if (rlm_model) {
    data <- anova_rlm(data, formula = "ratio ~ 1 + drug + batch + set",
              outdir = outdir_rlm)
  } else {
    data <- anova_lm(data, formula = "ratio ~ 1 + drug + batch + set",
              outdir = outdir_rlm, chunk = 200000)
  }
  data <- scaled_pval(data)
  for (folder in set_folders) {
    message(paste0("plotting for ", folder))
    plot_well(data[data$set %in% folder, ])
    plot_line(data[data$set %in% folder, ])
    plot_column(data[data$set %in% folder, ])
    data_chem <- move_none_well(data[data$set %in% folder, ])
    plot_well(data_chem, sufix = "_chem")
  }
}
