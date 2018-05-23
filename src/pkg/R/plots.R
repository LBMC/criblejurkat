#' plot the ratio distribution per well
#'
#' @param data data.frame object
#' @param sample_size (default: nrow(data) / 100) size of the subsample of data
#' to draw, to speedup the ploting proccess
#' @param sufix (default: "") sufix to append to pdf name
#' @return a ggplot2 object
#' @examples
#' \dontrun{
#' plot_plates_ratio(fsc_data)
#' }
#' @import ggplot2
#' @export plot_well_ratio
plot_well_ratio <- function(data, sample_size = nrow(data) / 100, sufix = "") {
  s_data <- data[sample(1:nrow(data), nrow(data)/100), ]
  mean_ratio <- by(s_data$ratio, s_data$code.well, mean)
  s_data$mean_ratio <- mean_ratio[s_data$code.well]
  p <- ggplot(data = s_data,
        aes(x = code.well, y = ratio, color = drug_status)) +
    geom_violin(alpha = 0) +
    geom_point(data = s_data,
              aes(x = code.well, y = mean_ratio, color = drug_status)) +
    facet_wrap(~code.well, scale = "free_x", ncol = 12) +
    theme_bw()
  print(p)
  outdir <- mk_outdir(fcs_data, "summary")
  ggplot2::ggsave(
    filename = paste0(outdir, "well_ratio", sufix, ".pdf"), plot = p,
    width = 29.7, height = 21, units = "cm", scale = 2
  )
  return(p)
}

#' plot the ratio distribution per column
#'
#' @param data data.frame object
#' @param sample_size (default: nrow(data) / 100) size of the subsample of data
#' to draw, to speedup the ploting proccess
#' @param file pdf file the save the object
#' @return a ggplot2 object
#' @examples
#' \dontrun{
#' plot_column_ratio(fsc_data)
#' }
#' @import ggplot2
#' @export plot_well_ratio
plot_column_ratio <- function(data, sample_size = nrow(data) / 100,
                              sufix = "") {
  s_data <- data[sample(1:nrow(data), nrow(data)/100), ]
  mean_ratio <- by(s_data$ratio, s_data$column, mean)
  s_data$mean_ratio <- mean_ratio[s_data$column]
  p <- ggplot(data = s_data,
        aes(x = column, y = ratio)) +
    geom_violin(alpha = 0) +
    geom_point(data = s_data,
              aes(x = column, y = mean_ratio)) +
    facet_wrap(~column, scale = "free_x", ncol = 12) +
    theme_bw()
  print(p)
  outdir <- mk_outdir(fcs_data, "summary")
  ggplot2::ggsave(
    filename = paste0(outdir, "column_ratio", sufix, ".pdf"), plot = p,
    width = 29.7, height = 11, units = "cm", scale = 2
  )
  return(p)
}

#' plot the ratio distribution per line
#'
#' @param data data.frame object
#' @param sample_size (default: nrow(data) / 100) size of the subsample of data
#' to draw, to speedup the ploting proccess
#' @param file pdf file the save the object
#' @return a ggplot2 object
#' @examples
#' \dontrun{
#' plot_line_ratio(fsc_data)
#' }
#' @import ggplot2
#' @export plot_well_ratio
plot_line_ratio <- function(data, sample_size = nrow(data) / 100, sufix = "") {
  s_data <- data[sample(1:nrow(data), nrow(data)/100), ]
  mean_ratio <- by(s_data$ratio, s_data$line, mean)
  s_data$mean_ratio <- mean_ratio[s_data$line]
  p <- ggplot(data = s_data,
        aes(x = line, y = ratio)) +
    geom_violin(alpha = 0) +
    geom_point(data = s_data,
              aes(x = line, y = mean_ratio)) +
    facet_wrap(~line, scale = "free_x", ncol = 8) +
    theme_bw()
  print(p)
  outdir <- mk_outdir(fcs_data, "summary")
  ggplot2::ggsave(
    filename = paste0(outdir, "line_ratio", sufix, ".pdf"), plot = p,
    width = 29.7, height = 11, units = "cm", scale = 2
  )
  return(p)
}
