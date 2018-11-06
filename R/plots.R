#' plot the ratio distribution per well
#'
#' @param data data.frame object
#' @param sample_size (default: nrow(data) / 100) size of the subsample of data
#' to draw, to speedup the ploting proccess
#' @param sufix (default: "") sufix to append to pdf name
#' @return a ggplot2 object
#' @examples
#' \dontrun{
#' plot_plates(fsc_data)
#' }
#' @importFrom ggplot2 ggplot aes geom_violin geom_point facet_wrap labs theme_bw scale_fill_gradient
#' @export plot_well
plot_well <- function(data, sample_size = nrow(data) / 100, sufix = "") {
  s_data <- data[sample(1:nrow(data), sample_size), ]
  for (x in c("ratio", "Y1.A", "B1.A")) {
    s_data$x <- s_data[[x]]
    mean_x <- by(s_data$x, s_data$code.well, mean)
    s_data$mean_x <- mean_x[s_data$code.well]
    p <- ggplot2::ggplot()
    if ("scaled_pval" %in% colnames(s_data)) {
      s_data$scaled_pval <- as.vector(s_data$scaled_pval)
      p <- p + ggplot2::geom_violin(data = s_data,
          ggplot2::aes(x = code.well, y = x, fill = scaled_pval)) +
          ggplot2::scale_fill_gradientn(colors = heat.colors(100),
                                        na.value = "grey50",
                                        breaks = c(0, 1),
                                        labels = c("low p-value",
                                                   "high p-value"),
                                        guide = "colourbar")
    } else {
      p <- p + ggplot2::geom_violin(data = s_data,
          ggplot2::aes(x = code.well, y = x, fill = drug_status))
    }
    p <- p + ggplot2::geom_point(data = s_data,
              ggplot2::aes(x = code.well, y = mean_x)) +
      ggplot2::facet_wrap(~code.well, scale = "free_x", ncol = 12) +
      ggplot2::labs(y = x) +
      ggplot2::theme_bw()
    outdir <- mk_outdir(data, "summary")
    ggplot2::ggsave(
      filename = paste0(outdir, "well_", x, sufix, ".pdf"), plot = p,
      width = 29.7, height = 21, units = "cm", scale = 2
    )
  }
}

#' plot the ratio distribution per well
#'
#' @param data data.frame object with a pval column
#' @return data data.frame object with a scaled_pval column
#' @examples
#' \dontrun{
#' data <- scaled_pval(data)
#' }
#' @importFrom scales rescale
#' @export scaled_pval
scaled_pval <- function(data) {
  data$scaled_pval <- data$pval
  data$scaled_pval[!is.na(data$pval) & data$scaled_pval < 10^-10] <- 10^-10
  data$scaled_pval[!is.na(data$pval)]<- log10(
    data$scaled_pval[!is.na(data$pval)]
  )
  print(summary(data$scaled_pval))
  data$scaled_pval[!is.na(data$pval)] <- scales::rescale(
    data$scaled_pval[!is.na(data$pval)]
  )
  print(summary(data$scaled_pval))
  return(data)
}

#' plot the ratio distribution per column
#'
#' @param data data.frame object
#' @param sample_size (default: nrow(data) / 100) size of the subsample of data
#' to draw, to speedup the ploting proccess
#' @param sufix (default: "") sufix to append to pdf name
#' @return a ggplot2 object
#' @examples
#' \dontrun{
#' plot_column(fsc_data)
#' }
#' @importFrom ggplot2 ggplot aes geom_violin geom_point facet_wrap labs theme_bw
#' @export plot_column
plot_column <- function(data, sample_size = nrow(data) / 100,
                              sufix = "") {
  s_data <- data[sample(1:nrow(data), nrow(data)/100), ]
  for (x in c("ratio", "Y1.A", "B1.A")) {
    s_data$x <- s_data[[x]]
    mean_x <- by(s_data$x, s_data$column, mean)
    s_data$mean_x <- mean_x[s_data$column]
    p <- ggplot2::ggplot(data = s_data,
          ggplot2::aes(x = column, y = x)) +
      ggplot2::geom_violin(alpha = 0) +
      ggplot2::geom_point(data = s_data,
                ggplot2::aes(x = column, y = mean_x)) +
      ggplot2::facet_wrap(~column, scale = "free_x", ncol = 12) +
      ggplot2::labs(y = x) +
      ggplot2::theme_bw()
    outdir <- mk_outdir(data, "summary")
    ggplot2::ggsave(
      filename = paste0(outdir, "column_", x, sufix, ".pdf"), plot = p,
      width = 29.7, height = 11, units = "cm", scale = 2
    )
  }
}

#' plot the ratio distribution per line
#'
#' @param data data.frame object
#' @param sample_size (default: nrow(data) / 100) size of the subsample of data
#' to draw, to speedup the ploting proccess
#' @param sufix (default: "") sufix to append to pdf name
#' @return a ggplot2 object
#' @examples
#' \dontrun{
#' plot_line(fsc_data)
#' }
#' @importFrom ggplot2 ggplot aes geom_violin geom_point facet_wrap labs theme_bw
#' @export plot_line
plot_line <- function(data, sample_size = nrow(data) / 100, sufix = "") {
  s_data <- data[sample(1:nrow(data), nrow(data)/100), ]
  for (x in c("ratio", "Y1.A", "B1.A")) {
    s_data$x <- s_data[[x]]
    mean_x <- by(s_data$x, s_data$line, mean)
    s_data$mean_x <- mean_x[s_data$line]
    p <- ggplot2::ggplot(data = s_data,
          ggplot2::aes(x = line, y = x)) +
      ggplot2::geom_violin(alpha = 0) +
      ggplot2::geom_point(data = s_data,
                ggplot2::aes(x = line, y = mean_x)) +
      ggplot2::facet_wrap(~line, scale = "free_x", ncol = 8) +
      ggplot2::labs(y = x)
      ggplot2::theme_bw()
    outdir <- mk_outdir(data, "summary")
    ggplot2::ggsave(
      filename = paste0(outdir, "line_", x, sufix, ".pdf"), plot = p,
      width = 29.7, height = 11, units = "cm", scale = 2
    )
  }
}
