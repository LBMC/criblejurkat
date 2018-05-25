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
#' @import ggplot2
#' @export plot_well
plot_well <- function(data, sample_size = nrow(data) / 100, sufix = "") {
  s_data <- data[sample(1:nrow(data), nrow(data)/100), ]
  for (x in c("ratio_norm", "Y1.A", "B1.A")) {
    s_data$x <- s_data[[x]]
    mean_x <- by(s_data$x, s_data$code.well, mean)
    s_data$mean_x <- mean_x[s_data$code.well]
    p <- ggplot()
    if ("signif" %in% colnames(s_data)) {
      p <- p + geom_violin(alpha = 0, data = s_data,
          aes(x = code.well, y = x, color = signif))
    } else {
      p <- p + geom_violin(alpha = 0, data = s_data,
          aes(x = code.well, y = x, color = drug_status))
    }
    p <- p + geom_point(data = s_data,
              aes(x = code.well, y = mean_x, color = drug_status)) +
      facet_wrap(~code.well, scale = "free_x", ncol = 12) +
      labs(y = x) +
      theme_bw()
    print(p)
    outdir <- mk_outdir(fcs_data, "summary")
    ggplot2::ggsave(
      filename = paste0(outdir, "well_", x, sufix, ".pdf"), plot = p,
      width = 29.7, height = 21, units = "cm", scale = 2
    )
  }
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
#' @export plot_column
plot_column <- function(data, sample_size = nrow(data) / 100,
                              sufix = "") {
  s_data <- data[sample(1:nrow(data), nrow(data)/100), ]
  for (x in c("ratio_norm", "Y1.A", "B1.A")) {
    s_data$x <- s_data[[x]]
    mean_x <- by(s_data$x, s_data$column, mean)
    s_data$mean_x <- mean_x[s_data$column]
    p <- ggplot(data = s_data,
          aes(x = column, y = x)) +
      geom_violin(alpha = 0) +
      geom_point(data = s_data,
                aes(x = column, y = mean_x)) +
      facet_wrap(~column, scale = "free_x", ncol = 12) +
      labs(y = x) +
      theme_bw()
    print(p)
    outdir <- mk_outdir(fcs_data, "summary")
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
#' @param file pdf file the save the object
#' @return a ggplot2 object
#' @examples
#' \dontrun{
#' plot_line_ratio(fsc_data)
#' }
#' @import ggplot2
#' @export plot_line
plot_line <- function(data, sample_size = nrow(data) / 100, sufix = "") {
  s_data <- data[sample(1:nrow(data), nrow(data)/100), ]
  for (x in c("ratio_norm", "Y1.A", "B1.A")) {
    s_data$x <- s_data[[x]]
    mean_x <- by(s_data$x, s_data$line, mean)
    s_data$mean_x <- mean_x[s_data$line]
    p <- ggplot(data = s_data,
          aes(x = line, y = x)) +
      geom_violin(alpha = 0) +
      geom_point(data = s_data,
                aes(x = line, y = mean_x)) +
      facet_wrap(~line, scale = "free_x", ncol = 8) +
      labs(y = x)
      theme_bw()
    print(p)
    outdir <- mk_outdir(fcs_data, "summary")
    ggplot2::ggsave(
      filename = paste0(outdir, "line_", x, sufix, ".pdf"), plot = p,
      width = 29.7, height = 11, units = "cm", scale = 2
    )
  }
}
