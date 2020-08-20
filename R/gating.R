#' remove debris from data
#'
#' @param fcs_data an object of class flowSet
#' @return an object of class flowSet
#' @examples
#' \dontrun{
#' fsc_data <- rm_debris(fsc_data)
#' }
#' @importFrom flowClust flowClust getEstimates split plot
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics abline par
#' @importFrom utils txtProgressBar capture.output
#' @export rm_debris
rm_debris <- function(fcs_data) {
  message("removing debris")
  fcs_nonDebris <- fcs_data
  outdir <- mk_outdir(fcs_data, "gating")
  grDevices::pdf(
    paste0(
      outdir,
      "nonDebris.pdf"
    ),
    width = 29.7, height = 21
  )
  graphics::par(mfrow = c(8, 12))
  pb <- utils::txtProgressBar(min = 0, max = length(fcs_data),
                             initial = 1, style = 3)
  for (i in 1:length(fcs_data)) {
    suppressMessages(
      res1 <- flowClust::flowClust(
        fcs_data[[i]],
        varNames = c("FSC.A", "SSC.A"),
        K = 1,
        B = 1000,
        level = 0.90,
        z.cutoff = 0
      )
    )
    cluster_location <- flowClust::getEstimates(res1)$locations
    cluster_prop <- flowClust::getEstimates(res1)$proportions

    nonDebris <- which(
      cluster_location[,1] == max(cluster_location[,1])
    )
    cluster_pop <- list(
      nonDebris = nonDebris
    )

    utils::capture.output(
      flowClust::plot(res1, data = fcs_data[[i]],
                      level = 0.90, z.cutoff = 0
      )
    )
    graphics::abline(v = cluster_location[nonDebris,1],
           h = cluster_location[nonDebris,2])
    fcs_nonDebris[[i]] <- flowClust::split(
      fcs_data[[i]],
      res1,
      population = cluster_pop,
      rm.outliers = TRUE,
    )$nonDebris
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)
  grDevices::dev.off()
  return(fcs_nonDebris)
}

#' remove non singlets
#'
#' @param fcs_data an object of class flowSet
#' @return an object of class flowSet
#' @examples
#' \dontrun{
#' fsc_data <- rm_nonsinglets(fsc_data)
#' }
#' @importFrom flowClust flowClust getEstimates split
#' @importFrom flowWorkspace GatingSet transformerList flowJo_biexp_trans getData
#' @importFrom flowCore transform
#' @importFrom openCyto add_pop
#' @importFrom ggcyto ggcyto geom_gate ggcyto_par_set
#' @importFrom ggplot2 geom_hex labs ggsave aes facet_wrap
#' @export rm_nonsinglets
rm_nonsinglets <- function(fcs_data) {
  message("removing non-singlets")
  outdir <- mk_outdir(fcs_data, "gating")
  bi_expTrans <- flowWorkspace::flowjo_biexp_trans()
  tl <- flowCore::transformList(
    c("Y1.A", "B1.A"), list(bi_expTrans$transform, bi_expTrans$transform)
  )
  fcs_singlets <- fcs_data
  flowWorkspace::sampleNames(fcs_singlets) <- gsub(
    ".*/([^/]*)", "\\1", flowWorkspace::sampleNames(fcs_singlets)
  )
  fcs_singlets <- flowWorkspace::flowSet_to_cytoset(fcs_singlets)
  fcs_singlets <- flowCore::transform(fcs_singlets, tl)
  fcs_singlets <- flowWorkspace::GatingSet(fcs_singlets)
  openCyto::gs_add_gating_method(
    fcs_singlets, alias = "singlets", pop = "+", parent = "root",
    dims = "FSC.A,FSC.H", gating_method = "singletGate",
    gating_args = "wider_gate=TRUE, maxit = 100"
  )
  p <- ggcyto::ggcyto(fcs_singlets,
                      mapping = ggplot2::aes(x = FSC.A,y = FSC.H)) +
    ggplot2::geom_hex(bins = 50) +
    ggcyto::geom_gate("singlets") +
    ggcyto::ggcyto_par_set(limits = "instrument") +
    ggplot2::labs(title = "Singlets gate") +
    ggplot2::facet_wrap(~well, ncol = 12)
  ggplot2::ggsave(
    filename = paste0(outdir, "singlets.pdf"), plot = p,
    width = 29.7, height = 21, units = "cm", scale = 2
  )
  fcs_singlets <- flowWorkspace::gs_pop_get_data(fcs_singlets, "singlets")
  fcs_singlets <- flowWorkspace::cytoset_to_flowSet(fcs_singlets)
  flowWorkspace::sampleNames(fcs_singlets) <- flowWorkspace::sampleNames(
    fcs_data
  )
  return(fcs_singlets)
}

#' remove outliers on fluorecence channels
#'
#' @param fcs_data an object of class flowSet
#' @return an object of class flowSet
#' @examples
#' \dontrun{
#' fsc_data <- rm_nonfluo(fsc_data)
#' }
#' @importFrom flowClust flowClust getEstimates split plot
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics par
#' @importFrom utils txtProgressBar capture.output
#' @export rm_nonfluo
rm_nonfluo <- function(fcs_data) {
  message("removing outliers")
  outdir <- mk_outdir(fcs_data, "gating")
  fcs_fluo <- fcs_data[seq(along = fcs_data)]
  grDevices::pdf(paste0(outdir, "fluo.pdf"), width = 29,7, height = 21)
  graphics::par(mfrow = c(8, 12))
  pb <- utils::txtProgressBar(min = 0, max = length(fcs_data),
                             initial = 1, style = 3)
  for (i in 1:length(fcs_data)) {
    res <- suppressMessages(
      flowClust::flowClust(
        fcs_data[[i]],
        varNames = c("Y1.A", "B1.A"),
        K = 1,
        B = 100
      )
    )
    utils::capture.output(
      flowClust::plot(
        res,
        data = fcs_data[[i]],
        level = 0.85,
        z.cutoff = 0
      )
    )
    res <- flowClust::split(
      fcs_data[[i]],
      res,
      rm.outliers = TRUE,
      population = list(
        fluo = 1
      )
    )
    fcs_fluo[[i]] <- res$fluo
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)
  grDevices::dev.off()
  return(fcs_fluo)
}

