#' remove debris from data
#'
#' @param fcs_data an object of class flowSet
#' @return an object of class flowSet
#' @examples
#' \dontrun{
#' fsc_data <- rm_debris(fsc_data)
#' }
#' @importFrom flowClust flowClust getEstimates split plot
#' @export rm_debris
rm_debris <- function(fcs_data) {
  fcs_nonDebris <- fcs_data
  outdir <- mk_outdir(fcs_data, "gating")
  pdf(
    paste0(
      outdir,
      "nonDebris.pdf"
    ),
    width = 80.3, height = 110.7
  )
  par(mfrow = c(8, 12))
  for (i in 1:length(fcs_data)) {
    res1 <- flowClust::flowClust(
      fcs_data[[i]],
      varNames=c("FSC.A", "SSC.A"),
      K=2,
      B=100
    )
    flowClust::plot(res1, data=fcs_data[[i]], level=0.95, z.cutoff=0)
    cluster_location <- flowClust::getEstimates(res1)$locations
    nonDebris <- which(
      cluster_location[,1] == max(cluster_location[,1])
    )
    Debris <- which(
      cluster_location[,1] == min(cluster_location[,1])
    )
    abline(v = cluster_location[nonDebris,1],
           h = cluster_location[nonDebris,2])
    fcs_nonDebris[[i]] <- flowClust::split(
      fcs_data[[i]],
      res1,
      population = list(
        nonDebris = nonDebris,
        Debris = Debris
      )
    )$nonDebris
  }
  dev.off()
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
#' @importFrom flowWorkspace GatingSet transformerList flowJo_biexp_trans
#' @importFrom flowCore transform
#' @importFrom openCyto add_pop
#' @importFrom ggcyto ggcyto geom_gate ggcyto_par_set
#' @importFrom ggplot2 geom_hex labs ggsave aes facet_wrap
#' @export rm_nonsinglets
rm_nonsinglets <- function(fcs_data) {
  outdir <- mk_outdir(fcs_data, "gating")
  wf <- flowWorkspace::GatingSet(fcs_data)
  bi_expTrans <- flowWorkspace::flowJo_biexp_trans()
  tl <- flowWorkspace::transformerList(c("Y1.A", "B1.A"), bi_expTrans)
  wf <- flowCore::transform(wf, tl)
  openCyto::add_pop(
    wf, alias = "singlets", pop = "+", parent = "root",
    dims = "FSC.A,FSC.H", gating_method = "singletGate",
    gating_args = "wider_gate=TRUE, maxit = 100"
  )
  p <- ggcyto::ggcyto(wf, mapping = ggplot2::aes(x = FSC.A,y = FSC.H)) +
    ggplot2::geom_hex(bins = 50) +
    ggcyto::geom_gate("singlets") +
    ggcyto::ggcyto_par_set(limits = "instrument") +
    ggplot2::labs(title = "Singlets gate") +
    ggplot2::facet_wrap(~well, ncol = 12)
  ggplot2::ggsave(
    filename = paste0(outdir, "singlets.pdf"), plot = p,
    width = 29.7, height = 21, units = "cm", scale = 2
  )
  fcs_singlets <- getData(wf, "singlets")
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
#' @importFrom flowClust flowClust getEstimates split
#' @export rm_nonfluo
rm_nonfluo <- function(fcs_data) {
  outdir <- mk_outdir(fcs_data, "gating")
  fcs_fluo <- fcs_data
  pdf(paste0(outdir, "fluo.pdf"), width = 80.3, height = 110.7)
  par(mfrow = c(12, 8))
  for (i in 1:length(fcs_data)) {
    res1 <- flowClust::flowClust(fcs_data[[i]], varNames=c("Y1.A", "B1.A"), K=1, B=100)
    plot(res1, data=fcs_data[[i]], level=0.85, z.cutoff=0)
    fcs_fluo[[i]] <- flowClust::split(
      fcs_data[[i]],
      res1,
      population = list(
        fluo = 1
      )
    )$fluo
  }
  dev.off()
  return(fcs_fluo)
}


