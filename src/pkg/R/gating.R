#' remove debris from data
#'
#' @param fcs_data an object of class flowSet
#' @return an object of class flowSet
#' @examples
#' \dontrun{
#' fsc_data <- rm_debris(fsc_data)
#' }
#' @importFrom flowClust flowClust getEstimates split
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
  par(mfrow = c(16, 6))
  for (i in 1:length(x)) {
    res1 <- flowClust::flowClust(
      fcs_data[[i]],
      varNames=c("FSC.A", "SSC.A"),
      K=2,
      B=100
    )
    plot(res1, data=fcs_data[[i]], level=0.95, z.cutoff=0)
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
#' @importFrom flowWorkspace GatingSet asinhtGml2_trans transformerList
#' @importFrom flowCore transform
#' @importFrom openCyto add_pop
#' @importFrom ggcyto ggcyto
#' @export rm_nonsinglets
rm_nonsinglets <- function(fcs_data) {
  outdir <- mk_outdir(fcs_data, "gating")
  wf <- flowWorkspace::GatingSet(fcs_data)
  asinhTrans <- flowWorkspace::asinhtGml2_trans()
  tl <- flowWorkspace::transformerList(c("Y1.A", "B1.A"), asinhTrans)
  wf <- flowCore::transform(wf, tl)
  openCyto::add_pop(
    wf, alias = "singlets", pop = "singlets", parent = "root",
    dims = "FSC.A,FSC.H", gating_method = "singletGate",
    gating_args = "wider_gate=TRUE"
  )
  p <- ggcyto::ggcyto(wf, mapping = aes(x = FSC.A,y = FSC.H)) +
    geom_hex(bins = 50) +
    geom_gate("singlets") +
    ggcyto_par_set(limits = "instrument") +
    labs(title = "Singlets gate")
  ggsave(
    filename = paste0(outdir, "singlets.pdf"), plot = p,
    width = 29.7, height = 21, units = "cm", scale = 2
  )
  # gating fluo
  fcs_singlets <- getData(wf, "singlets")
  return(fcs_singlets)
}

#' remove outliers on fluorecence channels
#'
#' @param fcs_data an object of class flowSet
#' @return an object of class flowSet
#' @examples
#' \dontrun{
#' fsc_data <- rm_nonfluop(fsc_data)
#' }
#' @importFrom flowClust flowClust getEstimates split
#' @export rm_nonfluop
rm_nonfluop<- function(fcs_data) {
  outdir <- mk_outdir(fcs_data, "gating")
  fcs_fluo <- fcs_data
  pdf(paste0(outdir, "fluo.pdf"), width = 80.3, height = 110.7)
  par(mfrow = c(16, 6))
  for (i in 1:length(x)) {
    res1 <- flowClust::flowClust(fcs_data[[i]], varNames=c("Y1.A", "B1.A"), K=1, B=100)
    plot(res1, data=fcs_data[[i]], level=0.85, z.cutoff=0)
    fcs_fluo[[i]] <- split(
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


