#' remove debris from data
#'
#' @param fcs_data an object of class flowSet
#' @return an object of class flowSet
#' @examples
#' \dontrun{
#' fsc_data <- rm_debris(fsc_data)
#' }
#' @importFrom flowClust flowClust getEstimates split
#' @export load_data
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
