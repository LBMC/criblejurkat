#' perform all the analysis for a CribleJurkat analysis
#'
#' @param data_path path to the data folder
#' @param results_path path to the results folder
#' @return TRUE if everythings ran correclty
#' @examples
#' \dontrun{
#' CribleJurkat("data/2018_01_08_20171212/", 'results/2018_01_08_20171212/')
#' }
#' @export load_data
CribleJurkat <- function(data_path, results_path) {
  infos <- load_infos(data_path)
  data <- load_data(data_path)
}
