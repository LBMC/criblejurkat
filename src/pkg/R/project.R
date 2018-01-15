#' project data object
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export scdata
#' @format An \code{\link{R6Class}} generator object
#' @keywords infos counts
scdata <- R6::R6Class("scdata",
  private = list(
    genes = NULL, # a vector of genes names
    cells = NULL, # a vector of cells id
    counts = NULL, # a matrix of counts
    features = NULL, # a data.frame of cells features
    # private method to get position of cells and genes in rows and columns
    get_rc_counts = function(cells = NULL, genes = NULL, v = F) {
      c_num <- 1:self$getngenes
      if (!is.null(genes)) {
        c_num <- which(private$genes %in% genes)
      }
      r_num <- 1:self$getncells
      if (!is.null(cells)) {
        r_num <- which(private$cells %in% cells)
      }
      return(list(c_num = c_num, r_num = r_num))
    },
    # private method to get position of cells and features in rows and columns
    get_rc_features = function(cells = NULL, features = NULL, v = F) {
      c_num <- 1:self$getnfeatures
      if (!is.null(features)) {
        c_num <- which(colnames(private$features) %in% features)
      }
      r_num <- 1:self$getncells
      if (!is.null(cells)) {
        r_num <- which(private$cells %in% cells)
      }
      return(list(c_num = c_num, r_num = r_num))
    },
    # private method to get position of cells in common between features and
    # counts
    get_common_cells = function(features, counts, v = F) {
      id_features <- as.vector(features$id)
      id_counts <- rownames(counts)
      common_cells <- intersect(id_features, id_counts)
      if (v) {
        print(
          paste0(length(common_cells),
            " cells in common between infos and counts"))
        print("cells present in infos, but not counts")
        print(setdiff(id_features, id_counts))
        print("cells present in counts, but not infos")
        print(setdiff(id_counts, id_features))
      }
      r_features <- which(
        id_features %in% common_cells & !duplicated(id_features)
      )
      r_counts <- which(id_counts %in% common_cells & !duplicated(id_counts))
      return(list(features = r_features, counts = r_counts))
    },
    order_by_cells = function(v = F) {
      if (nrow(private$counts) == length(private$cells)) {
        private$counts <- private$counts[order(private$cells), ]
      } else {
          stop(paste0(
            "error: number of cells (",
            length(private$cells),
            ") don't match number of counts rows (",
            nrow(private$counts),
          ")"))
      }
      if (nrow(private$features) == length(private$cells)) {
        private$features <- private$features[order(self$getfeature("id")), ]
      } else {
          stop(paste0(
            "error: number of cells (",
            length(private$cells),
            ") don't match number of feature rows (",
            nrow(private$features),
          ")"))
      }
      private$cells <- rownames(private$counts)
    },
    transpose_count = function(counts, v = F) {
      counts <- as.matrix(counts)
      if (nrow(private$features) < nrow(counts)) {
        if (v) {
          print("transposing counts...")
        }
        return(t(counts))
      }
      return(counts)
    },
    set_na_to_zero = function(v = F) {
      if (v) {
        print("cells with NA's counts:")
        print(self$getcells[rowSums(is.na(self$getcounts)) != 0])
        print("genes with NA's counts:")
        print(self$getgenes[colSums(is.na(self$getcounts)) != 0])
      }
      private$counts[is.na(private$counts)] <- 0
    },
    display_dim = function(infos, counts, v = F) {
      if (v) {
        print(
          paste0("dim of infos :",
            nrow(infos),
            " x ",
            ncol(infos)
          )
        )
        print(
          paste0("dim of counts :",
            nrow(counts),
            " x ",
            ncol(counts)
          )
        )
      }
    }
    ),
  public = list(
    initialize = function(infos = NA, counts = NA, v = F) {
      private$features <- as.data.frame(infos)
      private$counts <- private$transpose_count(counts, v = v)
      private$genes <- colnames(private$counts)
      private$cells <- rownames(private$counts)
      private$display_dim(private$features, private$counts, v = v)
      if ("id" %in% colnames(private$features)) {
        if (v) {
          print("id column found in infos")
        }
        if (length(intersect(private$cells, as.vector(private$features$id)))
          == 0) {
          print(setdiff(private$cells, as.vector(private$features$id)))
          stop("error: id's in infos don't match cells name in counts")
        }
        rr_num <- private$get_common_cells(
          private$features,
          private$counts,
          v = v
        )
        private$features <- private$features[rr_num$features, ]
        private$counts <- private$counts[rr_num$counts, ]
        private$cells <- rownames(private$counts)
        private$display_dim(private$features, private$counts, v = v)
        private$order_by_cells(v = v)
        if (any(private$getfeature['id'] != private$cells) ){
          stop("error : features order don't match counts order")
        }
      } else {
        if (nrow(private$features) != nrow(private$counts)) {
          stop(
            paste0("error: number of cells (",
              nrow(private$features),
              ") differ between infos and counts (",
              nrow(private$counts)))
        }
      }
      private$set_na_to_zero(v = v)
      self$summary(v = v)
    },
    add = function(infos = NA, counts = NA, v = F) {
      features <- as.data.frame(infos)
      counts <- private$transpose_count(counts, v = v)
      genes <- colnames(counts)
      cells <- rownames(counts)
      private$display_dim(features, counts, v = v)
      if (length(intersect(private$cells, cells)) != 0) {
        stop("error: trying to add cells already present")
      }
      if (length(intersect(private$genes, genes)) != length(private$genes)) {
        print(private$cells)
        print(as.vector(private$features$id))
        stop("error: genes set don't match existing genes set")
      }
      if ("id" %in% colnames(features)) {
        if (v) {
          print("id column found in infos")
        }
        if (length(intersect(cells, as.vector(features$id))) == 0) {
          print(setdiff(private$cells, as.vector(private$features$id)))
          stop("error: id's in infos don't match cells name in counts")
        }
      } else {
        if (nrow(features) != nrow(counts)) {
          stop(
            paste0("error: number of cells (",
              nrow(private$features),
              ") differ between infos and counts (",
              nrow(private$counts)))
        }
      }
      rr_num <- private$get_common_cells(features, counts, v = v)
      features <- features[rr_num$features, ]
      counts <- counts[rr_num$counts, ]
      genes <- colnames(counts)
      cells <- rownames(counts)
      not_here <- !(cells %in% private$cells)
      private$counts <- rbind(private$counts, counts[not_here, ])
      private$features <- rbind(private$features, features[not_here, ])
      private$cells <- rownames(private$counts)
      private$order_by_cells(v = v)
      private$set_na_to_zero(v = v)
      if (any(private$getfeature['id'] != private$cells) ){
        stop("error : features order don't match counts order")
      }
    },
    summary = function(v = F) {
      if (v) {
        cat(paste0("ncol: ", ncol(private$counts), ".\n"))
        cat(paste0("nrow: ", nrow(private$counts), ".\n"))
        cat(paste0("features: ", ncol(private$features), ".\n"))
      }
    },
    # accessors methods
    getfeature = function(feature) {
      return(private$features[[feature]])
    },
    addfeature = function(feature) {
      private$features[[feature]] <- NA
    },
    setfeature = function(feature, value) {
      private$features[[feature]] <- value
    },
    getgene = function(gene) {
      c_num <- which(private$genes %in% gene)
      return(private$counts[, c_num])
    },
    getcountsw = function(cells = NULL, genes = NULL) {
      rc_num <- private$get_rc_counts(cells = cells, genes = genes)
      return(private$counts[rc_num$r_num, rc_num$c_num])
    },
    getfeaturesw = function(cells = NULL, features = NULL) {
      rc_num <- private$get_rc_features(cells = cells, features = features)
      return(private$features[rc_num$r_num, rc_num$c_num])
    },
    getcountso = function(cells = NULL, genes = NULL) {
      rc_num <- private$get_rc_counts(cells = cells, genes = genes)
      return(private$counts[-rc_num$r_num, -rc_num$c_num])
    },
    getfeatureso = function(cells = NULL, features = NULL) {
      rc_num <- private$get_rc_features(cells = cells, features = features)
      return(private$features[-rc_num$r_num, -rc_num$c_num])
    },
    copy = function(
      cells = NULL, genes = NULL, features = NULL, b_cells = NULL){
      if (is.null(b_cells)){
        b_cells <- T
      }
      if (is.null(cells) | length(cells) > 1){
        return(
          scdata$new(
            infos = self$getfeaturesw(
              cells = cells, features = features)[b_cells, ],
            counts = self$getcountsw(
              cells = cells, genes = genes)[b_cells, ]
          )
        )
      } else {
        return(
          scdata$new(
            infos = self$getfeaturesw(
              cells = cells, features = features),
            counts = self$getcountsw(
              cells = cells, genes = genes)
          )
        )
      }
    },
    select = function(
      cells = NULL, genes = NULL, features = NULL, b_cells = NULL) {
      return(
        self$copy(
          cells = cells,
          genes =  genes,
          features = features,
          b_cells = b_cells
      ))
    },
    order = function(cells = NULL, genes = NULL) {
      if (!is.null(cells)) {
        private$features <- private$features[cells, ]
        private$counts <- private$counts[cells, ]

      }
      if (!is.null(genes)) {
        private$counts <- private$counts[, genes]
      }
      private$genes <- colnames(private$counts)
      private$cells <- rownames(private$counts)
    },
    transform = function(FUN = function(x){x}) {
      private$counts <- FUN(private$counts)
    },
    update = function(counts = NULL, features = NULL) {
      if (!is.null(counts)) {
        private$counts <- counts
        private$genes <- colnames(private$counts)
        private$cells <- rownames(private$counts)
      }
      if (!is.null(features)) {
        private$features <- features
      }
    }
  ),
  active = list(
    # accessors methods without arguments
    getcounts = function() {
      return(private$counts)
    },
    getfeatures = function() {
      return(private$features)
    },
    getcells = function() {
      return(private$cells)
    },
    getgenes = function() {
      return(private$genes)
    },
    getncells = function() {
      return(length(private$cells))
    },
    getngenes = function() {
      return(length(private$genes))
    },
    getnfeatures = function() {
      return(ncol(private$features))
    }
  )
)

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
function CribleJurkat(data_path, results_path) {

}
