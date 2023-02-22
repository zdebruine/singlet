#' @title (Shared) Local Nearest-neighbor graph construction
#' 
#' @description Computes the \code{k.param} nearest neighbors within a spatial radius for a given dataset. Can also optionally (via \code{compute.SNN}), construct a shared nearest neighbor graph by calculating the neighborhood overlap (Jaccard index) between every cell and it's \code{k.param} nearest neighbors. Local KNN calculations are exact.
#' 
#' @details
#' IMPORTANT: You must make sure that your \code{radius} is given in the same units as your 
#' \code{spatial.reduction} coordinates, and that your \code{spatial.reduction} gives fixed 
#' coordinates. This means distance on x-coordinates must be equal to distance on y-coordinates.
#'  Many spatial assays store distances in x and y scaled between 0 and 1, which is NOT going 
#'  to work. You must use \code{\link{RescaleSpatial}} to convert back to a fixed coordinate 
#'  system. If your radius is 5, this function will look for neighbors within a distance of 5 
#'  from a given point as determined by your spatial coordinates.
#'  
#' @param object An object
#' @param ... not implemented
#' @param k.param Defines k for the k-nearest neighbor algorithm
#' @param spatial.reduction Spatial coordinates to use as input for building the (S)NN. Ensure that radius is given in the same units as spatial coordinates, and that spatial coordinates are fixed on both axes (not scaled).
#' @param compute.SNN also compute the shared nearest neighbor graph
#' @param prune.SNN Sets the cutoff for acceptable Jaccard index when computing the neighborhood overlap for the SNN construction. Any edges with values less than or equal to this will be set to 0 and removed from the SNN graph. Essentially sets the stringency of pruning (0 = no pruning, 1 = prune everything).
#' @param prune.KNN Sets the cutoff for acceptable distance when computing the neighborhood for the Local KNN graph construction. Any edges with values less than or equal to this will be set to 0 and removed from the KNN graph.  Essentially sets the stringency of pruning (0 = no pruning, 1 = prune everything when distance is "jaccard" or "cosine", otherwise whatever the equivalent is in the distance specified).
#' @param verbose print output to the console
#' @param return.dist return distances to nearest neighbors rather than a binary result
#' @param nn.metric Distance metric for nearest neighbors search. Options include: jaccard, cosine, euclidean, manhattan, hamming, and kl (kullback-leibler divergence).
#' @param use.dist use distance instead of similarity (i.e. find k-furthest-neighbors). Useful for edge detection. Applies only to \code{metric = c("jaccard", "cosine")}.
#' @param reduction Reduction to use as input for building the (S)NN
#' @param dims Dimensions of the reduction to use as input (\code{NULL} = use all dimensions in reduction)
#' @param graph.name Naming parameter for stored (S)NN graph. Default is \code{<reduction>_local_(s)nn}. To store both the neighbor graph and the shared nearest neighbor graph, you must supply a vector containing two names to the \code{graph.name} parameter. The first element in the vector will be used to store the nearest neighbor graph, and the second element will be used to store the shared nearest neighbor graph. If only one name is supplied, only the nearest neighbor graph is stored.
#' @param threads number of threads to use for parallelization
#' @returns an object (Seurat object with graph, or just a graph)
#' @export
#' @rdname FindLocalNeighbors
#' @aliases FindLocalNeighbors.Seurat
#'
FindLocalNeighbors.Seurat <- function(
    object, 
    k.param = 20,
    spatial.radius = 4,
    spatial.reduction = "spatial",
    reduction = "nmf",
    nn.metric = "jaccard",
    use.dist = FALSE,
    compute.SNN = TRUE,
    prune.SNN = 1/15, 
    prune.KNN = 1/10,
    return.dist = FALSE,
    verbose = FALSE, 
    dims = NULL, 
    graph.name = NULL, 
    threads = 0,
    ...){

    if(!(spatial.reduction %in% names(object@reductions))) 
      stop(spatial.reduction, "was not found in 'object@reductions'")
  
    if(!(reduction %in% names(object@reductions)))
      stop(reduction, "was not found in 'object@reductions'")
  
    sp <- object@reductions[[spatial.reduction]]@cell.embeddings
    
    if(max(sp[,1]) == 1 & max(sp[,2] == 1))
      warning("maximum value in spatial.reduction is 1 in both dimensions. Double-check that your coordinates are fixed and on the same scale in both dimensions.")
    
    h <- t(object@reductions[[reduction]]@cell.embeddings)
    
    if(ncol(h) != nrow(sp))
      stop("there were", ncol(h), "samples in reduction", reduction,"but", nrow(sp), "samples in spatial reduction", spatial.reduction)

    if(!is.null(dims)){
      if(max(dims) > nrow(h)) stop("you requested up to", max(dims), "dims but there are only", nrow(h), "dimensions in reduction", reduction)
      h <- h[dims, ]
    }
    
    if(spatial.radius > max(sp[,1]) | spatial.radius > max(sp[,2])){
      stop("your spatial radius of", spatial.radius, "is greater than the maximum value of", max(sp[,1], "or", max(sp[,2]),"in your spatial coordinates reduction", spatial.reduction))
    }
    
    if(!(nn.metric %in% c("jaccard", "euclidean", "manhattan", "hamming", "kl", "cosine"))) 
      stop("specified nn.metric =", nn.metric, "is not one of c('jaccard', 'euclidean', 'manhattan', 'hamming', 'kl', or 'cosine')")

    if(use.dist & !(nn.metric %in% c("jaccard", "cosine")))
      stop("it doesn't make sense to use dissimilarity (use.dist = FALSE) on a distance metric not strictly bounded between 0 and 1. Try using 'cosine' or 'jaccard' distance instead.")

    if(prune.KNN >= 1) stop("prune.knn must be less than 1. You are currently asking to prune everything. You should frankly feel rather stupid.")
    if(prune.SNN >= 1) stop("prune.snn must be less than 1. You are currently asking to prune everything. You should frankly feel rather stupid.")

    if(compute.SNN && !is.null(graph.name) && length(graph.name) != 1) 
      stop("compute.SNN = TRUE. If you read the docs you would know you need to provide two a 'graph.name' vector with two names:  one for the KNN, and one for the SNN.")

    if(is.null(graph.name)){
      if(compute.SNN) graph.name <- c(paste0(reduction, "_local_nn"), paste0(reduction, "_local_snn"))
      else graph.name <- c(paste0(reduction, "_local_nn"))
    }
    
    # call to C++ routine    
    object@graphs[[graph.name[[1]]]] <- c_LKNN(h, sp[,1], sp[,2], k.param, spatial.radius, nn.metric, !use.dist, prune.KNN, verbose, threads)
    if(!return.dist) object@graphs[[graph.name[[1]]]]@x <- rep(1, length(object@graphs[[graph.name[[1]]]]@x))
    
    if(compute.SNN) object@graphs[[graph.name[[2]]]] <- c_SNN(object@graphs[[graph.name[[1]]]], prune.SNN, threads)
    
    object
}


#' @rdname FindLocalNeighbors
#' @name FindLocalNeighbors
#' @export
#'
FindLocalNeighbors <- function(object, ...) {
  UseMethod("FindLocalNeighbors")
}

#' @rdname FindLocalNeighbors
#' @name FindLocalNeighbors
#' @export
#'
.S3method("FindLocalNeighbors", "Seurat", FindLocalNeighbors.Seurat)
