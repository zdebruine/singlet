#' Rescale spatial coordinates
#' 
#' Convert coordinates in the "spatial" reduction to natural numbers rather than values between 0 and 1. This allows for intuitive graph construction based on the radius surrounding any given cell (i.e. a radius of one corresponds to all cells next to the cell of interest)
#' 
#' @param object Seurat object
#' @param reduction the name of the spatial reduction to use
#' @export
#' @return Seurat object with rescaled spatial coordinates
#' @aliases RescaleSpatial.Seurat
#' @rdname RescaleSpatial
RescaleSpatial.Seurat <- function(object, reduction = "spatial"){
  df <- object@reductions[[reduction]]@cell.embeddings
  df[,1] <- df[,1] - min(df[,1])
  df[,2] <- df[,2] - min(df[,2])
  df[,1] <- df[,1] / max(df[,1])
  df[,2] <- df[,2] / max(df[,2])
  df[,1] <- df[,1] * 1 / median(diff(sort(unique(df[,1]))))
  df[,2] <- df[,2] * 1 / median(diff(sort(unique(df[,2]))))
  df <- round(df)
  object@reductions[[reduction]]@cell.embeddings <- df
  object
}


#' @rdname RunGCNMF
#' @name RunGCNMF
#' @export
#'
RescaleSpatial <- function(object, ...) {
  UseMethod("RescaleSpatial")
}

#' @rdname RescaleSpatial
#' @name RescaleSpatial
#' @export
#'
.S3method("RescaleSpatial", "Seurat", RescaleSpatial.Seurat)
