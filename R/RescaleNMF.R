#' @title Rescaling of cell and gene features
#'
#' @param object is a seurat object  
#' @param features are the list of cell or genes of interest.
#' @param reduction use to specifiy the reduction .
#' @param lambda is the threshold parameter between 0 and 1 which helps improve the relationship of v and either H or W.
#' @param reduction.name names in which the new rescaled seurate object is stored. 
#' @param ... other unknown parameters can be added to then function.
#' @export
#' @examples 
#' RescaleNMF(object=pbmc3k,features=c("AGRN","dvl1"),lambda=0.4)
setGeneric("RescaleNMF", function(object,...) {
  standardGeneric("RescaleNMF")
})
setMethod("RescaleNMF", signature("Seurat"), function(object,features=c("AGRN","DVL1"),reduction="nmf",lambda=0.5,reduction.name="nmf_rescaled",...) { 
  nmf_model <- object@reductions[[reduction]]
  is_w <- sum(features %in% rownames(nmf_model@feature.loadings))
  is_h <- sum(features %in% rownames(nmf_model@cell.embeddings))
  if(is_w > is_h){
    
    w <- nmf_model@feature.loadings
    h <- nmf_model@cell.embeddings
    #v <- colSums(abs(w[which(rownames(w) %in% features), ]))
    #v <- v / mean(v)
    
    wo <- colSums(abs(w[which(rownames(w) %in% features), ]))
    wi <- colSums(abs(w))
    d  <- wo / wi
    d  <- lambda*d + (1- lambda)* rep(1,length(d))
    d  <- (d - mean(d))+ 1
    
    #v <- lambda * v + (1 - lambda) * rep(1, length(v))
    nmf_model@feature.loadings <- as.matrix(w %*% Matrix::Diagonal(x=d))
    nmf_model@cell.embeddings <- as.matrix(h %*% Matrix::Diagonal(x=d))
    object@reductions[[reduction.name]] <- nmf_model
    object
    
  } else {
    
    w <- nmf_model@feature.loadings
    h <- nmf_model@cell.embeddings
    #v <- colSums(abs(h[which(rownames(h) %in% features), ]))
    #v <- v / mean(v)
    #v <- lambda * v + (1 - lambda) * rep(1, length(v))
    ho <- colSums(abs(h[which(rownames(h) %in% features), ]))
    hi <- colSums(abs(h))
    d  <- ho / hi
    d  <- lambda*d + (1- lambda)* rep(1,length(d))
    d  <- (d - mean(d))+ 1
    
    nmf_model@feature.loadings <- as.matrix(w %*% Matrix::Diagonal(x=d))
    nmf_model@cell.embeddings <- as.matrix(h %*% Matrix::Diagonal(x=d))
    object@reductions[[reduction.name]] <- nmf_model
    object 
    
    }
  })






