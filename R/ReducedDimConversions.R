#' Convert RcppML NMF to Seurat DimReduc
#'
#' @param nmf_model RcppML::nmf object
#' @param assay Assay associated with the reduction
#' @param reduction.key Prefix for latent dimensions
#'
#' @return Seurat DimReduc object
#'
#' @import Seurat
#' @export
nmfToDimReduc = function(nmf_model,assay,reduction.key){
  # colnames(nmf_model@w) <- gsub("nmf", reduction.key, colnames(nmf_model@w))
  # rownames(nmf_model@h) <- gsub("nmf", reduction.key, rownames(nmf_model@h))
  dims <- paste0(reduction.key, seq_len(ncol(nmf_model@w)))
  colnames(nmf_model@w) <- dims
  rownames(nmf_model@h) <- dims
  
  reduct <- new("DimReduc",
                cell.embeddings = t(nmf_model@h),
                feature.loadings = nmf_model@w,
                assay.used = assay,
                stdev = nmf_model@d,
                key = reduction.key,
                global = FALSE
  )
  
  return(reduct)
}


#' Convert LEM to RcppML NMF and places it into sce metadata
#'
#' @param sce SingleCellExperiment containing an LEM from fletchR results
#' @param LEM_slot Name of reducedDim containing the LEM ("iterativeLSI")
#' @param metadata_slot Name of metadata slot for storing the NMF;
#'                      NULL uses LEM_slot
#'
#' @return SingleCellExperiment with NMF stored in metadata
#'
#' @details Converts the LEM SVD or retrieves its NMF decomposition,
#'          preserving additional LEM metadata in the misc slot.
#'
#' @seealso irlbaToNmf
#' @import SingleCellExperiment
#' @import RcppML
#' @export
fletchRLEMToNMF = function(sce,
                           LEM_slot = "iterativeLSI",
                           metadata_slot = NULL) {
  if (is.null(metadata_slot)) {
    metadata_slot = LEM_slot
  }
  
  stopifnot(
    "metadata_slot already exists." =
      !(metadata_slot %in% names(metadata(sce)))
  )
  
  LEM = SingleCellExperiment::reducedDim(sce, LEM_slot)
  
  if ("svd" %in% names(LEM@metadata)) {
    
    reduct <- irlbaToNmf(LEM@metadata$svd)
    
  } else if ("nmf" %in% names(LEM@metadata)) {
    
    reduct <- LEM@metadata$nmf
    
  } else {
    stop("No 'svd' or 'nmf' found in LEM metadata.")
  }
  
  # Preserve other LEM metadata
  for (slot in setdiff(names(LEM@metadata), c("svd", "nmf"))) {
    reduct@misc[[slot]] <- LEM@metadata[[slot]]
  }
  
  reduct@misc[["OrigLEMObject"]] <- LEM_slot
  
  metadata(sce)[[metadata_slot]] <- reduct
  
  return(sce)
}

#' Convert irlba SVD to RcppML NMF
#'
#' @param svd SVD result containing u, d, and v
#'
#' @return RcppML::nmf object with w = u, d = d, and h = t(v)
#'
#' @details Wraps SVD factors in an NMF object without enforcing
#'          nonnegativity.
#'
#' @import RcppML
#' @export
irlbaToNmf = function(svd) {
  stopifnot(
    all(c("u", "d", "v") %in% names(svd)),
    ncol(svd$u) == length(svd$d),
    ncol(svd$v) == length(svd$d)
  )
  
  methods::new(
    "nmf",
    w = svd$u,
    d = svd$d,
    h = t(svd$v),
    misc = list()
  )
}