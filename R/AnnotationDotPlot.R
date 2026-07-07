

plotAnnotation = function(dataset,datafield,transparency=NULL){
  assertthat::assert_that(!transparency %in% c("FDR","logodds","both") & !is.numeric(transparency) & !is.null(transparency),msg = "Invalid transparency values. Use 'FDR','logodds','both', or enter custom values.")
  assertthat::assert_that(length(transparency != nrow(dataset)),msg = "Length of transparency not equal number of rows in dataset")
  
  transparency_FDR = rev(normalize(dataset$negative_log10_fdr))
  transparency_logodds = normalize(dataset$lods)
  transparency_mean = (transparency_FDR+transparency_logodds)/2
  
  #how to do switch statement in R? :(
  if (transparency=="both"){
    transparency = normalize(transparency_mean)
  }else if (transparency=="FDR"){
    transparency = transparency_FDR
  }else if (transparency=="logodds"){
    transparency = transparency_logodds
  }else{
    transparency = 1
  }
  
  p <- ggplot(dataset, 
                      aes(x = factor, 
                          y = field, 
                          color = negative_log10_fdr, 
                          size = lods)) + 
    geom_point(alpha = transparency) + 
    scale_color_viridis_c(direction = -1,option = "B", end = 0.9) +
    theme_minimal() + 
    labs(y = datafield, 
         x = "NMF factor", 
         color = "FDR\n(-log10)", 
         size = "Association\n(log-odds)") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    guides(alpha = "none") + 
    NULL
  return(p)
}

normalize <- function(x, na.rm = TRUE) {
  return((x- min(x)) /(max(x)-min(x)))
}
