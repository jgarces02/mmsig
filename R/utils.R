#' @export
rotatedAxisElementText = function(angle,position='x'){
  angle     = angle[1];
  position  = position[1]
  positions = list(x=0,y=90,top=180,right=270)
  if(!position %in% names(positions))
    stop(sprintf("'position' must be one of [%s]",paste(names(positions),collapse=", ")),call.=FALSE)
  if(!is.numeric(angle))
    stop("'angle' must be numeric",call.=FALSE)
  rads  = (angle - positions[[ position ]])*pi/180
  hjust = 0.5*(1 - sin(rads))
  vjust = 0.5*(1 + cos(rads))
  element_text(angle=angle,vjust=vjust,hjust=hjust)
}

#' @export
scale_fill_sigs <- function(...){
  ggplot2:::manual_scale('fill',
                         values = setNames(c(RColorBrewer::brewer.pal(8, "Dark2"), "lightblue", "#FB8072"),
                                           c("SBS1", "SBS2", "SBS5", "SBS8", "SBS9", "SBS13", "SBS18", "SBS-MM1", "SBS35", "SBS84")),
                         ...)
}

#' @export
vcfCleaner <- function(snvs){

  snvs <- snvs[snvs$ref %in% c("A","C","G","T") & snvs$alt %in% c("A","C","G","T"),] # remove indels
  snvs$chr <- as.character(snvs$chr)
  snvs <- snvs[snvs$chr %in% paste0("chr", c(1:22, "X")),]
  snvs$pos <- as.numeric(snvs$pos)
  snvs <- snvs[stats::complete.cases(snvs),]
  return(snvs)
}
