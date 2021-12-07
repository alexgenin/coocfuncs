# 
# 
# This file contains functions that will compute and test coocurrence metrics 
#   on continous data.
# 

#'@export
cont_overlap <- function(df, tol = 0, 
                         groups = "transect", attribut = "attribut") { 
  
  # Create groups if absent
  if ( is.null(groups) || is.na(groups) ) { 
    groups <- "transect"
    df[ ,groups] <- 1
  }
  
  # Prepare attriboots 
  attribs <- as.factor(df[ ,attribut])
  attribs.names <- levels(attribs)
  attribs <- as.integer(attribs)-1 # adjust for indexing 
  
  # 
  N <- nrow(df) # Number of individuals
  Nsp <- length(unique(attribs)) # Number of species
  spcoocs <- cont_overlap_core(N, Nsp, df[ ,"xi"], df[ ,"xe"], attribs, 
                               as.factor(df[ ,groups]), tol = tol)
  colnames(spcoocs) <- rownames(spcoocs) <- attribs.names
  
  return(spcoocs)
}
