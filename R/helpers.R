# 
# 
# Left-join two dfs
# 

ljoin <- function(a, b, ...) { 
  plyr::join(a, b, type = "left", match = "first", ...)
}


tabularize <- function(coocmat, full = FALSE, trim = FALSE, 
                       no_inter_value = c(NA, 0), name = "ses") { 
  
  if ( is.null(dim(coocmat)) ) { 
    return(NA)
  }
  if ( is.null(colnames(coocmat)) ) { 
    spnames <- seq.int(ncol(coocmat))
  } else { 
    spnames <- colnames(coocmat)
  }
  
  tab <- data.frame(expand.grid(spi = spnames, spj = spnames), as.vector(coocmat))
  names(tab) <- c('spi', 'spj', name)
  
  if ( full ) { 
    tab <- data.frame(tab, 
                      pval = as.vector(attr(coocmat, "pvals")), 
                      obs  = as.vector(attr(coocmat, "obs")))
  }
  
  if ( trim ) { 
    tab <- tab[!tab[ ,"ses"] %in% no_inter_value, ]
  }
  
  return(tab)
}

matricize <- function(tab, no_inter_value = NA_real_, value.var = "ses", 
                      fun.aggregate = mean) { 
  require(reshape2)
  
  acast(tab, spi ~ spj, value.var = value.var, fill = no_inter_value, 
        fun.aggregate = mean)
} 
# This file contains cod that will convert an adjacency matrix to df 
#   and vice-versa.
# 

to_adjmat <- function(df, i, j, edgelist, fill = 0) { 
  
  acast(df, list(as.quoted(substitute(i)), 
                 as.quoted(substitute(j))),
        value.var = as.character(substitute(edgelist)), 
        fill = fill)
  
}

