
# Bin a transect data.frame
bin_data <- function(df, 
                     binsize = 5, 
                     bins = NULL, 
                     xmax = max(df[ ,"xe"]),
                     return_bin = FALSE, 
                     logical = TRUE, 
                     spnames = 'attribut') { 
  
  spp <- unique(df[ ,spnames])
  
  # Bin and compute number of assoc
  if ( is.null(bins) ) { 
    bins <- seq(0, xmax, by = binsize)
  } else { 
    xmax <- max(bins)
  }
  
  spfact <- as.factor(df[ ,spnames])
  
  result <- bin_data_core(df[ ,"xi"], df[ ,"xe"], spfact, bins, binsize)
  
  if (logical) { 
    result <- result > 0 
  }
  result <- as.data.frame(result)
  names(result) <- levels(spfact)
  
  if (return_bin) { 
    return( data.frame(bin = bins, result )) 
  } else { 
    return(result)
  }
}

sample_transect <- function(df, N, binsize = 5, xmax = max(df[ ,"xe"]), 
                            return_bin = FALSE, 
                            spnames = "attribut", 
                            transectid = "transect") { 
  
  if (is.null(transectid)) { 
    transectid <- "transect_default"
    df[ ,transectid] <- "1"
  }
  
  result <- sample_transect_core(N, df[ ,"xi"], df[ ,"xe"], 
                                 as.integer(as.factor(df[ ,spnames])), 
                                 as.integer(as.factor(df[ ,transectid])), 
                                 xmax = xmax, binsize = binsize)
  
  colnames(result) <- c('transect', 'bin', levels(as.factor(df[ ,spnames])))
  
  if (!return_bin) { 
    result <- result[ ,c(-1, -2)]
  }
  
  return(result)
}

shuffle <- function(df, xmax = NULL, 
                    method = "free", 
                    groups = NULL, 
                    perms = NULL, 
                    randomize_groups = TRUE) {
  
  if ( is.null(groups) ) { 
    df[ ,"groups"] <- rep("a", nrow(df))
    groups <- "groups"
  } 
  
  if ( ! groups %in% names(df) ) { 
    stop("Missing group variable in transect data: ", groups)
  }
  
  # Randomize the groups if asked for it
  if ( randomize_groups ) { 
    df[ ,groups] <- sample(df[ ,groups])
  }
  
  for ( grp in unique(df[ ,groups]) ) { 
    
    # Init 
    ingrp <- df[ ,groups] == grp
    
    if ( is.null(xmax) ) { 
      xmax <- max(df[ , "xe"])
    }
    
    if ( method == "free" ) { 
      shuffletab <- cbind(df[ingrp, 'xi'], df[ingrp, 'xe'])
      newxs <- shuffle_core(shuffletab, xmax)
      
    } else if ( method == "perms" ) { 
      if ( ! perms %in% names(df) ) { 
        stop("the column ", perms, " must exist in the data")
      }
      
      # Shuffle with permeabilities
      shuffletab <- as.matrix(df[ingrp, c("xi", "xe", perms)])
      newxs <- shuffle_perms(shuffletab, xmax = xmax)
      
    } else { 
      stop('Method not implemented')
    }
    
    df[ingrp, "xi"] <- newxs[ ,1]
    df[ingrp, "xe"] <- newxs[ ,2]
  }
  
  # Drop groups column if needed, i.e. when we did not pass a groups column name
  if ( is.null(groups) ) { 
    df <- df[ , ! names(df) %in% "groups"]
  }
  
  return(df)
}

merge_same_species <- function(df, spnames = "attribut") { 
  ddply(df, spnames, merge_if_overlap)
}

merge_if_overlap <- function(df) { 
  df <- df[order(df[ ,"xi"]), ]
  
  removed_this_pass <- 1
  while ( removed_this_pass > 0 ) { 
    removed_this_pass <- 0
    
    i <- 1
    while (i < nrow(df)) { 
      j <- i+1
      if ( overlap(df[i, "xi"], df[j,"xi"], df[i,"xe"], df[j,"xe"]) > 0 ) { 
        # We merge 
        df[j, "xi"] <- min(df[c(i,j), "xi"])
        df[j, "xe"] <- max(df[c(i,j), "xe"])
        if ( "height" %in% names(df) ) { 
          df[j, "height"] <- max(df[c(i,j), "height"])
        }
        df <- df[-i, ]
        removed_this_pass <- removed_this_pass + 1
      }
      
      i <- i + 1
    }
  }
  
  return(df)
}
