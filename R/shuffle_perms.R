# 
# Shuffle a transect, taking into account permeabilities
# 

shuffle_perms <- function(stab, xmax) { 
  
  transect_found <- FALSE
  
  while ( ! transect_found ) { 
    
    newstab <- stab * NA
    
    ntot <- nrow(stab) # total number of indivs
    
    inds_remain <- sample(seq(1, ntot)) # individuals to add in order
    inds_placed <- c() 
    
    # Add first individual
    to_add <- inds_remain[1]
    size <- (stab[to_add, 2] - stab[to_add, 1])
    newxi <- runif(1, 0, xmax - size)
    newxe <- newxi + size
    # Adjust the number of indivs and the list of indivs remaining
    newstab[to_add, ] <- c(newxi, newxe, stab[to_add, 3])
    inds_placed <- c(inds_placed, to_add)
    inds_remain <- inds_remain[-1]
    
    
    try <- 0
    while ( try < 500 && length(inds_remain) > 0 ) { 
  #     cat("Indivs: \n")
  #     print(inds_placed)
  #     print(inds_remain)
      try <- try + 1 
      
      # Get the line in stab to add
      to_add <- inds_remain[1] 
      
      # stab[to_add, ]
      
      # Get position at which we want the new individual
      size <- (stab[to_add, 2] - stab[to_add, 1])
      newxi <- runif(1, 0, xmax - size)
      newxe <- newxi + size
      
      # Compute total permeability at that place
      p_going_through <- 1
      for ( i in inds_placed ) { 
        # Compute overlap between the new position and current individual
        ov <- overlap(newstab[i, 1], newxi, 
                      newstab[i, 2], newxe)
        # Normalize that by the length of the individual we want to place, and 
        # multiply by the permeability of the individual already placed
        # cat(ov * stab[i , "perms"] / size, "\n")
        # P going through indiv = 1 - p being intercept 
        # P being intercept = ov/size * ( 1 - per2 )
        p_intercept <- ( 1 - stab[i, 3] ) * ( ov > 0 ) # / size ) 
        p_going_through <- p_going_through * ( 1 - p_intercept ) 
      }
  #     
  #     cat(p_going_through, "\n")
      # Test if we add the individual
  #     print(p)
      if ( runif(1) < p_going_through ) { 
  #       browser()
        # Adjust the number of indivs and the list of indivs remaining
        inds_remain <- inds_remain[-1]
        inds_placed <- c(inds_placed, to_add)
        newstab[to_add, ] <- c(newxi, newxe, stab[to_add, 3])
        try <- 1 
  #       browser()
      }
    }
      
    if ( length(inds_remain) == 0 ) { 
      transect_found <- TRUE
      cat("*")
    } else { 
      cat('.')
    }
    
  #     cat(length(inds_remain), "\n")
    
  
  } 
  return(newstab)
}
