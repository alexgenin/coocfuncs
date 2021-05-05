# 
# 
# 
# 

#'@export
with_shuffling <- function(metricf, 
                           transect, 
                           shuffling_method = "free", 
                           shuffling_args = list(), 
                           groups = NULL, 
                           nnull = 499, 
                           xmax = NULL, 
                           ...) { 
  
  # Compute observed value 
  obs <- metricf(transect, groups = groups, ...)
  
  # Compute a distribution of null values 
  nulls <- lapply(seq.int(nnull), function(n) { 
    shuf <- do.call(shuffle, 
                    c(list(transect, 
                           method = shuffling_method, 
                           groups = groups), 
                      shuffling_args, 
                      xmax = xmax
                    ))
    metricf(shuf, groups = groups, ...) 
  })
  nulls <- simplify2array(nulls)
  
  # Compute statistics
  ses <- matrix(NA_real_, ncol = ncol(nulls), nrow = nrow(nulls))
  colnames(ses) <- rownames(ses) <- colnames(obs)
  pvals <- effs <- nullq05 <- nullq95 <- nullmean <- nullsd <- ses
  for ( i in seq.int(ncol(obs)) ) { 
    for ( j in seq.int(ncol(obs)) ) { 
      nullmean[i, j] <- mean(nulls[i, j, ])
      nullsd[i, j] <- sd(nulls[i, j, ])
      pvals[i, j] <- rank(c(obs[i, j], nulls[i, j, ]))[1] / (nnull + 1)
      ses[i, j]   <- ( obs[i, j] - nullmean[i, j] ) / nullsd[i, j]
      effs[i, j]  <- obs[i ,j] / nullmean[i, j]
      nullq05[i, j] <- quantile(nulls[i, j, ], .05)
      nullq95[i, j] <- quantile(nulls[i, j, ], .95)
    }
  }
  
  # Return values. 
  obstab <- tabularize(obs, name = "obs")
  av <- as.vector # temporary alias to be more concise
  data.frame(obstab, nullmean = av(nullmean), nullsd = av(nullsd), 
             ses = av(ses), pvals = av(pvals), effs = av(effs), 
             nullq05 = av(nullq05), nullq95 = av(nullq95))
}

# Summarise a vector of values where the first one is the observed, and 
# the rest is the null distrib
sumnull <- function(X, na.rm = FALSE) { 
  data.frame(null_mean   = mean(X[-1], na.rm = na.rm), 
             null_sd     = sd(X[-1], na.rm = na.rm),
             null_se     = sd(X[-1], na.rm = na.rm) / sqrt(sum(!is.na(X))),
             null_q05    = quantile(X[-1], .05, na.rm = na.rm), 
             null_q25    = quantile(X[-1], .25, na.rm = na.rm), 
             null_median = median(X[-1], na.rm = na.rm), 
             null_q75    = quantile(X[-1], .75, na.rm = na.rm), 
             null_q95    = quantile(X[-1], .95, na.rm = na.rm), 
             null_N      = sum(!is.na(X[-1])), 
             pval = rank(X)[1] / (length(X) - 1), 
             ses  = ( X[1] - mean(X[-1], na.rm = na.rm) ) / 
                      sd(X[-1], na.rm = na.rm), 
             obs = X[1])
}
