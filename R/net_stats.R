
# Fit a multiplex SBM to a network
fit_sbm <- function(graph_dat) { 
    
  # Diversity of roles as an sbm fit. Functional diversity as associations and 
  # functional diversity as traits. We fit an multplex SBM to use the 
  # information on both + and - inters. 
  adjmats <- list(plus =  matricize(graph_dat, 0, "weight") > 0, 
                  minus = matricize(graph_dat, 0, "weight") < 0)

  #FIXME: fix the fitting problem by adding a link to the two 
  # adjacency matrix. The fit fails if the layers of the networks do 
  # not share at least a link between two nodes. So we add a random 
  # link to the network to make sure this condition is satisfied. 
  # This may be fixed at a later stage.
  i <- 1; j <- 2; ok <- FALSE
  print(dim(adjmats[[1]]))
  while ( ! ok ) { 
    # If no link on both layers, we add a link there. 
    if ( i != j && 
          isFALSE(adjmats[[1]][i, j]) && 
          isFALSE(adjmats[[2]][i, j]) ) { 
      adjmats[[1]][i, j] <- adjmats[[1]][j, i] <- 
        adjmats[[2]][i, j] <- adjmats[[2]][j, i] <- TRUE
        ok <- TRUE
    }
    i <- sample(nrow(adjmats[[1]]), 1)
    j <- sample(nrow(adjmats[[1]]), 1)
  }

  require(blockmodels)
  try_result <- tryNULL({ 
    sbm <- BM_bernoulli_multiplex("SBM_sym", adjmats, verbosity = 0, 
                                  plotting = "")
    sbm$estimate(1)
    ngrps <- which.max(sbm$ICL)
    iclmax <- sbm$ICL[ngrps]
  })
  if (is.null(try_result) ) { 
    ngrps <- NA
    iclmax <- NA
  }
  
  return(list(ngrps = ngrps, iclmax = iclmax))
}


# Function that retrieves network statistics 
net_stats <- function(dat, pthresh, long_form = FALSE, 
                      no_self = TRUE, do_spinglass = FALSE) { 
  require(igraph)
  require(entropart)
  require(purrr)
  
  if ( any( duplicated(dat[ ,c("spi", "spj")]) ) ) { 
    warning("Computing net stats on networks with duplicate links !")
  }
  
  if ( ! "cover" %in% names(dat) ) { 
    warning("Species covers are missing. They are necessary for some of the ", 
            "network metrics.")
  }
  
  # Remove self-interactions
  if ( no_self ) { 
    dat <- subset(dat, spi != spj)
  }
  
  # Association signs
  # ---------------------------------------------------------------------
  
  # Number of associations 
  qlow <- pthresh/2
  qtop <- 1 - qlow
  ktot <- mean( (dat$pvals > qtop ) | (dat$pvals < qlow) )
  kplus <- mean(dat$pvals > qtop)
  kminus <- mean(dat$pvals < qlow)
  ratio <- (kplus - kminus) / (kplus + kminus)
  # Format data to work on the graph
  dat %>% 
    mutate(weight = ifelse(pvals > qtop, 1, ifelse(pvals < qlow, -1, 0))) %>% 
    dplyr::select(spi, spj, weight) -> 
      graph_dat
  
  # Average degree 
  undir_graph <- 
    graph_from_data_frame(subset(graph_dat, weight != 0)[ ,c('spi', 'spj')], 
                          directed = FALSE)
  avdeg <- mean(degree(undir_graph))
  maxdeg <- max(degree(undir_graph))
  spreaddeg <- max(degree(undir_graph)) - min(degree(undir_graph))
  
  # Maximum positive and maximum negative degree
  adjmat <- matricize(graph_dat, 0, "weight")
  negdegs <- apply(adjmat, 1, function(X) sum(X<0))
  posdegs <- apply(adjmat, 1, function(X) sum(X>0))
  
  # Kurtosis plus and minus
  kurt_pos <- moments::kurtosis(negdegs)
  kurt_neg <- moments::kurtosis(posdegs)
  
  # Evenness plus and minus. Note: here we compute the evennes of P(k) ! It's 
  # a mistake to give directly the vector of species with their degrees, as 
  # vegan::diversity will read the degrees as abundances, which is very much 
  # *not* what we want it to do. We want cells to represent the number of 
  # times we get the degree k in the network (incl. zero)
  ksneg <- sapply(unique(negdegs), function(d) { 
    a <- sum(negdegs == d); names(a) <- d
    a
  })
  even_neg <-  vegan::diversity(ksneg, "shannon") / log( length(ksneg) )
  
  kspos <- sapply(unique(posdegs), function(d) { 
    a <- sum(posdegs == d); names(a) <- d
    a
  })
  even_pos <-  vegan::diversity(kspos, "shannon") / log( length(kspos) )
  
  # Estrada's measure of heterogeneity 
  estrada_pos <- estrada_heterogen(adjmat == 1)
  estrada_neg <- estrada_heterogen(adjmat == -1)
  
  # Abundance-weighted association signs: this reflect the dominance of 
  # negative over positive interactions in a community, assuming dominant 
  # species have a larger effect than rare species. 
  # ---------------------------------------------------------------------
  if ( ! is.null(dat$cover) ) { 
    
    # FIXME: this is brittle ! 
    cover_dat <- ddply(dat, ~ spi, summarise, cover = cover[1]) 
    # covers <- ddply(dat, ~ spi + spj, summarise, cover = mean(cover))
    
    # Add spi and spj cover to pairwise information
    cover_i <- mutate(cover_dat, coveri = cover)[ ,c("spi", "coveri")]
    cover_j <- mutate(cover_dat, spj = spi, coverj = cover)[ ,c("spj", "coverj")]
    graph_dat_cover <- join(graph_dat, cover_i, by = "spi", match = "first")
    graph_dat_cover <- join(graph_dat_cover, cover_j, by = "spj", match = "first")
    
    wKp <- with(graph_dat_cover, 
                weighted.mean(weight > 0, w = coveri + coverj))
    wKn <- with(graph_dat_cover, 
                weighted.mean(weight < 0, w = coveri + coverj))
    wKtot <- with(graph_dat_cover, 
                weighted.mean(weight != 0, w = coveri + coverj))
    wratio <- (wKp - wKn) / (wKp + wKn)
    
  } else { 
    wKp <- wKn <- wKtot <- wratio <- NA
  }
  
  # Interaction diversity. We use the index by Leinster and Cobbold (2012) 
  # that is based on a distance matrix. 
  # ---------------------------------------------------------------------
  if ( !is.null(dat$cover) ) { 
    maxdist <- 2 * length(unique(graph_dat$spi)) # maximum possible distance
    Z <- matricize(graph_dat, 0, "weight") %>% 
          dist() %>% as.matrix() %>% 
          (function(x) (max(x) - x)/max(x))
    
    if ( any(is.na(Z)) || any(Z < 0) || all(Z < 0) ) { 
      interaction_diver <- 0
    } else { 
      cover_dat <- ddply(dat, ~ spi, summarise, cover = cover[1]) %>% 
                    mutate(cover = cover / sum(cover))
      interaction_diver <- Dqz(cover_dat[ ,'cover'], q = 0, Z)
    }
  } else { 
    interaction_diver <- NA
  } 
  
  # Mean interactional distance 
  mean_dist <- matricize(graph_dat, 0, "weight") %>% 
                 dist() %>% mean()
  
  # Network structure indexes 
  # ---------------------------------------------------------------------
  
  # Make clustering based on spinglass, and get modularity index 
  if ( do_spinglass ) { 
    as.intf <- function(x, l) as.integer(factor(x, levels = l))
    # Crate igraph object to use spinglass 
    ig <- graph_dat %>% 
      subset( as.intf(spi, l = unique(spi)) <= as.intf(spj, l = unique(spi)) ) %>% 
      graph_from_data_frame(directed = FALSE) # this uses the field "weight"
    clust <- cluster_spinglass(ig)
    sg_modularity <- clust$mod
  } else { 
    sg_modularity <- NA_real_
  }
  
  # Make igraph object for use with signnet
  a <- tryNULL({ 
    ig_dat <- subset(graph_dat, weight != 0)
    ig <- graph_from_data_frame(ig_dat, directed = FALSE)
    E(ig)$sign <- ig_dat[ ,"weight"]
    # Compute balance scores
    bal_score_frustr <- signnet::balance_score(ig, method = "frustration") # long
    bal_score_walk   <- signnet::balance_score(ig, method = "walk")
    bal_score_tria   <- signnet::balance_score(ig, method = "triangle")
  })
  if ( is.null(a) ) { 
    bal_score_frustr <- NaN
    bal_score_walk   <- NaN
    bal_score_tria   <- NaN
  }
  
  ans <- list(kplus = kplus, 
              kminus = kminus, 
              ktot = ktot, 
              wktot = wKtot, 
              ratio = ratio, 
              wkplus = wKp, 
              wkminus = wKn, 
              wratio = wratio, 
              average_degree = avdeg, 
              max_degree = maxdeg, 
              max_pos_degree = max(posdegs), 
              max_neg_degree = max(negdegs), 
              spread_deg = spreaddeg, 
              pos_spread_deg = max(posdegs) - min(posdegs),  
              neg_spread_deg = max(negdegs) - min(negdegs), 
              var_pos = var(posdegs), 
              var_neg = var(negdegs), 
              kurt_neg = kurt_neg, 
              kurt_pos = kurt_pos, 
              even_neg = even_neg, 
              even_pos = even_pos, 
              estrada_pos = estrada_pos, 
              estrada_neg = estrada_neg, 
              pos_sum_invecdf = invecdf_get_sum(posdegs), 
              neg_sum_invecdf = invecdf_get_sum(negdegs), 
              bal_score_frustr = bal_score_frustr, 
              bal_score_walk = bal_score_walk, 
              bal_score_tria = bal_score_tria, 
              inter_diver = interaction_diver, 
              mean_dist = mean_dist, 
              spinglass_mod = sg_modularity)
  
  if (long_form) { 
    ans <- data.frame(metric = names(ans), 
                      value  = unlist(ans))
  } else { 
    ans <- as.data.frame(ans)
  }
  
  return(ans)
}

estrada_heterogen <- function(m) { 
  if ( ! is.logical(m) || sum(m) == 0 ) { 
    return(NA)
  }
  m <- m[apply(m, 1, any), apply(m, 2, any)]
  degs <- apply(m, 1, sum)
  # Get edge list and remove directionality as we are dealing with an undirected
  # graph
  edgel <- which(m, arr.ind = TRUE)
  dups <- duplicated(apply(edgel, 1, function(X) paste(sort(X), collapse = "")))
  edgel <- edgel[ ! dups, , drop = FALSE]
  rhos <- plyr::laply(split(edgel, seq.int(nrow(edgel))), function(edge) { 
    ( 1 / sqrt(degs[edge[1]]) - 1 / sqrt(degs[edge[2]]) ) ^ 2
  })
  rho <- sum(rhos)
  N <- nrow(m) # number of nodes
  rho_n <- rho / ( N - 2 * sqrt(N - 1) )
  rho_n
}

invecdf_get_sum <- function(degs) { 
  if ( ! is.finite(min(degs)) || 
       ! is.finite(max(degs)) ) { 
    return(NA_real_)
  }
  uniqdegs <- seq(min(degs), max(degs), by = 1L)
  inv.ecdf <- map_dbl(uniqdegs, function(d) { 
    mean(degs <= d)
  }) 
  sum(inv.ecdf)
}
