
context("Test network metrics computation")

test_that("Estrada's heterogeneity is correct", { 
  
  # A star graph has heterogeneity one 
  a <- matrix(FALSE, nrow = 10, ncol = 10)
  a[1, ] <- a[ ,1] <- TRUE
  diag(a) <- FALSE
  expect_true( abs(estrada_heterogen(a) - 1) < 1e-8 )
  
  # A regular network has heterogeneity zero 
  a <- matrix(FALSE, nrow = 10, ncol = 10)
  for (i in 1:10) { 
    a[i, i %% 10+1] <- TRUE
  }
  expect_true( abs(estrada_heterogen(a) - 0) < 1e-8 )
  
  # A random graph has heterogeneity close to zero 
  hets <- replicate(19, { 
    size <- 300
    a <- matrix(runif(size^2) < .4, nrow = size, ncol = size)
    estrada_heterogen(a)
  })
  expect_true({ 
    abs( mean(hets) - 0 ) < 1e-2
  })
  
})
