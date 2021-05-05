
test_that("Continous overlap metrics work", { 
  
  ex <- data.frame(xi = c(0, 0), 
                   xe = c(0.5, 2), 
                   attribut = c("A", "B"), 
                   groups = c("grp1", "grp1"))
  
  expect_true({ 
    abs(cont_overlap(ex, groups = "groups")[1, 2] - 0.5) <  0.0001
  })
  
  # Check groups taken into account
  ex <- data.frame(xi = c(0, 0), 
                   xe = c(0.5, 2), 
                   attribut = c("A", "B"), 
                   groups = c("grp1", "grp2"))
  
  expect_true({ 
    abs(cont_overlap(ex, groups = "groups")[1, 2] - 0) < 0.0001
  })
    
})
