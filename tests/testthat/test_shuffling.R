
test_that("Shuffling works", { 
  
  ex <- data.frame(xi = c(0, 1), 
                   xe = c(1, 2), 
                   attribut = c("A", "B"), 
                   groups = c("grp1", "grp1"))
  
  expect_true({ 
    abs(cont_overlap(ex, groups = "groups")[1, 2] - 0) <  0.0001
  })
  
  # Test that 
  ex_a <- with_shuffling(cont_overlap, ex, shuffling_method = "free", 
                         nnull = 1999)
  # The two species segregate in space
  expect_true({ 
    ex_a[2, "pvals"] < 0.001
  })
  
  # Check that groups are taken into account
  ex <- data.frame(xi = c(0,   0.1, 0.3), 
                   xe = c(0.5, 2, 0.4), 
                   attribut = c("A", "B", "C"), 
                   grp = c("grp1", "grp2", "grp2"))
    
  ex_b <- with_shuffling(cont_overlap, ex, shuffling_method = "free", 
                         nnull = 1999, groups = "grp")
                         
  ex_c <- with_shuffling(cont_overlap, ex, shuffling_method = "free", 
                         nnull = 1999, groups = "grp", 
                         shuffling_args = list(randomize_groups = TRUE))
  
  # Without shuffling groups: pval ~= 0.5 between A & B species as they 
  # never coocur in the null model
  expect_true({ 
    abs(ex_b[2, "pvals"] - 0.5) < 1e-2
  })
  
  # With shuffling groups: pval != 0.5 between the two species as they 
  # coocur sometimes in the null model
  expect_true({ 
    abs(ex_c[2, "pvals"] - 0.5) > 1e-2
  })
    
})
