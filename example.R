# 
# Example script for the analysis of continuous transects 
# 

# Load package 
devtools::load_all(".") # path to package

# Example data. Mind the column names which must match in your data 
example_data <- data.frame(
  # Start of plant indivs 
  xi = c(0, 1, 3, 8), 
  # End of plant indivs 
  xe = c(1.5, 1.5, 4, 10), 
  # Species name or attribute (rock, sand, etc.)
  attribut = c("sp1", "sp2", "sp1", "sp3"), 
  # Height of the plant/attribute (optional)
  height = c(1, 1, 3, 4))

# Visualize the transect data. You can give heights to the individuals, 
# otherwise it's random for the sake of representation 
ggtrans(example_data, height = example_data[ ,"height"])

# Compute observed pairwise associations 
obs_coocs <- with_shuffling(cont_overlap, example_data, 
                            # Number of randomizations 
                            # (>= 2000 is good for final runs)
                            nnull = 499)

# Compute statistics from the associations above. The most interesting ones 
# are kplus/kminus/ktot (proportion of pos/neg/total associations, and 
# pos_spread_deg/neg_spread_deg, the heterogeneity (~= variance) for the positive 
# and negative associations, see the oikos paper for the definition).
obs_net_stats <- net_stats(obs_coocs, 
                           # Significance threshold (0.25 is a good start)
                           pthresh = 0.25)

# Now we compute the expected number of pos/neg associations given the covers and 
# abundances of species in your transect
null_expectation <- ldply(seq.int(99), function(nrep) { 
  null_associations <- with_shuffling(cont_overlap, shuffle(example_data), 
                                      nnull = 99)
  null_net_stats <- net_stats(null_associations, pthresh = 0.25)
  data.frame(nrep = nrep, null_net_stats)
}, .progress = "time")

# We compare that to the observed value 
effect_sizes <- ldply(c("kplus", "kminus", "ktot"), function(metric) { 
  stats <- sumnull(c(obs_net_stats[ ,metric], null_expectation[ ,metric]))
  data.frame(metric = metric, stats)
})

ggplot(effect_sizes, aes(x = metric)) + 
  geom_linerange(aes(ymin = null_q25, ymax = null_q75), size = 1) + 
  geom_linerange(aes(ymin = null_q05, ymax = null_q95)) + 
  geom_point(aes(y = null_mean)) + 
  geom_point(aes(y = obs), color = "red") + 
  labs(x = "Network metric", 
       y = "Value", 
       caption = "Black point/bars represent the null distribution, red point the observed value")
