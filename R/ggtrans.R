
#'@export
ggtrans <- function(trans, height = rnorm(nrow(trans), 2, .5), 
                    xmin = min(trans$xi), xmax = max(trans$xe)) { 
  trans[ ,"attribut"] <- as.factor(trans[ ,"attribut"])
  ggplot(trans) + 
    geom_rect(aes(xmin = xi, xmax = xe, ymin = 0, ymax = height, 
                  fill = attribut), alpha = .5) + 
    coord_cartesian(xlim = c(xmin, xmax))
}
