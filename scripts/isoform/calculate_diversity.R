# Function to calculate isoform diversity

calc.shannon <- function(x) { 
  # get top two proportions
  top2 <- head(sort(x, decreasing = TRUE), 2)
  # apply diversity calculation
  -sum(top2[top2 > 0] * log(top2[top2 > 0]))
}
