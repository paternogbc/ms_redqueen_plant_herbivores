# Predict from log(x) linear regression------------------------------------------
pred_log2x <- function(a, x, b){
  y <- a + log2(x)*(b)
  return(y)
}

# Predict from log(x) linear regression------------------------------------------
pred_x <- function(a, x, b){
  y <- a + (x)*(b)
  return(y)
}

# Tidy table for phylolm models-------------------------------------------------
tidy.phylolm <- function(x){
  as_tibble(
    cbind(
    term = rownames(summary(x)[[2]]),
    as.data.frame(summary(x)[[2]])
  )
  )
}
