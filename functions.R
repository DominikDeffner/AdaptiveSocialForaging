#Some use
softmax <- function(x){
  sapply(1:length(x), function(y)
    exp(x[y]) / sum(exp(x))
  )
}
