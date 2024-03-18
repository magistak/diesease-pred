discret_by_people <- function(v, probs){
  if(sum(probs) != 1){
    stop("Sum of probs must be 1")
  }
  breaks <- quantile(v,probs=c(0,cumsum(probs)), na.rm = T) 
  cuts <- cut(v, 
              breaks = breaks,
              include.lowest = TRUE,
              right = FALSE)
  as.factor(cuts)
}