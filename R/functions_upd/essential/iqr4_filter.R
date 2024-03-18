iqr4_filter <- function(v){
  l <- quantile(v, probs = c(0,0.25,0.5,0.75,1), na.rm = TRUE)[c(2,3,4)]
  v>l[2]-4*(l[2]-l[1]) & v<l[2]+4*(l[3]-l[2])
}