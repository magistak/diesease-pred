spline_string <- function(vec, sub_vec, df = 4){
  if(!is.character(vec) | sum(sub_vec %in% vec) != length(sub_vec)){
    stop()
  }
  purrr::map_chr(vec,function(e){
    ifelse(e %in% sub_vec,sprintf("rcs(%s, %d)", e, df),e)
  })
}