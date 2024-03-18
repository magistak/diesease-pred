source_functions <- function(libs = c("cox", "essential")){
  paths <- paste0("R/functions_upd/", libs)
  purrr::walk(paths, function(path){
  files <- list.files(path)
    purrr::walk(files, function(file){
      file = paste0(path, "/", file)
      source(file)
    })
  }
  )
}
