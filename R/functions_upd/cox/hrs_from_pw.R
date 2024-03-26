hrs_from_pw <- function(pw_cox){
  met_coefs <- coefficients(pw_cox)[c("met_x", "met_y")]
  exp(c(met_coefs[1], met_coefs[1]+met_coefs[2])) |> unname()
}