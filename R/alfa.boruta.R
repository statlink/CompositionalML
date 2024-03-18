alfa.boruta <- function(y, x,  a = seq(-1, 1, by = 0.1), runs = 100 ) {
  mod <- list()
  if ( min(x) == 0 )  a <- a[a > 0]
  for ( k in 1:length(a) ) {
    z <- Compositional::alfa(x, a[k])$aff
    z <- as.data.frame(z)
    mod[[ k ]] <- Boruta::Boruta(y ~ z, pValue = 0.05, maxRuns = runs)
  }
  names(mod) <- paste("alpha=", a, sep = "")
  mod
}
