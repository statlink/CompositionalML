alfa.ppr <- function(xnew, y, x, a = seq(-1, 1, by = 0.1), nterms = 1:10) {

  if ( min(x) == 0 )  a <- a[a > 0]
  mod <- list()
  xnew <- as.matrix(xnew)
  if ( dim(xnew)[1] == 1 )  xnew <- t(xnew)

  for ( k in 1:length(a) ) {
    z <- Compositional::alfa(x, a[k])$aff
    z <- as.data.frame(z)
    znew <- Compositional::alfa(xnew, a[k])$aff
    znew <- as.data.frame(znew)
    colnames(znew) <- colnames(z)
    mod[[ k ]] <- .alfappr(znew, y, z, nterms = nterms)
  }

  names(mod) <- paste("alpha=", a, sep = "")
  mod
}




.alfappr <- function(xnew, y, x, nterms = 1:10) {
  len <- length(nterms)
  est1 <- matrix(nrow = dim(xnew)[1], ncol = len)
  colnames(est1) <- paste("nterms=", nterms, sep = "")
  mod <- list()
  for ( i in 1:len ) {
    b <- ppr(y ~., data = x, nterms = nterms[i])
    est1[, i] <- predict(b, newdata = xnew)
    mod[[ i ]] <- b
  }
  list(mod = mod, est1 = est1)
}
