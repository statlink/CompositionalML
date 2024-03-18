alfappr.tune <- function(y, x, a = seq(-1, 1, by = 0.1), nterms = 1:10, ncores = 1,
                         folds = NULL, nfolds = 10, seed = NULL, graph = FALSE) {

  if ( min(x) == 0 )  a <- a[a > 0]

  if (ncores > 1) {
    runtime <- proc.time()
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    folds <- Compositional::makefolds(y, nfolds = nfolds, stratified = FALSE, seed = seed )
    nfolds <- length(folds)

    ww <- foreach(k = 1:length(a), .combine = cbind, .export = c(".alfapprtune", "ppr",
                  "colmses", "colmeans"), .packages = "Rfast") %dopar% {
      z <- Compositional::alfa(x, a[k])$aff
      per <- as.numeric( .alfapprtune(y, z, nterms = nterms, folds = folds)$perf )
      return(per)
    }
    parallel::stopCluster(cl)
    per <- ww
    runtime <- proc.time() - runtime

  } else {
    if ( is.null(folds) )  folds <- Compositional::makefolds(y, nfolds = nfolds, stratified = FALSE, seed = seed)
    nfolds <- length(folds)
    per <- matrix(nrow = length(a), ncol = 2)
    runtime <- proc.time()
    for ( k in 1:length(a) ) {
      z <- Compositional::alfa(x, a[k])$aff
      z <- as.data.frame(z)
      per[k, ] <- as.numeric( .alfapprtune(y, z, nterms = nterms, folds = folds)$perf )
    }
    runtime <- proc.time() - runtime
  }

  if (graph) {
    plot(a, per[, 2], type = "b", ylim = c( min(per[, 2]), max(per[, 2]) ), ylab = "Estimated performance",
         xlab = expression( paste(alpha, " values") ), cex.lab = 1.2, cex.axis = 1.2, pch = 16, col = "green")
    abline(v = a, col = "lightgrey", lty = 2)
    abline(h = seq(min(per[ ,2]), max(per[, 2]), length = 10), col = "lightgrey", lty = 2)
  }

  rownames(per) <- paste("alpha=", a, sep = "")
  colnames(per) <- c("nterms", "performance")
  ind <- which.min(per[, 2])

  list(per = per, performance = per[ind, 2], best_a = a[ind], runtime = runtime)
}





.alfapprtune <- function(y, x, nterms = c(2:3), folds = NULL) {

  nfolds <- length(folds)
  per <- matrix( nrow = nfolds, ncol = length(nterms) )

  runtime <- proc.time()
  for ( k in 1:nfolds ) {
    ytrain <- y[ -folds[[ k ]] ]
    ytest <- y[ folds[[ k ]] ]
    xtrain <- x[-folds[[ k ]], ]
    xtest <- x[folds[[ k ]], ]
    colnames(xtest) <- colnames(xtrain)
    st <- .alfappr(xtest, ytrain, xtrain, nterms)$est1
    per[k, ] <- Rfast2::colmses(ytest, st)
  }  ##  end  for (k in 1:nfolds) {

  runtime <- proc.time() - runtime
  per <- cbind(nterms, Rfast::colmeans(per) )
  colnames(per) <- c("nterms", "mse")
  ind <- which.min(per[, 2])

  list(per = per, perf = per[ind, ])
}
