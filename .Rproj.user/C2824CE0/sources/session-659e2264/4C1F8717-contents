alfarf.tune <- function(y, x, a = seq(-1, 1, by = 0.1), size = c(1, 2, 3), depth = c(0, 1), splits = 2:5,
                        R = 500, ncores = 1, folds = NULL, nfolds = 10, stratified = TRUE, seed = NULL, graph = FALSE) {

  if ( min(x) == 0 )  a <- a[a > 0]
  config <- as.matrix( expand.grid(size = size, depth = depth, splits = splits, R = R) )

  if (ncores > 1) {
    runtime <- proc.time()
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    if ( is.factor(y) ) {
      if ( is.null(folds) )  folds <- Compositional::makefolds(y, nfolds = nfolds, stratified = stratified, seed = seed )
      nfolds <- length(folds)
    } else {
      if ( is.null(folds) )  folds <- Compositional::makefolds(y, nfolds = nfolds, stratified = FALSE, seed = seed )
      nfolds <- length(folds)
    }
    ww <- foreach(k = 1:length(a), .combine = cbind, .export = c(".rftune", "rf", "colaccs",
                                                                 "colmses", "colmeans"), .packages = c("ranger", "Rfast", "Rfast2") ) %dopar% {
                                                                   z <- Compositional::alfa(x, a[k])$aff
                                                                   per <- as.numeric( .rftune(y, z, config = config, folds = folds)$perf )
                                                                   return(per)
                                                                 }
    parallel::stopCluster(cl)
    per <- ww
    runtime <- proc.time() - runtime

  } else {

    if ( is.factor(y) ) {
      if ( is.null(folds) )  folds <- Compositional::makefolds(y, nfolds = nfolds, stratified = stratified, seed = seed)
      nfolds <- length(folds)
    } else {
      if ( is.null(folds) )  folds <- Compositional::makefolds(y, nfolds = nfolds, stratified = FALSE, seed = seed)
      nfolds <- length(folds)
    }

    per <- matrix(nrow = length(a), ncol = 3)
    runtime <- proc.time()
    for ( k in 1:length(a) ) {
      z <- Compositional::alfa(x, a[k])$aff
      per[k, ] <- as.numeric( .rftune(y, z, config = config, folds = folds)$perf )
    }
    runtime <- proc.time() - runtime
  }

  if (graph) {
    plot(a, per[, 3], type = "b", ylim = c( min(per[, 3]), max(per[, 3]) ), ylab = "Estimated performance",
         xlab = expression( paste(alpha, " values") ), cex.lab = 1.2, cex.axis = 1.2, pch = 16, col = "green")
    abline(v = a, col = "lightgrey", lty = 2)
    abline(h = seq(min(per[ ,3]), max(per[, 3]), length = 10), col = "lightgrey", lty = 2)
  }

  rownames(per) <- paste("alpha=", a, sep = "")
  colnames(per) <- c("gamma", "cost", "performance")

  if ( is.factor(y) ) {
    ind <- which.max(per[, 3])
  } else ind <- which.min(per[, 3])

  list(per = per, performance = per[ind, 3], best_a = a[ind], runtime = runtime)
}





.rftune <- function(y, x, config, folds = NULL) {

  nfolds <- length(folds)
  p <- dim(config)[1]
  per <- matrix(nrow = nfolds, ncol = p)

  if ( is.factor(y) ) {
    for ( k in 1:nfolds ) {
      ytrain <- y[ -folds[[ k ]] ]
      ytest <- y[ folds[[ k ]] ]
      xtrain <- x[-folds[[ k ]], ]
      xtest <- as.data.frame( x[folds[[ k ]], ] )
      colnames(xtest) <- colnames(xtrain)
      st <- .rfmodel(xtest, ytrain, xtrain, config = config)$est
      per[k, ] <- Rfast2::colaccs( as.numeric(ytest), st)
    }
    per <- cbind(config, Rfast::colmeans(per) )
    ind <- which.max(per[, 5])
    colnames(per) <- c( colnames(config), "acc")

  } else {
    for ( k in 1:nfolds ) {
      ytrain <- y[ -folds[[ k ]] ]
      ytest <- y[ folds[[ k ]] ]
      xtrain <- x[-folds[[ k ]], ]
      xtest <- as.data.frame( x[folds[[ k ]], ] )
      colnames(xtest) <- colnames(xtrain)
      st <- .rfmodel(xtest, ytrain, xtrain, config = config)$est
      per[k, ] <- Rfast2::colmses( as.numeric(ytest), st)
    }
    per <- cbind(config, Rfast::colmeans(per) )
    ind <- which.min(per[, 5])
    colnames(per) <- c( colnames(config), "mse")
  }

  list(per = per, perf = per[ind, ])
}
