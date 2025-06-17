#' gridcv_mec
#'
#' Reprende gridcv, mais utilise les foncions gcv et gcvlv au lieu de gridcv et gridcvlv
#'
#' @import data.table
#' @export


.mat <- function(X, prefix = NULL) {
  if(is.vector(X))
    X <- matrix(X, ncol = 1)
  if(!is.matrix(X))
    X <- as.matrix(X)
  if(is.null(row.names(X)))
    row.names(X) <- seq_len(dim(X)[1])
  if(is.null(prefix)) {
    if(is.null(colnames(X)))
      colnames(X) <- paste("x", seq_len(dim(X)[2]), sep = "")
  }
  else
    colnames(X) <- paste(prefix, seq_len(dim(X)[2]), sep = "")
  X
}


gcv <- function(X, Y, segm, score, fun, pars, verb = TRUE) {
  ## pars = List of named vectors (arguments) involved in the calculation of the score
  Y <- .mat(Y, "y")
  nrep <- length(segm)
  res_rep <- vector("list", length = nrep)
  nco <- length(pars[[1]])
  for(i in seq_len(nrep)) {
    if(verb)
      cat("/ rep=", i, " ", sep = "")
    listsegm <- segm[[i]]
    nsegm <- length(listsegm)
    zres <- vector("list", length = nsegm)
    for(j in seq_len(nsegm)) {
      if(verb)
        cat("segm=", j, " ", sep = "")
      s <- sort(listsegm[[j]])
      if((is.list(X)==TRUE)&(is.matrix(X)==FALSE)&(is.data.frame(X)==FALSE)){
        zres[[j]] <- gridscore(
          lapply(1:length(X), function(k) X[[k]][-s, , drop = FALSE]),
          Y[-s, , drop = FALSE],
          lapply(1:length(X), function(k) X[[k]][s, , drop = FALSE]),
          Y[s, , drop = FALSE],
          score = score, fun = fun, pars = pars)
      }else{
        zres[[j]] <- gridscore(
          X[-s, , drop = FALSE], Y[-s, , drop = FALSE],
          X[s, , drop = FALSE], Y[s, , drop = FALSE],
          score = score, fun = fun, pars = pars)
      }
    }
    zres <- setDF(rbindlist(zres))
    res_rep[[i]] <- cbind(rep = rep(i, nsegm * nco),
                          segm = sort(rep(1:nsegm, nco)), zres)
  }
  res_rep <- setDF(rbindlist(res_rep))
  if(verb)
    cat("/ End. \n\n")
  namy <- colnames(Y)
  nampar <- names(pars)
  res <- aggregate(res_rep[, namy, drop = FALSE],
                   by = res_rep[, nampar, drop = FALSE], FUN = mean)
  list(val = res, val_rep = res_rep)
}

gcvlv <- function(X, Y, segm, score, fun, nlv, pars = NULL, verb = TRUE) {
  ## pars must not contains nlv
  Y <- .mat(Y, "y")
  nrep <- length(segm)
  res_rep <- vector("list", length = nrep)
  nlv <- seq(min(nlv), max(nlv))
  le_nlv <- length(nlv)
  for(i in seq_len(nrep)) {
    if(verb)
      cat("/ rep=", i, " ", sep = "")
    listsegm <- segm[[i]]
    nsegm <- length(listsegm)
    zres <- vector("list", length = nsegm)
    for(j in seq_len(nsegm)) {
      if(verb)
        cat("segm=", j, " ", sep = "")
      s <- sort(listsegm[[j]])
      if((is.list(X)==TRUE)&(is.matrix(X)==FALSE)&(is.data.frame(X)==FALSE)){
        zres[[j]] <- gscorelv_mec(
          lapply(1:length(X), function(k) X[[k]][-s, , drop = FALSE]),
          Y[-s, , drop = FALSE],
          lapply(1:length(X), function(k) X[[k]][s, , drop = FALSE]),
          Y[s, , drop = FALSE],
          score = score, fun = fun, nlv = nlv, pars = pars)
      }else{
        zres[[j]] <- gscorelv_mec(
          X[-s, , drop = FALSE], Y[-s, , drop = FALSE],
          X[s, , drop = FALSE], Y[s, , drop = FALSE],
          score = score, fun = fun, nlv = nlv, pars = pars)
      }
    }

    yp_list <- lapply(zres, function(x) x$yp)
    ypred <- do.call(rbind, yp_list)

    zres = lapply(zres, `[[`, "res")
    zres <- setDF(rbindlist(zres))
    ## Case where pars is empty
    if(is.null(pars)) {
      res_rep[[i]] <- cbind(rep = rep(i, nsegm * le_nlv),
                            segm = sort(rep(1:nsegm, le_nlv)), zres)
    }
    ## End
    else {
      nco <- length(pars[[1]])
      res_rep[[i]] <- cbind(rep = rep(i, nsegm * le_nlv * nco),
                            segm = sort(rep(1:nsegm, le_nlv * nco)), zres)
    }
  }
  res_rep <- setDF(rbindlist(res_rep))
  if(verb)
    cat("/ End. \n\n")
  namy <- colnames(Y)
  nampar <- c("nlv", names(pars))
  res <- aggregate(res_rep[, namy, drop = FALSE],
                   by = res_rep[, nampar, drop = FALSE], FUN = mean)
  list(val = res, val_rep = res_rep, ypred=ypred)
}

gridcvlb <- function(X, Y, segm, score, fun, lb, pars = NULL, verb = TRUE) {
  ## pars = List of named vectors (arguments) involved in the calculation of the score
  ## Must not contains lb
  Y <- .mat(Y, "y")
  nrep <- length(segm)
  res_rep <- vector("list", length = nrep)
  lb <- sort(unique(lb))
  le_lb <- length(lb)
  for(i in seq_len(nrep)) {
    if(verb)
      cat("/ rep=", i, " ", sep = "")
    listsegm <- segm[[i]]
    nsegm <- length(listsegm)
    zres <- vector("list", length = nsegm)
    for(j in seq_len(nsegm)) {
      if(verb)
        cat("segm=", j, " ", sep = "")
      s <- sort(listsegm[[j]])
      if((is.list(X)==TRUE)&(is.matrix(X)==FALSE)&(is.data.frame(X)==FALSE)){
        zres[[j]] <- gridscorelb(
          lapply(1:length(X), function(k) X[[k]][-s, , drop = FALSE]),
          Y[-s, , drop = FALSE],
          lapply(1:length(X), function(k) X[[k]][s, , drop = FALSE]),
          Y[s, , drop = FALSE],
          score = score, fun = fun, lb = lb, pars = pars)
      }else{
        zres[[j]] <- gridscorelb(
          X[-s, , drop = FALSE], Y[-s, , drop = FALSE],
          X[s, , drop = FALSE], Y[s, , drop = FALSE],
          score = score, fun = fun, lb = lb, pars = pars)
      }
    }
    zres <- setDF(rbindlist(zres))
    ## Case where pars is empty
    if(is.null(pars)) {
      res_rep[[i]] <- cbind(rep = rep(i, nsegm * le_lb),
                            segm = sort(rep(1:nsegm, le_lb)), zres)
    }
    ## End
    else {
      nco <- length(pars[[1]])
      res_rep[[i]] <- cbind(rep = rep(i, nsegm * le_lb * nco),
                            segm = sort(rep(1:nsegm, le_lb * nco)), zres)
    }
  }
  res_rep <- setDF(rbindlist(res_rep))
  if(verb)
    cat("/ End. \n\n")
  namy <- colnames(Y)
  nampar <- c("lb", names(pars))
  res <- aggregate(res_rep[, namy, drop = FALSE],
                   by = res_rep[, nampar, drop = FALSE], FUN = mean)
  list(val = res, val_rep = res_rep)
}

