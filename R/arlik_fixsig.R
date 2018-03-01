#' @export
arlik.fixsig <- function(par, z, U, d, V, W = NULL, est.rho = TRUE) {

  p <- ncol(U)

  if (!est.rho) {
    par[1] <- 0
  }
  rho <- -1 + 2*exp(par[1])/(1 + exp(par[1]))
  tau.sq <- exp(par[2])

  Sig <- diag(p)

  if (est.rho) {
    for (i in 1:p) {
      for (j in 1:p) {
        Sig[i, j] <- rho^(abs(i - j))
      }
    }
  }
  Sig <- Sig*tau.sq

  if (est.rho) {
    if (is.null(W)) {
      W <- tcrossprod(diag(d), V)
    }
    Var <- tcrossprod(tcrossprod(W, Sig), W) + 1*diag(p)
    Var.ei <- eigen(Var)
    Var.in <- crossprod(t(Var.ei$vectors), tcrossprod(diag(1/Var.ei$values), Var.ei$vectors))
    if (min(Var.ei$values) > 0) {
      nl.lik <- sum(log(Var.ei$values)) + tcrossprod(crossprod(z, Var.in), t(z))
    } else {
      nl.lik <- NA
    }
  } else {
    nl.lik <- sum(log(d^2*tau.sq + 1)) + sum(z^2/(d^2*tau.sq + 1))
  }
  return(as.numeric(nl.lik))
}

#' @export
arlik.b.fixsig <- function(par, b, XtX.inv, est.rho = TRUE) {

  p <- nrow(XtX.inv)

  if (!est.rho) {
    par[1] <- 0
  }
  rho <- -1 + 2*exp(par[1])/(1 + exp(par[1]))
  tau.sq <- exp(par[2])

  Sig <- diag(p)

  if (est.rho) {
    for (i in 1:p) {
      for (j in 1:p) {
        Sig[i, j] <- rho^(abs(i - j))
      }
    }
  }
  Sig <- Sig*tau.sq


  Var <- Sig + 1*XtX.inv
  Var.ei <- eigen(Var)
  Var.in <- tcrossprod(tcrossprod(Var.ei$vectors, diag(1/Var.ei$values)), Var.ei$vectors)
  nl.lik <- sum(log(Var.ei$values)) + tcrossprod(crossprod(b, Var.in), t(b))

  return(as.numeric(nl.lik))
}

#' @export
arlik.y.fixsig <- function(par, y, X, est.rho = TRUE) {

  p <- ncol(X)
  n <- nrow(X)

  if (!est.rho) {
    par[1] <- 0
  }
  rho <- -1 + 2*exp(par[1])/(1 + exp(par[1]))
  tau.sq <- exp(par[2])

  Sig <- diag(p)

  if (est.rho) {
    for (i in 1:p) {
      for (j in 1:p) {
        Sig[i, j] <- rho^(abs(i - j))
      }
    }
  }
  Sig <- Sig*tau.sq

  Var <- tcrossprod(tcrossprod(X, Sig), X) + 1*diag(n)
  Var.ei <- eigen(Var)
  Var.in <- Var.in <- tcrossprod(tcrossprod(Var.ei$vectors, diag(1/Var.ei$values)), Var.ei$vectors)
  nl.lik <- sum(log(Var.ei$values[Var.ei$values > 0])) + tcrossprod(crossprod(y, Var.in), t(y))

  return(as.numeric(nl.lik))
}

#' @export
usear.fixsig <- function(y, X, est.rho = TRUE, use.y = TRUE) {

  XtX <- crossprod(X)
  p <- ncol(X)

  start <- rnorm(2)
  if (use.y |  (rcond(XtX) < .Machine$double.eps)) {
    fit <- optim(par = start, fn = arlik.y.fixsig,
                 y = y, X = X,
                 method = "L-BFGS-B", est.rho = est.rho)
  } else {
    XtX.inv <- solve(XtX)
    b <- crossprod(XtX.inv, crossprod(X, y))
    fit <- optim(par = start, fn = arlik.b.fixsig,
                 b = b, XtX.inv = XtX.inv,
                 method = "L-BFGS-B", est.rho = est.rho)
  }

  if (fit$conv == 0) {
    if (est.rho) {
      rho <-  -1 + 2*exp(fit$par[1])/(1 + exp(fit$par[1]))
    } else {
      rho <- 0
    }
    tau.sq <- exp(fit$par[2])
  } else {
    rho <- tau.sq <- NA
  }

  Sig <- matrix(nrow = p, ncol = p)

  for (i in 1:p) {
    for (j in 1:p) {
      Sig[i, j] <- rho^(abs(i - j))
    }
  }
  Sig <- Sig*tau.sq

  return(list("Sig" = Sig, "tau.sq" = tau.sq, "rho" = rho))

}
