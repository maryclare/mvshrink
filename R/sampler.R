# Need a function to compute log likelihood for s proposals
#' @export
cond.nb.s.log <- function(y, nb.s, fit, pr.nb.r.a, pr.nb.r.b) {
  nb.r <- exp(nb.s)
  log(choose(y  + nb.r - 1, y)) - nb.r*log(1 + exp(fit)) + pr.nb.r.a*nb.s - pr.nb.r.b*nb.r
}

# Need a function to compute the gig mode function
#' @export
m.fn <- function(alpha, lambda, psi, chi) {
  a <- (psi - 2*alpha)
  # (sqrt((1 - lambda)^2 + a*chi) - (1 - lambda))/a # Lower way of computing this quantity is more stable, based on Hormann and Ledyold (2014)
  (chi)/(sqrt((1 - lambda)^2 + a*chi) + (1 - lambda))
}

# Procedure:
# - Check l \leq m.fn(0, lambda = lambda, psi = psi, chi = chi).
#   - If yes then find psi that maximizes the top numerically over (0, psi/2)
#   - If no, then find the value \tilde{\alpha} such that l = m.fn(\tilde{\alpha}, lambda = lambda, psi = psi, chi = chi).
#     - If \tilde{\alpha} > psi/2, set \tilde{\alpha} = psi/2
#     - Then obtain max. accept over (0, \tilde{alpha}) and compare to max. acceptance given by acc. prob. at psi/2 if \tilde{\alpha} < psi/2

#' @export
M.fun <- function(alpha, lambda, psi, chi, l) {

  # For each value of alpha I want to compute the corresponding
  # value of m
  c <- (psi/chi)^(lambda/2)/(2*besselK(sqrt(psi*chi), lambda))

  m <- m.fn(alpha, lambda = lambda, psi = psi, chi = chi)
  m[alpha == psi/2] <- Inf
  M <- numeric(length(alpha))

  which.less <- which(l <= m)

  for (i in which.less) {
    m[i] <- ((lambda - 1) + sqrt((lambda - 1)^2 + chi*(psi - 2*alpha[i])))/(psi - 2*alpha[i])
    M[i] <- c*m[i]^(lambda - 1)*exp(-((psi - 2*alpha[i])*m[i] + chi/m[i])/2)*exp(-alpha[i]*l)/alpha[i]
  }

  if (length(which.less) > 0) {
    m[!(1:length(m) %in% which.less)] <- l
    M[!(1:length(m) %in% which.less)] <- c*l^(lambda - 1)*exp(-((psi - 2*alpha[i]/2)*l + chi/l)/2)*exp(-alpha[i]*l/2)/(alpha[i]/2)
  }
  return(M)
}

#' @export
rgiglt <- function(p, lambda, chi, psi, lower) {

  x <- numeric(p)
  if (length(lambda) == 1 & p > 1) {
    lambda <- rep(lambda, p)
  }
  if (length(chi) == 1 & p > 1) {
    chi <- rep(chi, p)
  }
  if (length(psi) == 1 & p > 1) {
    psi <- rep(psi, p)
  }
  if (length(lower) == 1 & p > 1) {
    lower <- rep(lower, p)
  }

  for (i in 1:p) {

    # Worked very poorly
    # u <- runif(1, 0, 1)
    # pr.low <- pgig(lower[i], lambda = lambda[i], psi = psi[i], chi = chi[i])
    # const <- u*(1 - pr.low) + pr.low
    # z <- qgig(const, lambda = lambda[i], psi = psi[i], chi = chi[i])

    # Had some trouble with very small alphas in some cases
    alpha <- seq(10^(-2), psi[i]/2, length.out = 1000)
    M <- M.fun(alpha, lambda = lambda[i], psi = psi[i],
               chi = chi[i], l = lower[i])

    max.ind <- which(1/M == max(1/M, na.rm = TRUE))
    if (length(max.ind) > 1) {
      max.ind <- length(M) - 1
    }
    alpha.max <- alpha[max.ind]
    M.max <- M[max.ind]

    # I tried adding a step that checks that the pr(accept) under rejection sampler > the
    # pr(accept) under naive sampler.  Seems to help!
    if (M.max != 0 & 1/M.max > 1) {

      z <- rexp(1, alpha.max) + lower[i]
      u <- runif(1, 0, 1)
      dtexp <- alpha.max*exp(-alpha.max*(z - lower[i]))
      while (u > GIGrvg::dgig(z, lambda = lambda[i], chi = chi[i], psi = psi[i])/(M.max*dtexp)) {
        z <- rexp(1, alpha.max) + lower[i]
        u <- runif(1, 0, 1)
        dtexp <- alpha.max*exp(-alpha.max*(z - lower[i]))
      }
    } else {

      z <- GIGrvg::rgig(1, lambda = lambda[i], chi = chi[i], psi = psi[i])
      while (z < lower[i]) {
        z <- GIGrvg::rgig(1, lambda = lambda[i], chi = chi[i], psi = psi[i])
      }

    }
    x[i] <- z
  }

  return (x)

}


#' Function for simulating from a right-truncated GIG
#' @export
rgigrt <- function(p, lambda, chi, psi, upper) {

  x <- numeric(p)
  if (length(lambda) == 1 & p > 1) {
    lambda <- rep(lambda, p)
  }
  if (length(chi) == 1 & p > 1) {
    chi <- rep(chi, p)
  }
  if (length(psi) == 1 & p > 1) {
    psi <- rep(psi, p)
  }
  if (length(upper) == 1 & p > 1) {
    upper <- rep(upper, p)
  }

  for (i in 1:p) {

    m <- m.fn(alpha = 0, lambda = lambda[i], chi = chi[i], psi = psi[i])


    if (upper[i] < m) {
      dgig.max.log <- GIGrvg::dgig(upper[i], lambda = lambda[i], chi = chi[i], psi = psi[i], log = TRUE)
    } else {
      dgig.max.log <- GIGrvg::dgig(m, lambda = lambda[i], chi = chi[i], psi = psi[i], log = TRUE)
    }

    log.upp <- log(upper[i])
    M.log <- dgig.max.log + log.upp # M satisfies dgig(x, lambda, chi, psi) <= M/upper for all X in range

    # Pr accept
    # p.acc.naive <- pgig(upper = upper[i],
    #                     lambda = lambda[i], chi = chi[i], psi = psi[i])
    #
    # p.acc.smart <- pgig(upper = upper[i],
    #                     lambda = lambda[i], chi = chi[i], psi = psi[i])/M
    # p.acc.smart > p.acc.naive -> 1/M > 1

    # I tried adding a step that checks that the pr(accept) under rejection sampler > the
    # pr(accept) under naive sampler but it doesn't seem to help, may not be implemented correctly
    if (M.log < 0) {
      z <- runif(1, 0, upper[i])
      u <- runif(1, 0, 1)

      while (log(u) > (GIGrvg::dgig(z, lambda = lambda[i], chi = chi[i], psi = psi[i], log = TRUE) - M.log - log.upp)) {
        z <- runif(1, 0, upper[i])
        u <- runif(1, 0, 1)
      }
    } else {

      z <- GIGrvg::rgig(1, lambda = lambda[i], chi = chi[i], psi = psi[i])
      while (z > upper[i]) {
        z <- GIGrvg::rgig(1, lambda = lambda[i], chi = chi[i], psi = psi[i])
      }
    }


    x[i] <- z
  }

  return (x)

}

# samp.beta <- function(XtX, Xty, s.sq, Omega.inv, sig.sq) {
#   V <- solve(XtX/sig.sq + Omega.inv*tcrossprod(1/sqrt(s.sq)))
#   m <- crossprod(V, Xty/sig.sq)
#   V.eig <- eigen(V/2 + t(V)/2)
#   V.rt <- V.eig$vectors%*%diag(sqrt(ifelse(V.eig$values > 0, V.eig$values, 0)))%*%t(V.eig$vectors)
#   return(m + V.rt%*%rnorm(length(s.sq)))
# }

# I think this is fine
#' @export
samp.beta <- function(XtX, Xty, s.sq, Omega.inv, sig.sq) {
  V.inv <- XtX/sig.sq + Omega.inv*tcrossprod(1/sqrt(s.sq))
  V.inv.eig <- eigen(V.inv/2 + t(V.inv)/2)
  V.rt <- tcrossprod(tcrossprod(V.inv.eig$vectors,
                                diag(sqrt(ifelse(1/V.inv.eig$values > 0, 1/V.inv.eig$values, 0)), nrow = length(V.inv.eig$values), ncol = length(V.inv.eig$values))),
                     V.inv.eig$vectors)
  V <- tcrossprod(tcrossprod(V.inv.eig$vectors,
                             diag(ifelse(1/V.inv.eig$values > 0, 1/V.inv.eig$values, 0), nrow = length(V.inv.eig$values), ncol = length(V.inv.eig$values))),
                  V.inv.eig$vectors)
  m <- crossprod(V, Xty/sig.sq)
  return(m + V.rt%*%rnorm(length(s.sq)))
}

# check.cond <- function(u, beta, Omega.inv, s.sq, j = NULL) {
#   A <- tcrossprod(beta/sqrt(s.sq))*Omega.inv
#   D <- A - diag(diag(A))
#   return(u <= exp(-sum(D)/2))
# }

#' @export
check.cond <- function(u, beta, Omega.inv, s.sq, j) {
  A <- tcrossprod(beta/sqrt(s.sq))*Omega.inv
  s.j <- sqrt(s.sq[j])
  D <- A - diag(diag(A))
  d.j <- -sum(s.j*D[j, ])
  D[j, ] <- rep(0, ncol(D))
  D[, j] <- rep(0, ncol(D))
  a.j <- log(u) + sum(D)/2
  # cat(ifelse(a.j > 0, "upper: ", "lower: "), d.j/a.j, "\n")
  return(ifelse(a.j > 0, s.j < d.j/a.j, s.j > d.j/a.j))
}

#' @export
get.lim <- function(u, beta, Omega.inv, s.sq, j) {
  A <- tcrossprod(beta/sqrt(s.sq))*Omega.inv
  s.j <- sqrt(s.sq[j])
  D <- A - diag(diag(A))
  d.j <- -sum(s.j*D[j, ])
  # cat("D[j, ]= ", D[j, ], "\n")
  D[j, ] <- rep(0, ncol(D))
  D[, j] <- rep(0, ncol(D))
  a.j <- log(u) + sum(D)/2
  # cat("s.j= ", s.j, "\n")
  # cat("d.j= ", d.j, "\n")
  # cat("a.j= ", a.j, "\n")

  return(list("ul" = ifelse(a.j > 0, "u", "l"),
              "lim" = ifelse(d.j/a.j < 0, 0, d.j/a.j)))
}

#' @export
samp.s.sq <- function(beta, Omega.inv, c, d, s.sq) {
  p <- length(beta)

  # Sample u's given s.sq's
  A <- tcrossprod(beta/sqrt(s.sq))*Omega.inv
  D <- A - diag(diag(A))
  upper <- ifelse(is.infinite(exp(-sum(D)/2)), 10^(308), exp(-sum(D)/2))
  u <- runif(1, 0, upper)

  for (j in 1:p) {
    # cat("j= ", j, "\n")
    lambda <- c - 1/2
    chi <- beta[j]^2*Omega.inv[j, j]
    # Check numerical isssue
    chi <- ifelse(chi < 10^(-14), 10^(-14), chi) # rgig gives values of infinity for chi < 10^(-14)
    psi <- 2*d

    lim <- get.lim(u = u, beta = beta, Omega.inv, s.sq = s.sq, j = j)
    # cat(ifelse(lim[["ul"]] == "u", "upper= ", "lower= "), lim[["lim"]], "\n")
    # cat("lambda = ", lambda, "\n")
    # cat("chi = ", chi, "\n")
    # cat("psi = ", psi, "\n")
    if (lim[["ul"]] == "l") {
      s.sq[j] <- rgiglt(1, lambda = lambda, chi = chi, psi = psi, lower = lim[["lim"]]^2)
    } else if (lim[["ul"]] == "u") {
      s.sq[j] <- rgigrt(1, lambda = lambda, chi = chi, psi = psi, upper = lim[["lim"]]^2)
    }

    # if (lim[["ul"]] == "l") {
    #   s.sq[j] <- rgig(1, lambda = lambda, chi = chi, psi = psi)
    #   while (!check.cond(u = u, beta = beta, Omega.inv, s.sq = s.sq, j = j)) {
    #     s.sq[j] <- rgig(1, lambda = lambda, chi = chi, psi = psi)
    #   }
    # } else {
    #
    #   s.sq[j] <- rgigrt(1, lambda = lambda, chi = chi, psi = psi, upper = lim[["lim"]]^2)
    #
    # }

    # cat("s.sq[j]=", s.sq[j], "\n")
    # if (lim[["ul"]] == "l") {
    #   if (lambda == 0 & chi < 10^(-6)) {
    #     s.sq[j] <- rexp(1, rate = psi/2)
    #     while (!check.cond(u = u, beta = beta, Omega.inv, s.sq = s.sq, j = j)) {
    #       s.sq[j] <- rexp(1, rate = psi/2)
    #     }
    #   } else {
    #     s.sq[j] <- rgig(1, lambda = lambda, chi = chi, psi = psi)
    #     while (!check.cond(u = u, beta = beta, Omega.inv, s.sq = s.sq, j = j)) {
    #       s.sq[j] <- rgig(1, lambda = lambda, chi = chi, psi = psi)
    #     }
    #   }
    # } else {
    #
    #   s.sq[j] <- rgigrt(1, lambda = lambda, chi = chi, psi = psi, upper = lim[["lim"]]^2)
    #
    # }

    # cat("s.sq[j]=", s.sq[j], "\n")
    # if (lim[["ul"]] == "l") {
    #   if (lambda == 0 & chi < 10^(-6)) {
    #     s.sq[j] <- rexp(1, rate = psi/2)
    #     while (!check.cond(u = u, beta = beta, Omega.inv, s.sq = s.sq, j = j)) {
    #       s.sq[j] <- rexp(1, rate = psi/2)
    #     }
    #   } else {
    #     s.sq[j] <- rgig(1, lambda = lambda, chi = chi, psi = psi)
    #     while (!check.cond(u = u, beta = beta, Omega.inv, s.sq = s.sq, j = j)) {
    #       s.sq[j] <- rgig(1, lambda = lambda, chi = chi, psi = psi)
    #     }
    #   }
    # } else {
    #
    #   s.sq[j] <- rgigrt(1, lambda = lambda, chi = chi, psi = psi, upper = lim[["lim"]]^2)
    #
    # }

    # if (lim[["ul"]] == "l") {
    #   s.sq[j] <- rgig(1, lambda = lambda, chi = chi, psi = psi)
    #   while (!check.cond(u = u, beta = beta, Omega.inv, s.sq = s.sq, j = j)) {
    #     s.sq[j] <- rgig(1, lambda = lambda, chi = chi, psi = psi)
    #   }
    # } else {
    #   s.sq[j] <- rgigrt(1, lambda = lambda, chi = chi, psi = psi, upper = lim[["lim"]]^2)
    #   # s.sq[j] <- rgig(1, lambda = lambda, chi = chi, psi = psi)
    #   # while (!check.cond(u = u, beta = beta, Omega.inv, s.sq = s.sq, j = j)) {
    #   #   s.sq[j] <- rgig(1, lambda = lambda, chi = chi, psi = psi)
    #   # }
    #   # s.sq[j] <- rgig(1, lambda = lambda, chi = chi, psi = psi)
    #   # while (sqrt(s.sq[j]) > lim[["lim"]]) {
    #   #   s.sq[j] <- rgig(1, lambda = lambda, chi = chi, psi = psi)
    #   # }
    #   # cat("s.sq=", s.sq[j], "\n")
    # }
  }

  return(s.sq)

}

#' @export
samp.s.sq.bridge <- function(beta, Omega.inv, q, eta, s.sq) {
  p <- length(beta)

  # Sample u's given s.sq's
  A <- tcrossprod(beta/sqrt(s.sq))*Omega.inv
  D <- A - diag(diag(A))
  u <- runif(1, 0, exp(-sum(D)/2))

  for (j in 1:p) {

    lambda <- beta[j]^2*Omega.inv[j, j]

    s.sq[j] <- 1/(2*copula::retstable(q/2, V0 = 1, h = lambda))
    while (!check.cond(u = u, beta = beta, Omega.inv, s.sq = s.sq, j = j)) {
      s.sq[j] <- 1/(2*copula::retstable(q/2, V0 = 1, h = lambda))

    }
  }

  return(s.sq)

}

#' @export
beta.em <- function(XtX, Xty, sig.sq, Omega.inv, es.i) {
  V <- solve(XtX/sig.sq + Omega.inv*es.i)
  return(crossprod(V, Xty/sig.sq))
}

#' @export
samp.Omega.inv <- function(Beta, pr.V.inv = diag(ncol(Beta)),
                             pr.df = ncol(Beta) + 2, str = "uns") {
  p <- ncol(Beta)
  if (str == "uns") {
    V.inv <- crossprod(Beta) + pr.V.inv
    df <- nrow(Beta) + pr.df
    return(rWishart(1, df, solve(V.inv))[, , 1])
  } else if (str == "het") {
    b <- apply(Beta, 2, function(x) {sum(x^2)})/(2) + diag(pr.V.inv)/2
    a <- rep(nrow(Beta), ncol(Beta))/2 + pr.df/2
    return(diag(rgamma(p, shape = a, rate = b)))
  } else if (str == "con") {
    b <- sum(apply(Beta, 2, function(x) {sum(x^2)}))/(2) + sum(diag(pr.V.inv))/2
    # I'm a little worried about the code below if 'Beta' is a matrix wtih more than 1 column,
    # Should check. I think it works!
    a <- sum(rep(nrow(Beta), ncol(Beta)))/2 + p*pr.df/2
    return(rgamma(1, shape = a, rate = b)*diag(p))
  }
}

#' @export
sample.uv <- function(old.v, sigma.sq.z,
                      Sigma.u.inv, Sigma.v.inv, XtX, Xty) {

  u <- samp.beta(sig.sq = sigma.sq.z, Omega.inv = Sigma.u.inv, XtX = XtX*tcrossprod(old.v),
                 Xty = Xty*old.v, s.sq = rep(1, length(old.v)))
  v <- samp.beta(sig.sq = sigma.sq.z, Omega.inv = Sigma.v.inv, XtX = XtX*tcrossprod(u),
                 Xty = Xty*u, s.sq = rep(1, length(old.v)))

  return(cbind(u, v))

}

#' Function for sampling from the posterior under multivariate normal, structured normal-gamma, structured power/bridge and structured product normal priors
#' @name mcmc.ssp
#' @description Function for sampling values of \eqn{beta} and possibly \eqn{Sigma} from the posterior distribution under multivariate normal, structured normal-gamma, structured power/bridge and structured product normal priors
#'
#'
#' @usage
#' samples <- mcmc.ssp(X = X, y = y, Sigma = diag(dim(X)[3]), sigma.sq = 1, prior = "sno", num.samp = 100)
#'
#' @param X \eqn{n} by \eqn{t} by \eqn{r} array of covariate values, \eqn{y[i] = vec(X_i)'vec(B) + z[i], i = 1,...,n}
#' @param y \eqn{n} by \eqn{1} response vector or \eqn{m} by \eqn{p} response matrix (if matrix, require \eqn{mp=n})
#' @param Sigma either a fixed \eqn{r} by \eqn{r} positive semidefinite covariance matrix (the provided \eqn{\Sigma} can always be used by the multivariate normal and SPN priors, but for the SNG and SPB priors, the provided value of \eqn{\Sigma} may not be achievable and the closest positive definite matrix will be used) or NULL, in which case values of \eqn{\Sigma} are resampled under either an inverse-Wishart prior for \eqn{\Sigma^{-1}}, or \eqn{\Sigma} is assumed to be diagonal with inverse-Wishart distributed elements or \eqn{\Sigma} is proportional to an identity matrix with inverse-gamma scale, the "str" parameter determines the prior used for \eqn{\Sigma}, prior hyperparameters given by the "pr.V.inv" and "pr.df" parameters
#' @param sigma.sq positive constant, error variance
#' @param prior string equal to "sno" for multivariate normal prior, "sng" for normal-gamma prior, "spb" for power/bridge prior, "spn" for normal product prior
#' @param c positive constant, required when "sng" prior is used, default is NULL
#' @param q positive constant, required when "spb" prior is used, default is NULL
#' @param Psi positive semi-definite matrix, used when prior="spn" and \eqn{\Sigma} is provided and fixed
#' @param num.samp integer; total number of posterior samples to return, defaults to 100
#' @param burn.in integer; total number of posterior samples to discard, defaults to 1
#' @param thin integer; number of posterior samples to be discarded between those returned, defaults to 0
#' @param print.iter logical; if TRUE, the function prints a count of the number of samples from the posterior, defaults to FALSE
#' @param str string equal to "uns", "het" or "con" which determines the prior used for Sigma. "con" forces \eqn{\Sigma} to be proportional to the identity matrix, "het" forces \eqn{\Sigma} to be diagonal, "str" allows \eqn{\Sigma} to have arbitrary structure
#' @param pr.V.inv prior scale matrix for inverse-Wishart prior, set by default to \eqn{I_r}, if a diagonal \eqn{\Sigma} is used provide a diagonal matrix with the inverse-gamma rate parameters for each diagonal element of \eqn{\Sigma}, if \eqn{\Sigma \propto I_r} is used, provide the rate parameter for the inverse-gamma prior on the scale divided by r
#' @param pr.df prior degrees of freedom for inverse-Wishart prior, set by default to \eqn{p + 2}, if a diagonal \eqn{\Sigma} is used provide the shape parameter for the inverse-gamma priors on the diagonal elements of \eqn{\Sigma}, if \eqn{\Sigma \propto I_r} is used provide the shape parameter for the inverse-gamma prior on the scale divided by r
#' @param pr.shape prior shape parameter for inverse-gamma prior on the noise variance, set to 3/2 by default
#' @param pr.rate prior rate parameter for inverse-gamma prior on the noise variance, set to 1/2 by default
#' @param W \eqn{n} by \eqn{l} matrix of unpenalized covariates to include in regression
#' @param reg string equal to "linear" for linear regression, "logit" for logistic regression, "nb" for negative binomial regression
#' @param pr.ssqi.V.inv If y has matrix structure, this is the scale matrix for an inverse wishart prior on column covariance matrix
#' @param pr.ssqi.df If y has matrix structure, this is the scale matrix for an inverse wishart on column covariance matrix
#' @param str.ssqi Structure to be used for covariance matrix of columns of Y
#' @param joint.beta logical; if TRUE, the function updates all elements of matrix of penalized regression coefficients jointly, otherwise the function updates the matrix of penalized regression coefficients one row at a time, set to TRUE by default
#' @param nb.r positive scalar; overdispersion parameter for negative binomial regression, defaults to NULL in which case a gamma prior is assumed for \eqn{r} and values of \eqn{r} are resampled
#' @param pr.nb.r.shape positive scalar; shape parameter for gamma prior on \eqn{r}, negative binomial overdisperson parameter. defaults to 1/2
#' @param pr.nb.r.rate positive scalar; rate parameter for gamma prior on \eqn{r}, negative binomial overdispersion parameter. defaults to 1/2
#' @param eps positive scalar or length p vector; variance for random-walk metropolis-hastings for updating log(r), defaults to 0.5
#' @source Based on the material in the working paper "Structured Shrinkage Priors"
#' @examples
#'
#' set.seed(1)
#' library(RColorBrewer)
#' # Simulate some data
#' n <- 30; t <- 5; r <- 4
#' X <- array(rnorm(n*t*r), dim = c(n, t, r))
#' Sigma <- (1 - 0.5)*diag(r) + 0.5*diag(r)
#' beta <- c(matrix(rnorm(t*r), nrow = t, ncol = r)%*%chol(Sigma))
#' y <- t(apply(X, 1, "c"))%*%beta + rnorm(n)
#'
#' # Fit a few models to this data, first fixing Sigma then not
#' mvn.fs <- mcmc.ssp(X = X, y = y, Sigma = Sigma, sigma.sq = 1, prior = "sno")
#' spn.fs <- mcmc.ssp(X = X, y = y, Sigma = Sigma, sigma.sq = 1, prior = "spn", Psi = sqrt(abs(Sigma)))
#' sng.fs <- mcmc.ssp(X = X, y = y, Sigma = Sigma, sigma.sq = 1, prior = "sng", c = 1)
#' spb.fs <- mcmc.ssp(X = X, y = y, Sigma = Sigma, sigma.sq = 1, prior = "spb", q = 1)
#'
#' # Take a look at some trace plots, compare very small sample posterior means
#' par(mfrow = c(2, 2))
#' plot(mvn.fs$betas[, 1], xlab = "", ylab = expression(beta[1]))
#' plot(spn.fs$betas[, 1], xlab = "", ylab = "")
#' plot(sng.fs$betas[, 1], xlab = "Sample", ylab = expression(beta[1]))
#' plot(spb.fs$betas[, 1], xlab = "Sample", ylab = "")
#'
#' par(mfrow = c(1, 3))
#' plot(colMeans(mvn.fs$betas), colMeans(spn.fs$betas), xlab = "MVN", ylab = "SPN")
#' abline(a = 0, b = 1, lty = 2)
#' plot(colMeans(mvn.fs$betas), colMeans(sng.fs$betas), xlab = "MVN", ylab = "SNG")
#' abline(a = 0, b = 1, lty = 2)
#' plot(colMeans(mvn.fs$betas), colMeans(spb.fs$betas), xlab = "MVN", ylab = "SPB")
#' abline(a = 0, b = 1, lty = 2)
#'
#' mvn.vs <- mcmc.ssp(X = X, y = y, Sigma = NULL, sigma.sq = NULL, prior = "sno", str = "uns")
#' spn.vs <- mcmc.ssp(X = X, y = y, Sigma = NULL, sigma.sq = NULL, prior = "spn", str = "uns")
#' sng.vs <- mcmc.ssp(X = X, y = y, Sigma = NULL, sigma.sq = NULL, prior = "sng", str = "uns", c = 1)
#' spb.vs <- mcmc.ssp(X = X, y = y, Sigma = NULL, sigma.sq = NULL, prior = "spb", str = "uns", q = 1)
#'
#' # Plot covariance matrix posterior means
#' cols <- brewer.pal(11, "RdBu")
#' cols <- cols[length(cols):1]
#' breaks <- seq(-1, 1, length.out = 12)
#' par(mfrow = c(2, 2))
#' image(matrix(rowMeans(apply(mvn.vs$Sigmas, 1, function(x) {cov2cor(x)})), nrow = r, ncol = r)[r:1, r:1], breaks = breaks, col = cols, main = "MVN")
#' image(matrix(rowMeans(apply(spn.vs$Sigmas, 1, function(x) {cov2cor(x)})), nrow = r, ncol = r)[r:1, r:1], breaks = breaks, col = cols, main = "SPN")
#' image(matrix(rowMeans(apply(sng.vs$Sigmas, 1, function(x) {cov2cor(x)})), nrow = r, ncol = r)[r:1, r:1], breaks = breaks, col = cols, main = "SNG")
#' image(matrix(rowMeans(apply(spb.vs$Sigmas, 1, function(x) {cov2cor(x)})), nrow = r, ncol = r)[r:1, r:1], breaks = breaks, col = cols, main = "SPB")
#'
#' @export
mcmc.ssp <- function(X, y, Sigma = NULL, sigma.sq = NULL, prior = "sng", c = NULL, q = NULL, Psi = NULL,
                     num.samp = 100, burn.in = 0, thin = 1, print.iter = FALSE, str = "uns",
                     pr.V.inv = diag(dim(X)[3]),
                     pr.df = dim(X)[3] + 2, pr.shape = 3/2, pr.rate = 1/2, W = NULL,
                     reg = "linear", pr.ssqi.V.inv = diag(dim(as.matrix(y))[2]),
                     pr.ssqi.df = dim(as.matrix(y))[2] + 2, str.ssqi = "uns",
                     joint.beta = TRUE, nb.r = NULL, pr.nb.r.shape = 1/2, pr.nb.r.rate = 1/2, eps = 0.5) {

  # X can be an n \times t \times r array, in which case beta is t \times r, allow unstructured
  # correlation along r dimension
  pr.ssqi.V.inv <- pr.ssqi.V.inv # Weirdly this is needed here to "initialize" pr.ssqi.V.inv before making y a vector
  if (length(eps) == 1 & !is.null(ncol(y))) {
    eps <- rep(eps, ncol(y))
  }

  n <- dim(X)[1]
  t <- dim(X)[2]
  r <- dim(X)[3]
  m <- dim(as.matrix(y))[1]
  p <- dim(as.matrix(y))[2]
  y.mat <- matrix(y, nrow = m, ncol = p)
  y <- c(y)
  null.W <- is.null(W)
  null.Sigma <- is.null(Sigma)
  null.sigma.sq <- is.null(sigma.sq)
  null.nb.r <- is.null(nb.r)
  if (null.nb.r) {
    nb.r <- rep(1, p)
  }
  if (!null.sigma.sq) {
    sigma.sq.inv <- solve(sigma.sq)
  }
  if (null.W) {
    W <- Matrix(0, nrow = n, ncol = 1)
  } else {
    W <- Matrix(W, sparse = TRUE)
  }
  if (reg == "logit") {
    offset <- rep(1/2, n)
  } else if (reg == "nb") {
    offset <- (y + nb.r%x%rep(1, m))/2
  } else {
    offset <- rep(0, n)
  }
  l <- dim(W)[2]

  # Transform Sigma to Omega, prep parameters for different priors
  if (prior != "spn") {
    if (prior == "sng") {
      es <- c^(-1/2)*gamma(c + 1/2)/gamma(c)
      tau <- 1
    } else if (prior == "spb") {
      es <- sqrt(pi/2)*gamma(2/q)/(sqrt(gamma(1/q)*gamma(3/q)))
      eta <- sqrt(gamma(3/q)/gamma(1/q))^q
      tau <- eta^(-1/q)
    } else if (prior == "sno") {
      es <- 1
      tau <- 1
    }
    els <- (1 - es^2)*diag(r) + es^2*matrix(1, nrow = r, ncol = r)
  }
  if (!null.Sigma) {
    if (prior != "spn") {
      Omega <- Sigma/els
      # Project Omega back into positive semidefinite cone
      Omega.eig <- eigen(Omega)
      Omega <- tcrossprod(tcrossprod(Omega.eig$vectors, diag(ifelse(Omega.eig$values > 0, Omega.eig$values, 0))), Omega.eig$vectors)
      Omega.inv <- tcrossprod(tcrossprod(Omega.eig$vectors, diag(ifelse(Omega.eig$values > 0, 1/Omega.eig$values, 0))), Omega.eig$vectors)
      # I've been lazy with notation, Omega actually needs to be (tr)\times(tr)
      if (joint.beta) {
        Omega <- Omega%x%diag(t)           # There isn't a way around doing this directly
        Omega.inv <- Omega.inv%x%diag(t)   # There isn't a way around doing this directly
      }
    } else {
      if (is.null(Psi)) {
        Psi <- sign(Sigma)*abs(Sigma)^(1/2)
      }
      Omega <- (Sigma/Psi)
      Omega[Psi == 0] <- 0
      Omega.eig <- eigen(Omega)
      Omega.inv <- tcrossprod(tcrossprod(Omega.eig$vectors, diag(ifelse(Omega.eig$values > 0, 1/Omega.eig$values, 0))), Omega.eig$vectors)
      if (joint.beta) {
        Omega <- Omega%x%diag(t)           # There isn't a way around doing this directly
        Omega.inv <- Omega.inv%x%diag(t)   # There isn't a way around doing this directly
      }

      Psi.eig <- eigen(Psi)
      Psi.inv <- tcrossprod(tcrossprod(Psi.eig$vectors, diag(ifelse(Psi.eig$values > 0, 1/Psi.eig$values, 0))), Psi.eig$vectors)
      if (joint.beta) {
        Psi <- Psi%x%diag(t)           # There isn't a way around doing this directly
        Psi.inv <- Psi.inv%x%diag(t)   # There isn't a way around doing this directly
      }
    }
  } else {
    if (joint.beta) {
      Omega <- diag(r*t)
      Omega.inv <- diag(r*t)
      if (prior == "spn") {
        Psi <- Psi.inv <- diag(r*t)
      }
    } else {
      Omega <- diag(r)
      Omega.inv <- diag(r)
      if (prior == "spn") {
        Psi <- Psi.inv <- diag(r)
      }
    }
  }
  if (null.sigma.sq & reg == "linear" & p > 1) {
    sigma.sq <- diag(p)
    sigma.sq.inv <- solve(sigma.sq)
  } else if (null.sigma.sq | reg == "logit" | reg == "nb") {
    sigma.sq <- 1
  }

  # Set up data
  if (joint.beta) {
    Z <- array(dim = c(n, 1, t*r))
    Z[, 1, ] <- t(apply(X, 1, "c"))
  } else {
    Z <- X
  }
  num.block <- dim(Z)[2]

  if (reg == "linear" & p > 1) {
    ZtZ <- array(dim = c(num.block, num.block, dim(Z)[3], dim(Z)[3]))
    Zty <- array(dim = c(num.block, dim(Z)[3]))
    ZtW <- array(dim = c(num.block, dim(Z)[3], l)); WtZ <- array(dim = c(num.block, l , dim(Z)[3]))

    for (i in 1:num.block) {
      Zty[i, ] <- crossprod(Z[, i, ], crossprod(sigma.sq.inv%x%diag(m), y))[, 1]
      ZtW[i, , ] <- as.matrix(Matrix::crossprod(Z[, i, ], Matrix::crossprod(sigma.sq.inv%x%diag(m), W)))
      WtZ[i, , ] <- t(ZtW[i, , ])
      for (j in 1:num.block) {
        ZtZ[i, j, , ] <- crossprod(Z[, i, ], crossprod(sigma.sq.inv%x%diag(m), Z[, j, ]))
      }
    }

    Wty <- as.matrix(Matrix::crossprod(W, Matrix::crossprod(sigma.sq.inv%x%diag(m), y)))
    WtW <- as.matrix(Matrix::crossprod(W, Matrix::crossprod(sigma.sq.inv%x%diag(m), W)))
  } else {
    ZtZ <- array(dim = c(num.block, num.block, dim(Z)[3], dim(Z)[3]))
    Zty <- array(dim = c(num.block, dim(Z)[3]))
    ZtW <- array(dim = c(num.block, dim(Z)[3], l)); WtZ <- array(dim = c(num.block, l , dim(Z)[3]))

    for (i in 1:num.block) {
      Zty[i, ] <- crossprod(Z[, i, ], y - offset)[, 1]

      ZtW[i, , ] <- as.matrix(Matrix::crossprod(Z[, i, ],  W))
      WtZ[i, , ] <- t(ZtW[i, , ])
      for (j in 1:num.block) {
        ZtZ[i, j, , ] <- crossprod(Z[, i, ], Z[, j, ])
      }
    }

    Wty <- as.matrix(Matrix::crossprod(W, as.matrix(y - offset)))
    WtW <- as.matrix(Matrix::crossprod(W))
  }


  fits <- matrix(NA, nrow = num.samp, ncol = n)
  deltas <- matrix(0, nrow = num.samp, ncol = ncol(W))
  betas <- matrix(nrow = num.samp, ncol = t*r)
  ss <- matrix(1, nrow = num.samp, ncol = t*r)
  Sigmas <- array(NA, dim = c(num.samp, r, r))
  sigma.sqs <- array(dim = c(num.samp, p, p))
  nb.rs <- array(dim = c(num.samp, p))
  acc.rs <- array(dim = c(num.samp, p))
  s.old <- rep(1, t*r)
  s.old.mat <- matrix(s.old, nrow = t, ncol = r)
  beta <- rep(0, t*r)
  beta.mat <- matrix(beta, nrow = t, ncol = r)
  delta <- rep(0, l)
  # Starting value for logistic regression scales
  if (reg == "logit" | reg == "nb") {
    ome <- rep(1, n)
    omeD <- Matrix(0, nrow = n, ncol = n, sparse = TRUE)
    diag(omeD) <- ome
  }

  for (i in 1:(burn.in + thin*num.samp)) {

    if (print.iter) {cat("i=", i, "\n")}

    if (reg == "linear" & p > 1 & null.sigma.sq) {

      for (k in 1:num.block) {
        Zty[k, ] <- crossprod(Z[, k, ], crossprod(sigma.sq.inv%x%diag(m), y))
        ZtW[k, , ] <- as.matrix(Matrix::crossprod(Z[, k, ], Matrix::crossprod(sigma.sq.inv%x%diag(m), W))); WtZ[k, , ] <- t(ZtW[k, , ])
        for (kk in 1:num.block) {
          ZtZ[k, kk, , ] <- crossprod(Z[, k, ], crossprod(sigma.sq.inv%x%diag(m), Z[, kk, ]))
        }
      }

      WtW <- as.matrix(Matrix::crossprod(W, Matrix::crossprod(sigma.sq.inv%x%diag(m), W)))
      Wty <- as.matrix(Matrix::crossprod(W, Matrix::crossprod(sigma.sq.inv%x%diag(m), y)))
    } else if (reg == "logit" | reg == "nb") {

      for (k in 1:num.block) {
        Zty[k, ] <- crossprod(Z[, k, ], y - offset)
        ZtW[k, , ] <- as.matrix(Matrix::crossprod(Z[, k, ], Matrix::crossprod(omeD, W))); WtZ[k, , ] <- t(ZtW[k, , ])
        for (kk in 1:num.block) {
          ZtZ[k, kk, , ] <- as.matrix(Matrix::crossprod(Z[, k, ], Matrix::crossprod(omeD, Z[, kk, ])))

        }
      }
      WtW <- as.matrix(Matrix::crossprod(W, Matrix::crossprod(omeD, W)))
      Wty <- as.matrix(Matrix::crossprod(W, as.matrix(y - offset)))
    }

    if (!null.W) {
      if (num.block == 1) {
        Xty <- Wty - crossprod(as.matrix(WtZ[1, , ]), beta)
      } else {
        Xty <- Wty
        for (k in 1:num.block) {
          Xty <- Wty - crossprod(as.matrix(WtZ[k, , ]), beta.mat[k, ])
        }
      }
      if (reg == "linear" & p > 1) {
        delta <- samp.beta(XtX = WtW, Xty = Xty, s.sq = rep(1, l),
                           Omega.inv = matrix(0, nrow = l, ncol = l), sig.sq = 1)
      } else {
        delta <- samp.beta(XtX = WtW, Xty = Xty, s.sq = rep(1, l),
                           Omega.inv = matrix(0, nrow = l, ncol = l), sig.sq = sigma.sq)
      }
    }

    if (prior != "spn") {

      if (reg == "linear" & p > 1) {
        if (num.block == 1) {
          beta <- samp.beta(XtX = ZtZ[1, 1, , ], Xty = Zty[1, ] - crossprod(t(ZtW[1, , ]), delta), s.sq = s.old^2,
                            Omega.inv = Omega.inv/tau^2, sig.sq = 1)
        } else {
          for (k in 1:num.block) {
            XtX <- ZtZ[k, k, , ]
            Xty <- Zty[k, ] - crossprod(t(ZtW[k, , ]), delta)
            for (kk in (1:num.block)[1:num.block != k]) {
              Xty <- Xty - crossprod(t(ZtZ[k, kk, , ]), beta.mat[kk, ])
            }
            beta.mat[k, ] <- samp.beta(XtX = XtX, Xty = Xty, s.sq = s.old.mat[k, ]^2,
                                       Omega.inv = Omega.inv/tau^2, sig.sq = 1)
          }
          beta <- c(beta.mat)
        }
      } else {
        if (num.block == 1) {
          beta <- samp.beta(XtX = ZtZ[1, 1, , ], Xty = Zty[1, ] - crossprod(t(ZtW[1, , ]), delta), s.sq = s.old^2,
                            Omega.inv = Omega.inv/tau^2, sig.sq = sigma.sq)
        } else {
          for (k in 1:num.block) {
            XtX <- ZtZ[k, k, , ]
            Xty <- Zty[k, ] - crossprod(t(ZtW[k, , ]), delta)
            for (kk in (1:num.block)[1:num.block != k]) {
              Xty <- Xty - crossprod(t(ZtZ[k, kk, , ]), beta.mat[kk, ])
            }
            beta.mat[k, ] <- samp.beta(XtX = XtX, Xty = Xty, s.sq = s.old.mat[k, ]^2,
                                       Omega.inv = Omega.inv/tau^2, sig.sq = sigma.sq)
          }
          beta <- c(beta.mat)
        }
      }
      if (prior == "sng") {
        if (num.block == 1) {
          s <- sqrt(samp.s.sq(beta = beta, Omega.inv = Omega.inv, c = c,
                              d = c, s.sq = s.old^2))
        } else {
          for (k in 1:num.block) {
            s.old.mat[k, ] <- sqrt(samp.s.sq(beta = beta.mat[k, ], Omega.inv = Omega.inv, c = c,
                                             d = c, s.sq = s.old.mat[k, ]^2))
          }
          s <- c(s.old.mat)
        }
      } else if (prior == "spb") {
        if (num.block == 1) {
          s <- sqrt(samp.s.sq.bridge(beta, Omega.inv = Omega.inv/tau^2, q = q,
                                     eta = eta, s.sq = s.old^2))
        } else {
          for (k in 1:num.block) {
            s.old.mat[k, ] <- sqrt(samp.s.sq.bridge(beta.mat[k, ], Omega.inv = Omega.inv/tau^2, q = q,
                                                    eta = eta, s.sq = s.old.mat[k, ]^2))
          }
           s <- c(s.old.mat)
        }
      } else if (prior == "sno") {
        s <- rep(1, t*r); s.old.mat <- matrix(1, nrow = t, ncol = r)
      }

      if (null.Sigma) {
        Omega.inv <- samp.Omega.inv(Beta = matrix(beta, nrow = t, ncol = r)/matrix(s, nrow = t, ncol = r),
                                    pr.V.inv = pr.V.inv,
                                    pr.df = pr.df, str = str)
        Omega <- solve(Omega.inv)
        Sigmas[(i - burn.in)/thin, , ] <- Omega*els
        if (num.block == 1) {
          Omega <- Omega%x%diag(t)           # There isn't a way around doing this directly
          Omega.inv <- Omega.inv%x%diag(t)   # There isn't a way around doing this directly
        }
      }
    } else {

      if (num.block == 1) {
        if (reg == "linear" & p > 1) {
          uv <- sample.uv(old.v = s.old, sigma.sq.z = 1, Sigma.u.inv = Omega.inv,
                          Sigma.v.inv = Psi.inv, XtX = ZtZ[1, 1, , ], Xty = Zty[1, ] - crossprod(t(ZtW[1, , ]), delta))
        } else {
          uv <- sample.uv(old.v = s.old, sigma.sq.z = sigma.sq, Sigma.u.inv = Omega.inv,
                          Sigma.v.inv = Psi.inv, XtX = ZtZ[1, 1, , ], Xty = Zty[1, ] - crossprod(t(ZtW[1, , ]), delta))
        }
        beta <- uv[, 1]*uv[, 2]
        s <- uv[, 2]
      } else {
        for (k in 1:num.block) {

          XtX <- ZtZ[k, k, , ]
          Xty <- Zty[k, ] - crossprod(t(ZtW[k, , ]), delta)
          for (kk in (1:num.block)[1:num.block != k]) {
            Xty <- Xty - crossprod(t(ZtZ[k, kk, , ]), beta.mat[kk, ])
          }

          if (reg == "linear" & p > 1) {
            uv <- sample.uv(old.v = s.old.mat[k, ], sigma.sq.z = 1, Sigma.u.inv = Omega.inv,
                            Sigma.v.inv = Psi.inv, XtX = XtX, Xty = Xty)
          } else {
            uv <- sample.uv(old.v = s.old.mat[k, ], sigma.sq.z = sigma.sq, Sigma.u.inv = Omega.inv,
                            Sigma.v.inv = Psi.inv, XtX = XtX, Xty = Xty)
          }
          beta.mat[k, ] <- uv[, 1]*uv[, 2]
          s.old.mat[k, ] <- uv[, 2]
        }
        beta <- c(beta.mat)
        s <- c(s.old.mat)
      }

      if (null.Sigma) {
        Omega.inv <- samp.Omega.inv(Beta = matrix(beta, nrow = t, ncol = r)/matrix(s, nrow = t, ncol = r),
                                    str = str,
                                    pr.V.inv = pr.V.inv,
                                    pr.df = pr.df)
        Omega <- solve(Omega.inv)
        Psi.inv <- samp.Omega.inv(Beta = matrix(s, nrow = t, ncol = r), str = str)
        Psi <- solve(Psi.inv)
        if (i > burn.in & (i - burn.in)%%thin == 0) {
          Sigmas[(i - burn.in)/thin, , ] <- Omega*Psi
        }
        if (num.block == 1) {
          Omega <- Omega%x%diag(t)           # There isn't a way around doing this directly
          Omega.inv <- Omega.inv%x%diag(t)   # There isn't a way around doing this directly
          Psi <- Psi%x%diag(t)           # There isn't a way around doing this directly
          Psi.inv <- Psi.inv%x%diag(t)   # There isn't a way around doing this directly
        }
      }
    }

    s.old <- s

    if (num.block == 1) {
      fit <- crossprod(t(Z[, 1, ]), beta) + as.matrix(Matrix::crossprod(Matrix::t(W), delta))
    } else {
      fit <- as.matrix(Matrix::crossprod(Matrix::t(W), delta))
      for (k in 1:num.block) {
        fit <- fit + crossprod(t(Z[, k, ]), beta.mat[k, ])
      }
    }
    fit.mat <- matrix(fit, nrow = m, ncol = p)
    if (null.sigma.sq & reg == "linear") {

      res <- y - fit
      if (p == 1) {
        sigma.sq <- 1/rgamma(1, shape = pr.shape + length(res)/2, rate = pr.rate + sum(res^2)/2)

      } else {
        sigma.sq.inv <- samp.Omega.inv(Beta = matrix(res, nrow = m, ncol = p),
                                       pr.V.inv = pr.ssqi.V.inv,
                                       pr.df = pr.ssqi.df, str = str.ssqi)
        sigma.sq <- solve(sigma.sq.inv)
      }
    }

    if (null.nb.r & reg == "nb") {

      acc <- rep(0, p)
      for (k in 1:p) {
        nb.s <- log(nb.r[k])
        nb.s.new <- nb.s + eps[k]*rnorm(1)

        old.lik <- sum(cond.nb.s.log(y = y.mat[, k], nb.s = nb.s, fit = fit.mat[, k], pr.nb.r.a = pr.nb.r.shape, pr.nb.r.b = pr.nb.r.rate))
        new.lik <- sum(cond.nb.s.log(y = y.mat[, k], nb.s = nb.s.new, fit = fit.mat[, k], pr.nb.r.a = pr.nb.r.shape, pr.nb.r.b = pr.nb.r.rate))
        ratio <- min(1, exp(new.lik - old.lik))
        if (runif(1) < ratio) {
         nb.s <- nb.s.new
         nb.r[k] <- exp(nb.s)
         acc[k] <- 1
        }
      }
      offset <- (y + nb.r%x%rep(1, m))/2
    }

    if (reg == "logit" | reg == "nb") {
      ome <- BayesLogit::rpg(n, offset*2, fit)
      diag(omeD) <- ome
    }

    if (i > burn.in & (i - burn.in)%%thin == 0) {
      betas[(i - burn.in)/thin, ] <- beta
      deltas[(i - burn.in)/thin, ] <- delta
      ss[(i - burn.in)/thin, ] <- s
      if (reg == "linear") {
        fits[(i - burn.in)/thin, ] <- fit
      } else if (reg == "logit") {
        fits[(i - burn.in)/thin, ] <- exp(fit)/(1 + exp(fit))
      } else {
        fits[(i - burn.in)/thin, ] <- exp(as.numeric(fit) + log(nb.r%x%rep(1, m)))
      }
      sigma.sqs[(i - burn.in)/thin, , ] <- sigma.sq
      nb.rs[(i - burn.in)/thin, ] <- nb.r
      if (null.nb.r & reg == "nb") {
        acc.rs[(i - burn.in)/thin, ] <- acc
      }
    }
  }

  return(list("betas" = betas, "ss" = ss, "Sigmas" = Sigmas, "sigma.sqs" = sigma.sqs,
              "deltas" = deltas, "fits" = fits, "rs" = nb.rs, "acc.rs" = acc.rs))

}

# Functions for later for resampling AR structured covariances
# sample.rho <- function(old, tau.sq,
#                        rho.old, pr, tune, C.inv.old) {
#
#   p <- length(old)
#
#   # Draw new value of rho
#   z.old <- log(((rho.old + 1)/2)/(1 - (rho.old + 1)/2))
#   z.new <- z.old + tune*rnorm(1)
#   rho.new <- -1 + 2*(exp(z.new)/(1 + exp(z.new)))
#
#   # Compute proposal probabilities
#   pr.old <- dnorm(z.old, z.new, sd = tune, log = TRUE) - log(1 - rho.old^2)
#   pr.new <- dnorm(z.new, z.old, sd = tune, log = TRUE) - log(1 - rho.new^2)
#
#   # .b gives the more standard way of computing the acc. prob, same answer
#   # pr.old.b <- dnorm(z.old, z.new, sd = tune, log = TRUE)
#   # pr.new.b <- dnorm(z.new, z.old, sd = tune, log = TRUE)
#
#   C.inv.new <- diag(p)
#   for (i in 1:p) {
#     if (i %in% c(1, p)) {
#       C.inv.new[i, i] <- (1 - rho.new^2)^(p - 2)
#     } else {
#       C.inv.new[i, i] <- (-1)^(p - 2)*(-1 + rho.new)^(p - 2)*(rho.new + 1)^(p - 2)*(1 + rho.new^2)
#     }
#     if (i < p) {
#       C.inv.new[i, i + 1] <- -rho.new*(1 - rho.new^2)^(p - 2)
#       C.inv.new[i + 1, i] <- -rho.new*(1 - rho.new^2)^(p - 2)
#     }
#   }
#   C.inv.new <- C.inv.new/(1 - rho.new^2)^(p - 1)
#
#   # Compute likelihood
#   ll.old <- -(p - 1)*log((1 - rho.old^2))/2 - tcrossprod(crossprod(old, C.inv.old), t(old))/(2*tau.sq) + dbeta((rho.old + 1)/2, pr, pr, log = TRUE)
#   ll.new <- -(p - 1)*log((1 - rho.new^2))/2 - tcrossprod(crossprod(old, C.inv.new), t(old))/(2*tau.sq) + dbeta((rho.new + 1)/2, pr, pr, log = TRUE)
#
#   # .b gives the more standard way of computing the acc. prob, same answer
#   # ll.old.b <- -(p - 1)*log((1 - rho.old^2))/2 - tcrossprod(crossprod(old, C.inv.old), t(old))/(2*tau.sq*sigma.sq.z) + dbeta((rho.old + 1)/2, pr, pr, log = TRUE) + log(exp(z.old)/(1 + exp(z.old))^2)
#   # ll.new.b <- -(p - 1)*log((1 - rho.new^2))/2 - tcrossprod(crossprod(old, C.inv.new), t(old))/(2*tau.sq*sigma.sq.z) + dbeta((rho.new + 1)/2, pr, pr, log = TRUE) + log(exp(z.new)/(1 + exp(z.new))^2)
#
#   ratio <- min(1, exp(ll.new + pr.old - (ll.old + pr.new)))
#
#   if (!runif(1) < ratio) {
#     rho.new <- rho.old
#     C.inv.new <- C.inv.old
#   }
#
#   return(list("rho" = rho.new,
#               "C.inv" = C.inv.new,
#               "acc" = (rho.new != rho.old)))
#
# }
#
# sample.tau.sq.inv <- function(old, C.inv,
#                               pr.a, pr.b) {
#   p <- length(old)
#   ssr <- tcrossprod(crossprod(old, C.inv), t(old))
#   b <- as.numeric(ssr)/2 + pr.b
#   a <- p/2 + pr.a
#   return(rgamma(1, shape = a, rate = b))
# }
