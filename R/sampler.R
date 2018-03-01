library(GIGrvg) # Much better than generalized hyperbolic pacakge
library(coda)
library(copula)

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
    # if (TRUE) {
    if (M.max != 0 & 1/M.max > 1) {

      z <- rexp(1, alpha.max) + lower[i]
      u <- runif(1, 0, 1)
      dtexp <- alpha.max*exp(-alpha.max*(z - lower[i]))
      while (u > dgig(z, lambda = lambda[i], chi = chi[i], psi = psi[i])/(M.max*dtexp)) {
        z <- rexp(1, alpha.max) + lower[i]
        u <- runif(1, 0, 1)
        dtexp <- alpha.max*exp(-alpha.max*(z - lower[i]))
      }
    } else {

      z <- rgig(1, lambda = lambda[i], chi = chi[i], psi = psi[i])
      while (z < lower[i]) {
        z <- rgig(1, lambda = lambda[i], chi = chi[i], psi = psi[i])
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
      dgig.max.log <- dgig(upper[i], lambda = lambda[i], chi = chi[i], psi = psi[i], log = TRUE)
    } else {
      dgig.max.log <- dgig(m, lambda = lambda[i], chi = chi[i], psi = psi[i], log = TRUE)
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

      while (log(u) > (dgig(z, lambda = lambda[i], chi = chi[i], psi = psi[i], log = TRUE) - M.log - log.upp)) {
        z <- runif(1, 0, upper[i])
        u <- runif(1, 0, 1)
      }
    } else {

      z <- rgig(1, lambda = lambda[i], chi = chi[i], psi = psi[i])
      while (z > upper[i]) {
        z <- rgig(1, lambda = lambda[i], chi = chi[i], psi = psi[i])
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
  V.rt <- tcrossprod(tcrossprod(V.inv.eig$vectors, diag(sqrt(ifelse(1/V.inv.eig$values > 0, 1/V.inv.eig$values, 0)))), V.inv.eig$vectors)
  V <- tcrossprod(tcrossprod(V.inv.eig$vectors, diag(ifelse(1/V.inv.eig$values > 0, 1/V.inv.eig$values, 0))), V.inv.eig$vectors)
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
    chi <- ifelse(chi < 10^(-14), 10^(-14), chi)
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

    s.sq[j] <- 1/(2*retstable(q/2, V0 = 1, h = lambda))
    while (!check.cond(u = u, beta = beta, Omega.inv, s.sq = s.sq, j = j)) {
      s.sq[j] <- 1/(2*retstable(q/2, V0 = 1, h = lambda))

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

#' @export
mcmc.ssp <- function(X, y, Sigma, sigma.sq, prior = "sng", c = NULL, q = NULL, m = 2,
                     num.samp, burn.in = 0, thin = 1, print.iter = FALSE, str = "uns") {

  # X can be an n \times t \times r array, in which case beta is t \times r, allow unstructured
  # correlation along r dimension

  n <- dim(X)[1]
  t <- dim(X)[2]
  r <- dim(X)[3]

  # Set up data
  Z <- t(apply(X, 1, "c"))

  ZtZ <- crossprod(Z)
  Zty <- crossprod(Z, y)

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
  if (!is.null(Sigma)) {
    if (prior != "spn") {
      Omega <- Sigma/els
      # Project Omega back into positive semidefinite cone
      Omega.eig <- eigen(Omega)
      Omega <- tcrossprod(tcrossprod(Omega.eig$vectors, diag(ifelse(Omega.eig$values > 0, Omega.eig$values, 0))), Omega.eig$vectors)
      Omega.inv <- tcrossprod(tcrossprod(Omega.eig$vectors, diag(ifelse(Omega.eig$values > 0, 1/Omega.eig$values, 0))), Omega.eig$vectors)
      # I've been lazy with notation, Omega actually needs to be (tr)\times(tr)
      Omega <- Omega%x%diag(t)           # There isn't a way around doing this directly
      Omega.inv <- Omega.inv%x%diag(t)   # There isn't a way around doing this directly
    } else {
      Omega <- abs(Sigma)^(1/m)
      Omega.eig <- eigen(Omega)
      Omega.inv <- tcrossprod(tcrossprod(Omega.eig$vectors, diag(ifelse(Omega.eig$values > 0, 1/Omega.eig$values, 0))), Omega.eig$vectors)
      Omega <- Omega%x%diag(t)           # There isn't a way around doing this directly
      Omega.inv <- Omega.inv%x%diag(t)   # There isn't a way around doing this directly

      Psi <- sign(Sigma)*abs(Sigma)^((m - 1)/m)
      Psi.eig <- eigen(Psi)
      Psi.inv <- tcrossprod(tcrossprod(Psi.eig$vectors, diag(ifelse(Psi.eig$values > 0, 1/Psi.eig$values, 0))), Psi.eig$vectors)
      Psi <- Psi%x%diag(t)           # There isn't a way around doing this directly
      Psi.inv <- Psi.inv%x%diag(t)   # There isn't a way around doing this directly

    }
  } else {
    Omega <- diag(r*t)
    Omega.inv <- diag(r*t)
    if (prior == "spn") {
      Psi <- Psi.inv <- diag(r*t)
    }
  }
  betas <- matrix(nrow = num.samp, ncol = ncol(Z))
  ss <- matrix(1, nrow = num.samp, ncol = ncol(Z))
  Sigmas <- array(NA, dim = c(num.samp, r, r))
  s.old <- rep(1, ncol(Z))

  for (i in 1:(burn.in + thin*num.samp)) {

    if (print.iter) {cat("i=", i, "\n")}

    if (prior != "spn") {
      beta <- samp.beta(XtX = ZtZ, Xty = Zty, s.sq = s.old^2,
                        Omega.inv = Omega.inv/tau^2, sig.sq = sigma.sq)
      if (prior == "sng") {
        s <- sqrt(samp.s.sq(beta = beta, Omega.inv = Omega.inv, c = c,
                            d = c, s.sq = s.old^2))
      } else if (prior == "spb") {
        s <- sqrt(samp.s.sq.bridge(beta, Omega.inv = Omega.inv/tau^2, q = q,
                                   eta = eta, s.sq = s.old^2))
      } else if (prior == "sno") {
        s <- rep(1, ncol(Z))
      }

      if (is.null(Sigma)) {
        Omega.inv <- samp.Omega.inv(Beta = matrix(beta, nrow = t, ncol = r)/matrix(s, nrow = t, ncol = r))
        Omega <- solve(Omega.inv)
        if (i > burn.in & (i - burn.in)%%thin == 0) {
          Sigmas[(i - burn.in)/thin, , ] <- cov2cor(Omega*els)
        }
        Omega <- Omega%x%diag(t)           # There isn't a way around doing this directly
        Omega.inv <- Omega.inv%x%diag(t)   # There isn't a way around doing this directly
      }
    } else {
      uv <- sample.uv(old.v = s.old, sigma.sq.z = sigma.sq, Sigma.u.inv = Omega.inv,
                      Sigma.v.inv = Psi.inv, XtX = ZtZ, Xty = Zty)
      beta <- uv[, 1]*uv[, 2]
      s <- uv[, 2]

      if (is.null(Sigma)) {
        Omega.inv <- samp.Omega.inv(Beta = matrix(beta, nrow = t, ncol = r)/matrix(s, nrow = t, ncol = r), str = str)
        Omega <- solve(Omega.inv)
        Psi.inv <- samp.Omega.inv(Beta = matrix(s, nrow = t, ncol = r), str = str)
        Psi <- solve(Psi.inv)
        if (i > burn.in & (i - burn.in)%%thin == 0) {
          Sigmas[(i - burn.in)/thin, , ] <- cov2cor(Omega*Psi)
        }
        Omega <- Omega%x%diag(t)           # There isn't a way around doing this directly
        Omega.inv <- Omega.inv%x%diag(t)   # There isn't a way around doing this directly
        Psi <- Psi%x%diag(t)           # There isn't a way around doing this directly
        Psi.inv <- Psi.inv%x%diag(t)   # There isn't a way around doing this directly
      }
    }

    s.old <- s
    if (i > burn.in & (i - burn.in)%%thin == 0) {
      betas[(i - burn.in)/thin, ] <- beta
      ss[(i - burn.in)/thin, ] <- s
    }
  }

  return(list("betas" = betas, "ss" = ss, "Sigmas" = Sigmas))

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
