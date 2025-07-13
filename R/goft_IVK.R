#' Goodness-of-fit tests for the cure rate in a mixture cure model
#'
#' @description
#' The aim of this function is to test whether the cure rate \eqn{p}, as a function of the covariates, satisfies a certain parametric model.
#'
#' @param x A numeric vector representing the covariate of interest.
#' @param time A numeric vector of observed survival times.
#' @param delta A numeric vector indicating censoring status (1 = event occurred, 0 = censored).
#' @param model A character string specifying the parametric model for the incidence part. Can be \code{"logit"}, \code{"probit"}, or \code{"cloglog"}.
#' @param theta0 Optional numeric vector with initial values for the model parameters. Default is \code{NULL}.
#' @param nsimb An integer indicating the number of bootstrap replicates.Default is \code{499}.
#' @param h Optional bandwidth value used for nonparametric estimation of the cure rate. Default is \code{NULL}.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{statistic}{Numeric value of the test statistic.}
#'   \item{p.value}{Numeric value of the bootstrap p-value for testing the null hypothesis.}
#'   \item{bandwidth}{The bandwidth used.}
#' }
#'
#' @details
#'  We want to test wether the cure rate \eqn{p}, as a function of covariates, satisfies a certain parametric model, such as, logistic, probit or cloglog model.
#'  The hypothesis are:
#' \deqn{
#'   \mathcal{H}_0 : p = p_{\theta} \quad \text{for some} \quad \theta \in \Theta
#'   \quad \text{vs} \quad
#'   \mathcal{H}_1 : p \neq p_{\theta} \quad \text{for all} \quad \theta \in \Theta,
#' }
#' where \eqn{\Theta} is a finite-dimensional parameter space and \eqn{p_{\theta}} is a known function up to the parameter vector \eqn{\theta}.
#'
#' The test statistic is based on a weighted \eqn{L_2} distance between a nonparametric estimator \eqn{\hat{p}(x)} and a parametric estimator \eqn{p_{\hat{\theta}}(x)} under \eqn{\mathcal{H}_0}:
#' \deqn{
#'   \mathcal{T}_n = n h^{1/2} \int \left(\hat{p}(x) - p_{\hat{\theta}}(x)\right)^2 \pi(x) dx,
#' }
#' where \eqn{\pi(x)} is a known weighting function, often chosen as the covariate density \eqn{f(x)}.
#'
#' A practical empirical version of the statistic is given by:
#' \deqn{
#'   \tilde{\mathcal{T}}_n = n h^{1/2} \frac{1}{n} \sum_{i = 1}^n \left(\hat{p}(x_i) - p_{\hat{\theta}}(x_i)\right)^2,
#' }
#' where the integral is replaced by a sample average.
#'
#' @references
#' Müller, U.U, & Van Keilegom, I. (2019).
#' Goodness-of-fit tests for the cure rate in a mixture cure model.
#' \emph{Biometrika}, 106, 211-227. \doi{10.1093/biomet/asy058}
#'
#' @examples
#'
#' ## Some artificial data
#' set.seed(123)
#' n <- 50
#' x <- runif(n, -2, 2) ## Covariate values
#' y <- rweibull(n, shape = .5*(x + 4)) ## True lifetimes
#' c <- rexp(n) ## Censoring values
#' p <- exp(2*x)/(1 + exp(2*x)) ## Probability of being susceptible
#' u <- runif(n)
#' t <- ifelse(u < p, pmin(y, c), c) ## Observed times
#' d <- ifelse(u < p, ifelse(y < c, 1, 0), 0) ## Uncensoring indicator
#' data <- data.frame(x = x, t = t, d = d)
#'
#' goft(x, t, d, model = "logit")
#' @importFrom smcure smcure coefsmcure
#' @importFrom stats na.omit optim pchisq pnorm runif
#' @importFrom utils capture.output
#' @export

goft <- function(x, time, delta, model = c("logit", "probit", "cloglog"), theta0 = NULL, nsimb = 499, h = NULL) {
  model <- match.arg(model)

  data <- data.frame(Y = time, Delta = delta, x = x)
  data <- data[order(data[,1], 1 - data[,2]), ]
  Y <- data[,1]
  Delta <- data[,2]
  X <- data[,3]
  n <- length(X)

  if (is.null(h)) {
    Ch <- (max(X) - min(X))
    h <- Ch * n^(-1/5)
    g <- Ch * n^(-0.11)
  }else{
    g <- h*n^(0.09)
  }

  ww <- weight_cpp(X, h)
  pnonp <- beran_cpp(as.matrix(data), ww)

  Lik_fun <- switch(model,
                    logit = function(par, X, pnonp) {
                      theta1 <- par[1]
                      theta2 <- par[2]
                      phi <- exp(theta1 + theta2 * X)
                      phi <- phi / (1 + phi)
                      phi <- pmin(pmax(phi, 1e-8), 1 - 1e-8)
                      sum((1 - pnonp) * log(phi) + pnonp * log(1 - phi))
                    },
                    probit = function(par, X, pnonp) {
                      theta1 <- par[1]
                      theta2 <- par[2]
                      phi <- pnorm(theta1 + theta2 * X)
                      phi <- pmin(pmax(phi, 1e-8), 1 - 1e-8)
                      sum((1 - pnonp) * log(phi) + pnonp * log(1 - phi))
                    },
                    cloglog = function(par, X, pnonp) {
                      theta1 <- par[1]
                      theta2 <- par[2]
                      eta <- theta1 + theta2 * X
                      phi <- 1 - exp(-exp(eta))
                      phi <- pmin(pmax(phi, 1e-8), 1 - 1e-8)
                      sum((1 - pnonp) * log(phi) + pnonp * log(1 - phi))
                    }
  )


  environment(Lik_fun) <- environment()
  if (is.null(theta0)) {
    invisible(capture.output(
      cure <- smcure(Surv(Y, Delta) ~ x,
                     cureform = ~x, data = data,
                     model = "ph", Var = FALSE, link = model)
    ))

    starting <- coefsmcure(cure)[c(1, 2)]
    starting <- starting * 1.1
  } else {
    starting <- theta0[1:2] * 1.1
  }

  mixture <- optim(par = starting, fn = Lik_fun, method = "BFGS", control = list(fnscale = -1, maxit = 25),
                   X = X, pnonp = pnonp)

  theta1hat <- mixture$par[1]
  theta2hat <- mixture$par[2]

  phihat <- switch(model,
                   logit = exp(theta1hat + theta2hat * X) / (1 + exp(theta1hat + theta2hat * X)),
                   probit = pnorm(theta1hat + theta2hat * X),
                   cloglog = 1 - exp(-exp(theta1hat + theta2hat * X)))

  ppar <- 1 - phihat
  TT <- n * sqrt(h) * mean((pnonp - ppar)^2)



  wwg <- weight_cpp(X, g)
  pnonp <- beran_cpp(as.matrix(data), wwg)
  pnonpp <- matrix(pnonp, n, n, byrow = TRUE)
  Shatg <- t(beranT_cpp(as.matrix(data), wwg))
  Shatg <- (Shatg - pnonpp) / (1 - pnonpp)
  Ghatg <- beranC_cpp(as.matrix(data), wwg)

  mixture <- optim(par = c(theta1hat, theta2hat), fn = Lik_fun, method = "BFGS", control = list(fnscale = -1, maxit = 15),
                   X = X, pnonp = pnonp)
  theta1hatg <- mixture$par[1]
  theta2hatg <- mixture$par[2]

  phihatg <- switch(model,
                    logit = exp(theta1hatg + theta2hatg * X) / (1 + exp(theta1hatg + theta2hatg * X)),
                    probit = pnorm(theta1hatg + theta2hatg * X),
                    cloglog = 1 - exp(-exp(theta1hatg + theta2hatg * X)))

  pparg <- 1 - phihatg

  TTstar <- rep(0, nsimb)
  countdiv <- 0
  b <- 1

  while ((b <= nsimb + countdiv)) {
    uu <- runif(n)
    Zstar <- (uu > pparg)
    Tstar <- rep(100000, n)

    uu1 <- matrix(runif(n), n, n, byrow = TRUE)
    uu1 = (uu1 > 1 - Shatg)
    uu1 <- apply(uu1, 2, sum) + 1
    uu1 <- uu1[Zstar == 1]
    Tstar[Zstar == 1] <- Y[uu1]

    uu2 <- matrix(runif(n), n, n, byrow = TRUE)
    uu2 <- (uu2 > Ghatg)
    uu2 <- apply(uu2, 2, sum) + 1
    Cstar <- Y[uu2]

    Ystar <- apply(cbind(Tstar,Cstar),1,min)
    Deltastar <- as.numeric(Tstar <= Cstar)

    datastar <- cbind(Ystar, Deltastar, X)
    datastar <- datastar[order(datastar[,1], 1 - datastar[,2]), ]
    X <- datastar[,3]

    ww = weight_cpp(X,h)
    pnonp = beran_cpp(as.matrix(datastar),ww)

    mixture <- optim(par = c(theta1hat, theta2hat), fn = Lik_fun, method = "BFGS", control = list(fnscale = -1, maxit = 15),
                     X = X, pnonp = pnonp)

    theta1star <- mixture$par[1]
    theta2star <- mixture$par[2]

    phistar <- switch(model,
                      logit = exp(theta1star + theta2star * X) / (1 + exp(theta1star + theta2star * X)),
                      probit = pnorm(theta1star + theta2star * X),
                      cloglog = 1 - exp(-exp(theta1star + theta2star * X)))

    ppar <- 1 - phistar
    TTstar[b - countdiv] <- n * sqrt(h) * mean((pnonp - ppar)^2)

    if (mixture$convergence == 1) countdiv <- countdiv + 1
    b <- b + 1
  }

  p.value <- mean(TTstar >= TT)

  out <- list(
    statistic = TT,
    p.value = p.value,
    bandwidth = h,
    model = model,
    B = nsimb
  )

  class(out) <- "goft"
  return(out)
}
