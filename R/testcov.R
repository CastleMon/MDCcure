#' @title Covariate Hypothesis Test of the Cure Probability based on Martingale Difference Correlation
#'
#' @description
#' Performs nonparametric hypothesis tests to evaluate the association between a covariate and the cure probability
#' in mixture cure models. Several test statistics are supported, including martingale difference correlation (MDC)-based tests
#' and an alternative GOFT test.
#'
#' @param x A numeric vector representing the covariate of interest.
#' @param time A numeric vector of observed survival times.
#' @param delta A binary vector indicating censoring status: \code{1} for event and \code{0} for censored.
#' @param h Bandwidth parameter for kernel smoothing. Either a positive numeric value, \code{NULL}, or the character string \code{"bootstrap"}.
#' If \code{NULL}, an optimal bandwidth is selected automatically. If \code{"bootstrap"}, the bandwidth is selected using the bootstrap method
#' proposed by López-Cheda et al. (2016).
#' @param method Character string specifying the test to perform. One of:
#' \itemize{
#'   \item \code{"MDCU"}: Martingale Difference Correlation with U-centering.
#'   \item \code{"MDCV"}: Martingale Difference Correlation with double-centering.
#'   \item \code{"FMDCU"}: Fast approximation of MDC with U-centering.
#'   \item \code{"GOFT"}: Goodness-of-fit test for the cure model.
#'   \item \code{"All"}: All of the above tests.
#' }
#' Default is \code{"FMDCU"}.
#' @param P Integer. Number of permutations or bootstrap replications used to compute the null distribution of the test statistic.
#' For methods \code{"MDCU"} or \code{"MDCV"}, this is the number of permutations. For the \code{"GOFT"} method, it is the number of bootstrap replications.
#' Defaults to 999.
#' @param parallel Logical. If \code{TRUE}, parallel computing is used to speed up computations. Default is \code{TRUE}.
#' @param ncores Integer. Number of cores to use for parallel computing. If \code{NULL}, it defaults to one less than the number of available cores.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{test_results}: A list with the results (e.g., test statistics and p-values) of the selected test(s).
#'   \item \code{nu_hat}: A numeric vector of estimated cure probabilities.
#' }
#'
#' @details
#' The function computes a statistic, based on the methodology proposed by Monroy-Castillo et al.,
#' to test whether a covariate \eqn{\boldsymbol{X}} has an effect on the cure probability.
#'
#' \deqn{
#'   \mathcal{H}_0 : \mathbb{E}(\nu | \boldsymbol{X}) \equiv 1 - p \quad \text{a.s.}
#'   \quad \text{vs} \quad
#'   \mathcal{H}_1 : \mathbb{E}(\nu | \boldsymbol{X}) \not\equiv 1 - p \quad \text{a.s.}
#' }
#'
#' The main problem is that the response variable (cure indicator \eqn{\nu}) is partially observed due to censoring.
#' This is addressed by estimating the cure indicator using the methodology of Amico et al. (2021).
#' We define \eqn{\tau = \sup_x \tau(x)}, with \eqn{\tau(x) = \inf\{t: S_0(t|x) = 0\}}.
#' We assume \eqn{\tau < \infty} and that follow-up is long enough so that \eqn{\tau < \tau_{G(x)}} for all \eqn{x}.
#' Therefore, individuals with censored observed times greater than \eqn{\tau} are considered cured (\eqn{\nu = 1}).
#'
#' Four tests are proposed: three are based on the martingale difference correlation (MDC).
#' For the MDCU and MDCV tests, the null distribution is approximated via a permutation procedure.
#' To provide a faster alternative, a chi-squared approximation is implemented for the MDCU test statistic (FMDCU).
#' Additionally, a modified version of the goodness-of-fit test proposed by Müller and Van Keilegom (2019) is included (GOFT).
#' The test statistic is given by:
#' \deqn{
#'    \widehat{\mathcal{T}}_n = nh^{1/2}\frac{1}{n}\sum_{i = 1}^{n}\left\{\hat{p}_h(X_i) - \hat{p}\right\}^2,
#' }
#' where \eqn{\hat{p}_h(X_i)} denotes the nonparametric estimator of the cure probability under the alternative hypothesis,
#' and \eqn{\hat{p}} denotes the nonparametric estimator of the cure probability under the null hypothesis.
#' The approximation of the critical value for the test is done using the bootstrap procedure given in Section 3 of Müller and Van Keilegom (2019).
#'
#' @references
#' Amico, M, Van Keilegom, I. & Han, B. (2021).
#' Assessing cure status prediction from survival data using receiver operating characteristic curves.
#' \emph{Biometrika}, 108(3), 727–740. \url{10.1093/biomet/asaa080}
#'
#' López-Cheda, A., Cao, R., Jácome, M. A., & Van Keilegom, I. (2016).
#' Nonparametric incidence estimation and bootstrap bandwidth selection in mixture cure models.
#' \emph{Computational Statistics & Data Analysis}, 100, 490–502. \url{10.1016/j.csda.2016.04.006}
#'
#' Müller, U.U, & Van Keilegom, I. (2019).
#' Goodness-of-fit tests for the cure rate in a mixture cure model.
#' \emph{Biometrika}, 106, 211-227. \url{10.1093/biomet/asy058}
#'
#' Shao, X., & Zhang, J. (2014). Martingale difference correlation and its use in high-dimensional variable screening.
#' \emph{Journal of the American Statistical Association}, 105, 144-165. \url{10.1080/01621459.2014.887012}
#'
#' @seealso
#'
#' \code{\link{mdc}}, \code{\link{mdd}}, \code{\link{mdc_test}}, \code{\link{testcov2}}
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
## Test of the significance of the covariate 'x'
#' testcov(x, t, d)
#'
#' @importFrom npcure probcure latency
#' @importFrom survival survfit Surv
#' @export


testcov <- function(x, time, delta, h = NULL, method = "FMDCU", P = 999, parallel = TRUE, ncores = -1) {
  data <- data.frame(x = x, Time = time, delta = delta)
  data <- na.omit(data)
  n <- length(data$x)
  pos <- which(data$delta == 0)

  cure <- min(survfit(Surv(Time, delta) ~ 1, data = data)$surv)
  tau <- max(data$Time[data$delta == 1])

  pos1 <- pos[which(data$Time[pos] <= tau)]
  if (length(pos1) == 0) {
    nu_hat <- rep(0, n)
    nu_hat[pos] <- 1
  } else {
    x0 <- data$x[pos1]
    s0 <- NULL

    if (is.character(h) && h == "bootstrap") {
      use_h <- FALSE
    } else {
      use_h <- TRUE
      h1 <- if (is.null(h)) {
        2 * (max(data$x) - min(data$x)) * n^(-1/5)
      } else {
        h
      }
    }

    s0 <- sapply(1:length(x0), function(j) {
      suppressWarnings({
        est <- latency(
          x, Time, delta, data,
          x0 = x0[order(x0)[j]],
          h = if (use_h) h1 else NULL,
          testimate = sort(data$Time)[which(order(data$Time) == pos1[order(x0)[j]])]
        )
        est$S
      })
    })

    s0 <- unlist(s0)

    nu_hat <- rep(0, n)
    nu_hat[pos] <- 1
    nu_hat[pos1[order(x0)]] <- cure / (cure + (1 - cure) * s0)
  }

  x <- as.matrix(x)
  result <- switch(
    method,
    "All" = {
      permMDC_U <- permutation_test_cpp_parallel(x, nu_hat, n_permutations = P, center = "U", parallel, ncores)

      permMDC_V <- permutation_test_cpp_parallel(x, nu_hat, n_permutations = P, center = "D", parallel, ncores)

      FMDCU_stat <- mdc_cpp(x, nu_hat, center = "U")
      FMDCU <- list(
        statistic = FMDCU_stat,
        p.value = 1 - pchisq(length(nu_hat) * FMDCU_stat + 1, df = 1)
      )

      goft <- GOFT_cpp(data, nsimb = P, h = h)

      list(MDCU = permMDC_U, MDCV = permMDC_V, FMDCU = FMDCU, GOFT = goft)
    },

    "MDCU" = permutation_test_cpp_parallel(x, nu_hat, n_permutations = P, center = "U", parallel, ncores),

    "MDCV" = permutation_test_cpp_parallel(x, nu_hat, n_permutations = P, center = "D", parallel, ncores),

    "FMDCU" = {
      FMDCU_stat <- mdc_cpp(x, nu_hat, center = "U")
      list(statistic = FMDCU_stat, p.value = 1 - pchisq(length(nu_hat) * FMDCU_stat + 1, df = 1))
    },
    "GOFT" = GOFT_cpp(data, nsimb = P, h = h),

    stop("Unknown method")
  )

  structure(
    list(result = result, nu_hat = nu_hat, method = method, P = P),
    class = "testcov"
  )

}
