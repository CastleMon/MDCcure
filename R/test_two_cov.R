#' Hypothesis test for association between covariate and cure indicator adjusted by a second covariate
#'
#' Performs a permutation-based test assessing the association between a primary covariate (`x`) and the cure indicator, while adjusting for a secondary covariate (`z`).
#' The test calculates the p-value via permutation using the partial martingale difference correlation.
#'
#' @param x Numeric vector. The primary covariate whose association with the latent cure indicator is tested.
#' @param time Numeric vector. Observed survival or censoring times.
#' @param z Numeric vector. Secondary covariate for adjustment.
#' @param delta Numeric vector. Censoring indicator (1 indicates event occurred, 0 indicates censored).
#' @param P Integer. Number of permutations used to compute the permutation p-value. Default is 999.
#' @param H Optional numeric. Bandwidth parameter (currently unused, reserved for future extensions).
#'
#' @return List with components:
#' \describe{
#'   \item{statistic}{Numeric. The test statistic value.}
#'   \item{p.value}{Numeric. The permutation p-value assessing the null hypothesis of no association between `x` and the latent cure indicator, adjusting for `z`.}
#' }
#'
#' @details
#' In order to test if the cure rate depends on the covariate \eqn{\boldsymbol{X}} given it depends on the covariate \eqn{\boldsymbol{Z}}. The hypotheses are
#' \deqn{
#'   \mathcal{H}_0 : \mathbb{E}(\nu | \boldsymbol{X}) \equiv 1 - p(\boldsymbol{X}) \quad \text{a.s.}
#'   \quad \text{vs} \quad
#'   \mathcal{H}_1 : \mathbb{E}(\nu | \boldsymbol{X}) \not\equiv 1 - p(\boldsymbol{X}) \quad \text{a.s.}
#' }
#' The proxy of the cure rate under the null hypothesis \eqn{\mathcal{H}_0} is obtained by:
#' \deqn{
#'  \mathbb{I}(T > \tau) + (1-\delta)\mathbb{I}(T \leq \tau) \, \frac{1 - p(\boldsymbol{Z})}{1 - p(\boldsymbol{Z}) + p(\boldsymbol{Z})S_0(T|\boldsymbol{X,Z})}.
#' }
#'
#' The statistic for testing the covariate hypothesis is based on partial martingale difference correlation and it is given by:
#'
#' \deqn{
#'  \text{pMDC}_n(\hat{\nu}_{\boldsymbol{H}}|\boldsymbol{X,Z})^2.
#' }
#'
#' The null distribution is approximated using a permutation test.
#'
#' @references
#'
#' Park, T., Saho, X. & Yao, S. (2015). Partial martingale difference correlation.
#' \emph{Electronic Journal of Statistics}, 9, 1492–1517. \url{10.1214/15-EJS1047}
#'
#' @seealso \code{\link{pmdc}} for the partial martingale difference correlation, \code{\link{pmdd}} for the partial martingale difference divergence,
#'  \code{\link{testcov}} for the test for one covariate.
#'
#'
#' @importFrom npcure probcure
#' @export



testcov2 <- function(x, time, z, delta, P = 999, H = NULL){
  call <- match.call()
  X_name <- deparse(call$x)
  Z_name <- deparse(call$z)

  n <- length(x)
  data <- data.frame(x = x, z = z, Time = time, delta = delta)
  data <- data[order(data$Time, 1 - data$delta), ]

  tau <- max(data$Time[data$delta == 1])
  pos <- which(data$delta == 0)
  pos1 <- pos[data$Time[pos] <= tau]

  if (length(pos1) == 0) {
    nu_hat <- rep(0, n)
    nu_hat[pos] <- 1
  } else {
    S0 <- sapply(pos1, function(idx) {
      x0 <- c(data$x[idx], data$z[idx])
      t <- data$Time[idx]
      latency_estimator_multivariate_cpp(t, x0, data$Time, as.matrix(data[, c("x", "z")]), data$delta, H)
    })

    if (any(is.null(S0))) return(NA)

    Ch_value <- 2 * (max(data$z) - min(data$z))
    h0 <- Ch_value * n^(-1/5)
    h.1 <- rep(h0, length(pos1))
    cure <- probcure(z, Time, delta, data[, c("z", "Time", "delta")], x0 = data$z[pos1], h = h.1)$q

    if (is.null(cure) || length(cure) != length(pos1)) return(NA)
    cure <- cure[order(data$z[pos1])]

    nu_hat <- rep(0, n)
    nu_hat[pos] <- 1
    nu_hat[pos1] <- cure / (cure + (1 - cure) * S0)
  }

  permMDC_U <- permutation_test_pmdc.cpp(as.matrix(data$x), as.matrix(nu_hat), as.matrix(data$z), n_permutations = P)
  res <- list(
    test = list(
      statistic = permMDC_U$statistic,
      p.value = permMDC_U$p.value,
      B = P
    ),
    nu_hat = nu_hat,
    varnames = list(x = X_name, z = Z_name)
  )
  class(res) <- "testcov2"
  return(res)
}
