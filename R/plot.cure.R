#' Plot Cure Probability: A Comparison of Nonparametric and Parametric Estimation
#'
#' This function generates a plot comparing nonparametric and parametric estimations of cure probability in a univariate setting.
#' The nonparametric estimate is displayed with 95% confidence bands, while the parametric estimate is based on a logit, probit or
#' complementary log-log link.
#' An optional covariate density curve can be added as a secondary axis.
#'
#' @param x A numeric vector containing the covariate values.
#' @param time A numeric vector representing the observed survival times.
#' @param delta A binary vector indicating the event status (1 = event, 0 = censored).
#' @param main.title Character string for the main title of the plot. If \code{NULL}, a default is used.
#' @param title.x Character string for the x-axis label. If \code{NULL}, a default is used.
#' @param model A character string indicating the assumed model. Options include \code{"logit"}, \code{"probit"}, and \code{"cloglog"}. Defaults to \code{"logit"}.
#' @param theta A numeric vector of length 2, specifying the coefficients for the logistic model to generate the parametric estimate.
#' @param legend.pos A character string indicating the position of the legend. Options include \code{"bottom"}, \code{"top"}, \code{"left"}, \code{"right"}, \code{"none"}, etc.
#' @param density Logical; if \code{TRUE}, adds a secondary y-axis with the covariate density curve.
#' @param hsmooth Numeric. Smoothing bandwidth parameter (h) for the cure probability estimator.
#' @param npoints Integer. Number of points at which the estimator is evaluated over the covariate range.
#'
#' @details The function estimates the cure probability nonparametrically using the \code{\link{probcure}} function
#' and overlays it with a parametric estimate obtained from a logistic regression model.
#' Confidence intervals (95%) are included for the nonparametric estimate. Optionally,
#' the density of the covariate can be shown as a shaded area with a secondary y-axis.
#'
#' @import ggplot2
#' @import ggtext
#' @import future
#' @import future.apply
#' @import npcure
#' @importFrom stats density quantile
#' @importFrom gridExtra grid.arrange
#' @importFrom grid unit
#' @importFrom utils head
#' @importFrom smcure smcure coefsmcure
#' @importFrom stats na.omit
#'
#' @return A ggplot object representing the cure probability plot.
#'
#' @seealso \link[npcure:probcure]{probcure}
#'
#' @examples
#' \dontrun{
#'   set.seed(123)
#'   x <- runif(200)
#'   time <- rexp(200)
#'   delta <- rbinom(200, 1, 0.7)
#'   theta <- c(0.5, -1)
#'   plot.cure.uni(x = x, time = time, delta = delta,
#'                 main.title = "Estimated Cure Probability",
#'                 theta = theta,
#'                 legend.pos = "bottom",
#'                 density = FALSE)
#' }
#'
#' @export

plotCure <- function(x, time, delta,
                          main.title = NULL,
                          title.x = NULL,
                          model = "logit",
                          theta = NULL,
                          legend.pos = "bottom",
                          density = TRUE,
                          hsmooth = 10,
                          npoints = 100) {

  if (is.null(main.title)) main.title <- "Estimated Cure Probability"
  if (is.null(title.x)) title.x <- "Covariate"
  title.y <- "Cure Probability"

  Data_clean <- data.frame(covariate = x, Time = time, status = delta)

  x0 <- seq(quantile(Data_clean$covariate, 0.05),
            quantile(Data_clean$covariate, 0.95),
            length.out = npoints)


  suppressWarnings(plan(multisession, workers = parallel::detectCores() - 1))

  q.covariate <- suppressWarnings(
    future_lapply(x0, function(xi) {
      tryCatch({
        probcure(covariate, Time, status, Data_clean,
                 x0 = xi, conflevel = 0.95,
                 bootpars = controlpars(hsmooth = hsmooth))
      }, error = function(e) NULL)
    })
  )

  q.covariate <- q.covariate[!sapply(q.covariate, function(x) is.null(x) || length(x) == 0)]

  covariate_data <- do.call(rbind, lapply(q.covariate, function(res) {
    data.frame(
      x0 = res$x0,
      q = res$q,
      lower = res$conf$lower,
      upper = res$conf$upper
    )
  }))

  if (is.null(theta)) {
    invisible(capture.output(
      cure <- smcure(Surv(Time, status) ~ covariate,
                     cureform = ~covariate, data = Data_clean,
                     model = "ph", Var = FALSE, link = model)
    ))

    theta <- coefsmcure(cure)[c(1, 2)]
  } else {
    theta <- theta[1:2]
  }

  eta <- theta[1] + theta[2] * covariate_data$x0
  q.smc <- switch(tolower(model),
                  "logit" = 1 / (1 + exp(eta)),
                  "probit" = pnorm(-eta),
                  "cloglog" = exp(-exp(eta)),
                  stop("Invalid model. Choose one of 'logit', 'probit', or 'cloglog'.")
  )

  covariate_data$q.sm <- q.smc

  if (density) {
    d.covariate <- density(Data_clean$covariate)
    density_data <- data.frame(x = d.covariate$x, y = d.covariate$y)
    grob_density <- ggplotGrob(
      ggplot(density_data, aes(x = x, y = y)) +
        geom_line(color = "darkslategrey") +
        geom_area(fill = "darkslategrey", alpha = 0.12) +
        theme_void()
    )
  }

  plot_cure <- ggplot(covariate_data, aes(x = x0)) +
    geom_line(aes(y = lower, linetype = "95% CI (Nonparametric)", color = "95% CI (Nonparametric)")) +
    geom_line(aes(y = upper, linetype = "95% CI (Nonparametric)", color = "95% CI (Nonparametric)")) +
    geom_line(aes(y = q, color = "Nonparametric", linetype = "Nonparametric"), size = 0.9) +
    geom_line(aes(y = q.sm, color = "Parametric", linetype = "Parametric"), size = 0.9) +
    scale_color_manual(
      name = NULL,
      values = c("95% CI (Nonparametric)" = "dodgerblue4",
                 "Nonparametric" = "dodgerblue4",
                 "Parametric" = "darkred")
    ) +
    scale_linetype_manual(
      name = NULL,
      values = c("95% CI (Nonparametric)" = "dashed",
                 "Nonparametric" = "solid",
                 "Parametric" = "solid")
    ) +
    labs(
      title = main.title,
      x = title.x,
      y = title.y
    ) +
    theme_bw() +
    theme(
      legend.position = legend.pos,
      legend.text = element_text(size = 13),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 13),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
    )

  if (density) {
    plot_cure <- plot_cure +
      annotation_custom(grob_density, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
      scale_y_continuous(
        name = title.y,
        limits = c(0, 1),
        sec.axis = sec_axis(
          transform = ~ . * (max(density_data$y) / 1),
          name = "Density"
        )
      ) +
      theme(
        axis.title.y.right = element_text(color = "darkslategrey"),
        axis.text.y.right = element_text(color = "darkslategrey"),
        axis.line.y.right = element_line(color = "darkslategrey"),
        axis.title.y = element_text(color = "dodgerblue4"),
        axis.text.y = element_text(color = "dodgerblue4"),
        axis.line.y = element_line(color = "dodgerblue4")
      )
  } else {
    plot_cure <- plot_cure +
      scale_y_continuous(name = title.y, limits = c(0, 1)) +
      theme(
        axis.title.y = element_text(color = "dodgerblue4"),
        axis.text.y = element_text(color = "dodgerblue4"),
        axis.line.y = element_line(color = "dodgerblue4")
      )
  }

  print(plot_cure)
}
