GOFT_cpp <- function(Data1, nsimb, h = NULL) {
  if (any(is.na(Data1[, 1]))) {
    Data <- Data1[!is.na(Data1[, 1]), ]
  } else {
    Data <- Data1
  }

  Y <- Data[, 2]
  X <- Data[, 1]
  Delta <- Data[, 3]
  n <- length(X)

  data <- cbind(Y, Delta, X)

  data <- data[order(data[, 1], 1 - data[, 2]), ]
  Y <- data[, 1]
  Delta <- data[, 2]
  X <- data[, 3]
  Data <- as.data.frame(data)

  if (is.null(h)) {
    ch <- (max(X) - min(X))
    h <- 2 * ch * n^(-1/5)
  }

  data_mat <- as.matrix(data)
  ww <- weight_cpp(data_mat[, 3], h)
  pnonp <- beran_cpp(data_mat, ww)

  cure <- min(survfit(Surv(Y, Delta) ~ 1, data = Data)$surv)
  ppar <- cure

  TT <- n * sqrt(h) * mean((pnonp - ppar)^2)

  g <- h * n^0.09
  wwg <- weight_cpp(X, g)
  pnonp_g <- beran_cpp(data_mat, wwg)
  pnonpp <- matrix(pnonp_g, n, n, byrow = TRUE)

  Shatg <- t(beranT_cpp(data_mat, wwg))
  Shatg <- (Shatg - pnonpp) / (1 - pnonpp)
  Ghatg <- beranC_cpp(data_mat, wwg)

  pparg <- cure
  TTstar <- numeric(nsimb)

  set.seed(123)

  for (b in 1:nsimb) {
    uu <- runif(n)
    Zstar <- (uu > pparg)

    Tstar <- rep(1e5, n)
    uu1 <- stats::runif(n)
    uu1 <- matrix(uu1, n, n, byrow = TRUE)
    uu1 <- (uu1 > 1 - Shatg)
    uu1 <- apply(uu1, 2, sum) + 1
    uu1 <- uu1[Zstar == 1]
    Tstar[Zstar == 1] <- Y[uu1]

    uu2 <- runif(n)
    uu2 <- matrix(uu2, n, n, byrow = TRUE)
    uu2 <- (uu2 > Ghatg)
    uu2 <- apply(uu2, 2, sum) + 1
    Cstar <- Y[uu2]

    Ystar <- pmin(Tstar, Cstar)
    Deltastar <- (Tstar <= Cstar) * 1

    datastar <- cbind(Ystar, Deltastar, X)
    datastar <- datastar[order(datastar[, 1], 1 - datastar[, 2]), ]
    Ystar <- datastar[, 1]
    Deltastar <- datastar[, 2]
    X <- datastar[, 3]

    Datastar <- as.data.frame(datastar)

    ww <- weight_cpp(X, h)
    pnonp_star <- beran_cpp(as.matrix(datastar), ww)

    ppar_star <- min(survfit(Surv(Ystar, Deltastar) ~ 1, data = Datastar)$surv)

    TTstar[b] <- n * sqrt(h) * mean((pnonp_star - ppar_star)^2)
  }

  TTstar <- sort(TTstar)
  p.value <- mean(TTstar >= TT)

  return(list(statistic = TT, p.value = p.value))
}
