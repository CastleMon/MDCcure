weight = function(xx,h)
{
  n = length(xx)
  xx = matrix(xx,n,n,byrow=TRUE)
  z = (xx-t(xx))/h
  z = replace(z,abs(z)>1,-1)
  #w = 15/16*(1-z^2)^2
  w = 3/4*(1-z^2)
  sumw = apply(w,1,sum)
  w = w/matrix(sumw,n,n,byrow=FALSE)
  return(w)
}

beran = function(data,w)
{
  n = nrow(data)
  k = cbind(rep(0,n),w[,1:(n-1)])
  cumsumk = t(apply(k,1,cumsum))
  cumw = 1-cumsumk
  w = replace(w,cumw==0,0)  # CHANGE 0 TO 1 ?
  cumw = replace(cumw,cumw==0,1)
  q = 1 - w/cumw
  delta = matrix(data[,2],n,n,byrow=TRUE)
  q = replace(q,delta==0,1)
  q = apply(q,1,prod)
  return(q)
}


beranT = function(data,w)
{
  n = nrow(data)
  k = cbind(rep(0,n),w[,1:(n-1)])
  cumsumk = t(apply(k,1,cumsum))
  cumw = 1-cumsumk
  w = replace(w,cumw==0,0)   # CHANGE 0 TO 1
  cumw = replace(cumw,cumw==0,1)
  q = 1 - w/cumw
  delta = matrix(data[,2],n,n,byrow=TRUE)
  q = replace(q,delta==0,1)
  q = apply(q,1,cumprod)
  return(q)
}


beranC = function(data,w)
{
  n = nrow(data)
  k = cbind(rep(0,n),w[,1:(n-1)])
  cumsumk = t(apply(k,1,cumsum))
  cumw = 1-cumsumk
  w = replace(w,cumw==0,0)   # CHANGE 0 TO 1
  cumw = replace(cumw,cumw==0,1)
  q = 1 - w/cumw
  delta = matrix(data[,2],n,n,byrow=TRUE)
  q = replace(q,delta==1,1)
  q = 1 - apply(q,1,cumprod)
  q[n,] = rep(1,n)
  return(q)
}

Lik_logit = function(par) {
  theta1 = par[1]
  theta2 = par[2]
  phi = exp(theta1+theta2*X)
  phi = phi/(1+phi)
  likelihood = sum((1-pnonp)*log(phi)+pnonp*log(1-phi))
  return(likelihood)
}

Lik_probit <- function(par) {
  theta1 <- par[1]
  theta2 <- par[2]
  phi <- pnorm(theta1 + theta2 * X)
  likelihood <- sum((1 - pnonp) * log(phi) + pnonp * log(1 - phi))
  return(likelihood)
}

Lik_cloglog <- function(par) {
  theta1 <- par[1]
  theta2 <- par[2]
  eta <- theta1 + theta2 * X
  phi <- 1 - exp(-exp(eta))
  likelihood <- sum((1 - pnonp) * log(phi) + pnonp * log(1 - phi))
  return(likelihood)
}

GOFT.cure <- function(Data1, nsimb, h = NULL){
  boot <- 1
  if (any(is.na(Data1[,1]))) {
    invalid_indices <- which(is.na(Data1[,1]))
    Data <- Data1[-c(invalid_indices), ]
  } else {
    Data <- Data1
  }

  Y <- Data[,2]
  X <- Data[,1]
  Delta <- Data[,3]
  n <- length(X)
  data <- cbind(Y, Delta, X)

  data = data[order(data[,1],1-data[,2]),]
  Y = data[,1]
  Delta = data[,2]
  X = data[,3]
  Data <- as.data.frame(data)

  if(is.null(h)){
    ch <- 2*(max(X) - min(X))
    h <- ch*n^(-1/5)
  }else{
    h <- h
  }

  ######
  # Calculation of the nonparametric estimator
  ######

  ww <- weight(data[,3],h)
  pnonp <- beran(data,ww)

  ######
  # Calculation of the parametric estimator
  ######

  cure <- min(survfit(Surv(Y, Delta) ~ 1, data = Data)$surv)

  ppar <- cure

  ######
  # Calculation of the gof test statistic
  ######

  TT = n*sqrt(h)*mean((pnonp-ppar)^2)

  if (boot==1){
    ######
    # Calculation of the bootstrap critical value
    ######

    g = h*n^0.09
    wwg = weight(X,g)
    pnonp = beran(data,wwg)
    pnonpp = matrix(pnonp,n,n,byrow=TRUE)
    Shatg = beranT(data,wwg)
    Shatg = (Shatg-pnonpp)/(1-pnonpp)
    Ghatg = beranC(data,wwg)

    pparg = cure

    TTstar = rep(0,nsimb)
    b = 1
    stop = 0

    while ((b <= nsimb)*(stop==0)){
      uu = stats::runif(n)
      Zstar = (uu>pparg)
      Tstar = rep(100000,n)
      uu1 = runif(n)
      uu1 = matrix(uu1,n,n,byrow=TRUE)
      uu1 = (uu1>1-Shatg)
      uu1 = apply(uu1,2,sum)+1
      uu1 = uu1[Zstar==1]
      Tstar[Zstar==1] = Y[uu1]

      uu2 = runif(n)
      uu2 = matrix(uu2,n,n,byrow=TRUE)
      uu2 = (uu2>Ghatg)
      uu2 = apply(uu2,2,sum)+1
      Cstar = Y[uu2]

      Ystar = apply(cbind(Tstar,Cstar),1,min)
      Deltastar = (Tstar <= Cstar)

      datastar = cbind(Ystar,Deltastar,X)
      datastar = datastar[order(datastar[,1],1-datastar[,2]),]
      Ystar = datastar[,1]
      Deltastar = datastar[,2]
      X = datastar[,3]

      Datastar <- as.data.frame(datastar)

      ww = weight(X,h)
      pnonp = beran(datastar,ww)

      ppar <- min(survfit(Surv(Ystar, Deltastar) ~ 1, data = Datastar)$surv)


      TTstar[b] = n*sqrt(h)*mean((pnonp-ppar)^2)
      b = b + 1

    }
    TTstar = sort(TTstar)
    p.value <- mean(TTstar >= TT)
  }
  return(list(statistic = TT, p.valUe = p.value))
}
