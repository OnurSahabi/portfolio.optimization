#' @export
p_estim <- function(data, target_return, rf = 0) {

  mu <- colMeans(data)
  Sigma <- cov(data)

  n <- length(mu)
  ones <- rep(1, n)

  # KKT sistemi (eşitlik kısıtları)
  KKT <- rbind(
    cbind(2 * Sigma, ones, mu),
    cbind(t(ones), 0, 0),
    cbind(t(mu),   0, 0)
  )

  rhs <- c(rep(0, n), 1, target_return)

  sol <- solve(KKT, rhs)
  w <- sol[1:n]
  names(w) <- colnames(data)

  port_return <- sum(w * mu)
  port_risk   <- sqrt(as.numeric(t(w) %*% Sigma %*% w))
  port_sharpe <- (port_return - rf) / port_risk

  return(list(
    weights = w,
    metrics = c(
      Return = port_return,
      Risk   = port_risk,
      Sharpe = port_sharpe
    )
  ))
}
