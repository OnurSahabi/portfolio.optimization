#' @export
p_estim <- function(data, risk = NULL, return = NULL, rf = 0, digits = 2) {

  if (!requireNamespace("CVXR", quietly = TRUE))
    stop("CVXR kurulu değil")

  if (xor(is.null(risk), is.null(return)) == FALSE)
    stop("Sadece risk veya sadece return vermelisin")

  mu <- colMeans(data)
  S  <- cov(data)
  L  <- chol(S)
  n  <- length(mu)

  w <- CVXR::Variable(n)

  obj <- if (is.null(return)) {
    CVXR::Maximize(t(mu - rf) %*% w)
  } else {
    CVXR::Minimize(CVXR::norm2(L %*% w))
  }

  cons <- list(
    w >= 0,
    sum(w) == 1
  )

  cons <- if (is.null(return)) {
    c(cons, list(CVXR::norm2(L %*% w) <= risk))
  } else {
    c(cons, list(t(mu) %*% w == return))
  }

  sol <- CVXR::solve(CVXR::Problem(obj, cons))

  w_hat <- as.numeric(sol$getValue(w))
  names(w_hat) <- colnames(data)

  Return  <- sum(mu * w_hat)
  Risk <- sqrt(drop(t(w_hat) %*% S %*% w_hat))
  Sharpe = (Return - rf) / Risk

  list(
    weights = round(w_hat, digits),
    metrics = round(c(Return, Risk, Sharpe), 6)
  )
}
