#' @export
p_estim <- function(data, return = NULL, risk = NULL, rf = 0, digits = 2) {

  if (xor(is.null(return), is.null(risk)) == FALSE)
    stop("Provide exactly one of return or risk.")

  mu  <- colMeans(data)
  S   <- cov(data)
  one <- rep(1, length(mu))

  invS <- chol2inv(chol(S))

  A <- drop(crossprod(one, invS %*% one))
  B <- drop(crossprod(one, invS %*% mu))
  C <- drop(crossprod(mu,  invS %*% mu))
  D <- A*C - B^2

  r <- if (!is.null(return)) {
    return
  } else {
    s2 <- risk^2
    (B + sqrt(max(B^2 - A*(C - D*s2), 0))) / A
  }

  w <- drop(invS %*% (((C - B*r)/D)*one + ((A*r - B)/D)*mu))
  names(w) <- colnames(data)

  R <- sum(w * mu)
  sd <- sqrt(drop(crossprod(w, S %*% w)))


  list(
    weights = round(w, digits),
    metrics = round(c(Return = R, Risk = sd, Sharpe = (R - rf)/sd), 6)
  )
}
