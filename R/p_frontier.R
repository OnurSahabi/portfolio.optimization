#' @export
p_frontier <- function(data, step = 0.01, plot = FALSE) {

  R  <- na.omit(as.matrix(data))
  mu <- colMeans(R)

  S  <- cov(R)
  S  <- (S + t(S)) / 2

  n <- length(mu)

  Dmat <- 2 * S
  dvec <- rep(0, n)

  # --- GMV ---
  w_gmv <- quadprog::solve.QP(
    Dmat, dvec,
    cbind(rep(1,n), diag(n)),
    c(1, rep(0,n)),
    meq = 1
  )$solution

  ret_gmv <- sum(w_gmv * mu)
  ret_max <- max(mu)

  targets <- seq(ret_gmv, ret_max, length.out = round(1/step))
  Amat    <- cbind(mu, rep(1,n), diag(n))

  frontier <- sapply(targets, function(r) {

    w <- quadprog::solve.QP(
      Dmat, dvec,
      Amat,
      c(r,1,rep(0,n)),
      meq = 2
    )$solution

    c(
      Return = sum(w * mu),
      Risk   = sqrt(drop(t(w) %*% S %*% w))
    )
  })

  result <- list(
    Risk   = frontier["Risk",],
    Return = frontier["Return",]
  )

  if (isTRUE(plot)) {
    graphics::plot(result$Risk, result$Return,
                   type="l", lwd=2, col="darkblue",
                   xlab=expression(Risk~(sigma)),
                   ylab=expression(Expected~Return~(mu)),
                   main="Markowitz Efficient Frontier")
    graphics::grid()
  }

  return(result)
}
