#' @export
p_optim <- function(data, rf = 0, digits = 2) {

  # ---- Input checks ----
  if (!is.matrix(data) && !is.data.frame(data))
    stop("data must be matrix or data.frame")

  if (ncol(data) < 2)
    stop("At least two assets required")

  data <- as.matrix(data)

  mu    <- colMeans(data)
  Sigma <- cov(data)
  Sigma <- (Sigma + t(Sigma)) / 2
  Sigma <- Sigma + diag(1e-10, ncol(Sigma))

  n      <- ncol(data)
  assets <- colnames(data)

  # ---- Portfolio functions ----
  port_return <- function(w) sum(w * mu)
  port_risk   <- function(w) sqrt(drop(t(w) %*% Sigma %*% w))
  port_sharpe <- function(w) port_return(w) / port_risk(w)

  # ---- Max Sharpe ----
  sharpe_obj <- function(w) {
    w <- w / sum(w)
    -port_sharpe(w)
  }

  sharpe <- optim(
    par    = rep(1/n, n),
    fn     = sharpe_obj,
    method = "L-BFGS-B",
    lower  = rep(0, n)
  )

  w_sharpe <- sharpe$par / sum(sharpe$par)
  names(w_sharpe) <- assets

  # ---- Min Variance ----
  Dmat <- 2 * Sigma
  dvec <- rep(0, n)

  Amat <- cbind(
    rep(1, n),
    mu,
    diag(n)
  )

  bvec <- c(
    1,
    rf,
    rep(0, n)
  )

  qp <- quadprog::solve.QP(
    Dmat = Dmat,
    dvec = dvec,
    Amat = Amat,
    bvec = bvec,
    meq  = 1
  )

  w_var <- qp$solution
  names(w_var) <- assets

  # ---- Results ----
  stats <- data.frame(
    Type   = c("Max Sharpe", "Min Risk"),
    Sharpe = round(c(port_sharpe(w_sharpe),
                     port_sharpe(w_var)), 6),
    Return = round(c(port_return(w_sharpe),
                     port_return(w_var)), 6),
    Risk   = round(c(port_risk(w_sharpe),
                     port_risk(w_var)), 6)
  )

  weights <- round(
    rbind(
      Max_Sharpe = w_sharpe,
      Min_Risk   = w_var
    ),
    digits
  )

  return(list(
    stats   = stats,
    weights = weights
  ))
}
