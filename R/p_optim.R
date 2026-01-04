#' @export
p_optim <-function(data, rf = 0, digits = 2) {

  mu <- colMeans(data)
  Sigma <- cov(data)
  n <- ncol(data)
  assets <- colnames(data)

  port_return <- function(w) sum(w * mu)
  port_risk   <- function(w) sqrt(as.numeric(t(w) %*% Sigma %*% w))
  port_sharpe <- function(w) port_return(w) / port_risk(w)

  # ---------- Max Sharpe ----------
  sharpe_obj <- function(w) {
    w <- w / sum(w)
    -port_sharpe(w)
  }

  sharpe <- optim(
    par    = rep(1, n),
    fn     = sharpe_obj,
    method = "L-BFGS-B",
    lower  = rep(0, n)
  )

  w_sharpe <- sharpe$par / sum(sharpe$par)
  names(w_sharpe) <- assets

  # ---------- Min Risk ----------
  if (rf == 0) {

    var_obj <- function(w) {
      w <- w / sum(w)
      as.numeric(t(w) %*% Sigma %*% w)
    }

    var <- optim(
      par    = rep(1, n),
      fn     = var_obj,
      method = "L-BFGS-B",
      lower  = rep(0, n)
    )

    w_var <- var$par / sum(var$par)

  } else {

    if (!requireNamespace("quadprog", quietly = TRUE))
      install.packages("quadprog")
    library(quadprog)

    Dmat <- 2 * Sigma
    dvec <- rep(0, n)

    Amat <- cbind(
      rep(1, n),
      -rep(1, n),
      mu,
      diag(n)
    )

    bvec <- c(1, -1, rf, rep(0, n))

    qp <- solve.QP(
      Dmat = Dmat,
      dvec = dvec,
      Amat = Amat,
      bvec = bvec,
      meq  = 0
    )

    w_var <- qp$solution
    w_var[w_var < 0] <- 0
    w_var <- w_var / sum(w_var)
  }

  names(w_var) <- assets

  # ---------- results ----------
  stats <- data.frame(
    Type   = c("Max Sharpe", "Min Risk"),
    Sharpe = round(c(port_sharpe(w_sharpe), port_sharpe(w_var)), 6),
    Return = round(c(port_return(w_sharpe), port_return(w_var)), 6),
    Risk   = round(c(port_risk(w_sharpe),   port_risk(w_var)), 6)
  )

  # ---------- weights ----------
  weights <- round(
    rbind(
      Max_Sharpe = w_sharpe,
      Min_Risk   = w_var
    ),
    digits
  )

  return(list(
    stats = stats,
    weights = weights
  ))
}



