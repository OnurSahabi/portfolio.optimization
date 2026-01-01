#' @export
p_optim <- function(data, rf = 0, digits = 2) {

  mu <- colMeans(data)
  Sigma <- cov(data)
  n <- ncol(data)
  assets <- colnames(data)

  # ---------- yardımcı ----------
  port_return <- function(w) sum(w * mu)
  port_risk   <- function(w) sqrt(as.numeric(t(w) %*% Sigma %*% w))
  port_sharpe <- function(w) (port_return(w) - rf) / port_risk(w)

  # ---------- Max Sharpe (AYNI) ----------
  sharpe_obj <- function(w) {
    w <- w / sum(w)
    -port_sharpe(w)
  }

  res_sharpe <- optim(
    par    = rep(1, n),
    fn     = sharpe_obj,
    method = "L-BFGS-B",
    lower  = rep(0, n)
  )

  w_sharpe <- res_sharpe$par / sum(res_sharpe$par)
  names(w_sharpe) <- assets

  # ---------- Min Risk ----------
  if (rf == 0) {

    # ---- eski davranış ----
    var_obj <- function(w) {
      w <- w / sum(w)
      as.numeric(t(w) %*% Sigma %*% w)
    }

    res_var <- optim(
      par    = rep(1, n),
      fn     = var_obj,
      method = "L-BFGS-B",
      lower  = rep(0, n)
    )

    w_var <- res_var$par / sum(res_var$par)

  } else {

    if (!requireNamespace("quadprog", quietly = TRUE))
      install.packages("quadprog")
    library(quadprog)

    Dmat <- 2 * Sigma
    dvec <- rep(0, n)

    # Amat^T w >= bvec
    Amat <- cbind(
      rep(1, n),        # sum(w) >= 1
      -rep(1, n),       # -sum(w) >= -1
      mu,               # mu'w >= rf
      diag(n)           # w >= 0
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
  stats_df <- data.frame(
    Type   = c("Max Sharpe", "Min Risk"),
    Sharpe = c(port_sharpe(w_sharpe), port_sharpe(w_var)),
    Return = c(port_return(w_sharpe), port_return(w_var)),
    Risk   = c(port_risk(w_sharpe),   port_risk(w_var))
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
    stats = stats_df,
    weights = weights
  ))
}



