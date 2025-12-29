#' @export
p_optim <- function(data, rf = 0, digits = 2) {
  
  mu <- colMeans(data)
  Sigma <- cov(data)
  n <- ncol(data)
  assets <- colnames(data)
  
  # ---------- yardımcı fonksiyonlar ----------
  port_return <- function(w) sum(w * mu)
  port_risk   <- function(w) sqrt(as.numeric(t(w) %*% Sigma %*% w))
  port_sharpe <- function(w) (port_return(w) - rf) / port_risk(w)
  
  # ---------- Max Sharpe ----------
  sharpe_obj <- function(w) {
    w <- w / sum(w)
    -( (sum(w * mu) - rf) / sqrt(t(w) %*% Sigma %*% w) )
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
  names(w_var) <- assets
  
  # ---------- results_df ----------
  results_df <- data.frame(
    Type   = c("Max Sharpe", "Min Risk"),
    Sharpe = c(port_sharpe(w_sharpe), port_sharpe(w_var)),
    Return = c(port_return(w_sharpe), port_return(w_var)),
    Risk   = c(port_risk(w_sharpe),   port_risk(w_var))
  )
  
  # ---------- weights ----------
  weights <- rbind(
    Max_Sharpe = w_sharpe,
    Min_Risk   = w_var
  )
  weights <- round(weights, digits)
  
  # ---------- yazdırma ----------
  return(list(
    results = results_df,
    weights = weights
  ))
}



