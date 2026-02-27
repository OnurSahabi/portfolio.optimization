#' Beta Decomposition and Treynor Ratio (Full and Tail-Based)
#'
#' Computes CAPM beta coefficients for multiple assets relative to a market index.
#' In addition to full-sample beta, the function estimates lower-tail,
#' upper-tail, and middle-region betas based on quantile thresholds.
#'
#' Optionally, portfolio beta and Treynor ratio can be computed
#' if portfolio weights are provided.
#'
#' @param data_var A numeric matrix or data.frame containing asset returns.
#' @param data_index A numeric matrix or data.frame containing index returns.
#' @param var_cols Column indices or names of asset return variables.
#' @param index_col Column index or name of the market index variable.
#' @param weights Optional numeric vector of portfolio weights.
#'   Must have same length as number of assets.
#' @param tau Numeric. Tail probability level for quantile thresholds
#'   (default = 0.05).
#' @param rf Numeric. Risk-free rate used in Treynor ratio (default = 0).
#'
#' @return A data.frame containing:
#' \describe{
#'   \item{beta}{Full-sample CAPM beta}
#'   \item{beta_mid}{Beta estimated in the middle quantile region}
#'   \item{beta_low}{Lower-tail beta}
#'   \item{beta_high}{Upper-tail beta}
#'   \item{treynor}{Treynor ratio}
#' }
#'
#' @details
#' Tail betas are estimated using conditional subsamples defined by
#' the lower and upper quantiles of the market index returns.
#' The Treynor ratio is computed as:
#' \deqn{(E[R_i] - r_f) / \beta_i}
#'
#' @examples
#' set.seed(1)
#' R_assets <- matrix(rnorm(300), ncol = 3)
#' R_index  <- matrix(rnorm(100), ncol = 1)
#'
#' colnames(R_assets) <- c("A", "B", "C")
#' colnames(R_index)  <- "MKT"
#'
#' p_beta(R_assets, R_index,
#'        var_cols = 1:3,
#'        index_col = 1)
#'
#' @importFrom stats cov var quantile colMeans
#' @export
p_beta <- function(data_var, data_index, var_cols, index_col,
                   weights = NULL, tau = 0.05, rf = 0) {

  index_vec <- as.vector(data_index[, index_col])
  var_mat <- data_var[, var_cols]

  # Full sample beta
  beta <- cov(var_mat, index_vec)[, 1] / var(index_vec)

  # Quantile thresholds
  q_low  <- quantile(index_vec, tau)
  q_high <- quantile(index_vec, 1 - tau)

  # Lower tail
  idx_low <- index_vec <= q_low
  beta_low <- cov(
    var_mat[idx_low, , drop = FALSE],
    index_vec[idx_low]
  )[, 1] / var(index_vec[idx_low])

  # Upper tail
  idx_high <- index_vec >= q_high
  beta_high <- cov(
    var_mat[idx_high, , drop = FALSE],
    index_vec[idx_high]
  )[, 1] / var(index_vec[idx_high])

  # Middle region
  idx_mid <- index_vec > q_low & index_vec < q_high
  beta_mid <- cov(
    var_mat[idx_mid, , drop = FALSE],
    index_vec[idx_mid]
  )[, 1] / var(index_vec[idx_mid])

  beta_mat <- rbind(
    beta,
    beta_mid,
    beta_low,
    beta_high
  )

  rownames(beta_mat) <- c(
    "beta",
    paste0("beta_mid_", 1 - 2*tau),
    paste0("beta_low_", tau),
    paste0("beta_high_", tau)
  )

  colnames(beta_mat) <- var_cols

  # Treynor ratio (assets)
  mean_ret_assets <- colMeans(var_mat, na.rm = TRUE)
  treynor_assets  <- (mean_ret_assets - rf) / beta

  treynor_row <- matrix(
    treynor_assets,
    nrow = 1,
    dimnames = list("treynor", var_cols)
  )

  if (!is.null(weights)) {

    beta_p <- as.numeric(beta_mat %*% weights)

    port_ret   <- as.numeric(var_mat %*% weights)
    mean_ret_p <- mean(port_ret, na.rm = TRUE)

    treynor_p <- (mean_ret_p - rf) / sum(beta * weights)

    beta_mat    <- cbind(beta_mat, Portfolio = beta_p)
    treynor_row <- cbind(treynor_row, Portfolio = treynor_p)
  }

  beta_mat <- rbind(beta_mat, treynor_row)

  beta_df <- round(as.data.frame(beta_mat), 4)

  return(beta_df)
}
