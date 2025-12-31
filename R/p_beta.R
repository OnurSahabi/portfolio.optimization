#' @export
p_beta <- function(data_var, data_index, var_cols, index_col, weights = NULL, tau = 0.05) {

  index_vec <- as.vector(data_index[, index_col])
  var_mat <- data_var[, var_cols]

  # all sample
  beta <- cov(var_mat, index_vec)[, 1] / var(index_vec)

  # kuyruk eşikleri
  q_low  <- quantile(index_vec, tau)
  q_high <- quantile(index_vec, 1 - tau)

  # low tail
  idx_low <- index_vec <= q_low
  beta_low <- cov(
    var_mat[idx_low, , drop = FALSE],
    index_vec[idx_low]
  )[, 1] / var(index_vec[idx_low])

  # high tail
  idx_high <- index_vec >= q_high
  beta_high <- cov(
    var_mat[idx_high, , drop = FALSE],
    index_vec[idx_high]
  )[, 1] / var(index_vec[idx_high])

  # MIDDLE (inside)
  idx_mid <- index_vec > q_low & index_vec < q_high
  beta_mid <- cov(
    var_mat[idx_mid, , drop = FALSE],
    index_vec[idx_mid]
  )[, 1] / var(index_vec[idx_mid])

  # table format
  beta_mat <- rbind(
    beta,
    beta_mid,
    beta_low,
    beta_high
  )

  rownames(beta_mat) <- c(
    "beta_all",
    paste0("beta_mid_", 1 - 2*tau),
    paste0("beta_low_", tau),
    paste0("beta_high_", tau)
  )

  colnames(beta_mat) <- var_cols

  if (!is.null(weights)) {
    beta_mat <- cbind(
      beta_mat,
      Portfolio = as.numeric(beta_mat %*% weights)
    )
  }

  beta_mat <- round(beta_mat, 4)
  beta_df  <- as.data.frame(beta_mat)

  return(beta_df)
}
