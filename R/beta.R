#' @export
beta <- function(data1, data2, stock_cols, index_col) {

  index_vec <- as.vector(data2[, index_col])

  betas <- cov(data1[, stock_cols], index_vec)[, 1] /
    var(index_vec)

  return(betas)
}
