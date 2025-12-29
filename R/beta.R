#' @export
beta <- function(stocks_ret, index_ret) {

  betas <- apply(stocks_ret, 2, function(x) {
    cov(x, index_ret) / var(index_ret)
  })

  return(betas)
}
