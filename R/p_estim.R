#' Constrained Portfolio Estimation via Convex Optimization
#'
#' Solves a long-only portfolio optimization problem using convex programming
#' via the \code{CVXR} package.
#'
#' The function supports two alternative problems:
#' \itemize{
#'   \item Maximize expected return subject to a risk constraint.
#'   \item Minimize portfolio risk subject to a target return constraint.
#' }
#'
#' @param data A numeric matrix or data.frame of asset returns
#'   (rows = time, columns = assets).
#' @param risk Optional numeric scalar. Maximum allowable portfolio risk.
#'   If supplied, the function maximizes expected return subject to this risk bound.
#' @param return Optional numeric scalar. Target portfolio return.
#'   If supplied, the function minimizes portfolio risk subject to this return level.
#' @param rf Numeric. Risk-free rate used in Sharpe ratio calculation (default = 0).
#' @param digits Integer. Number of decimal places used to round portfolio weights.
#'
#' @return A list containing:
#' \describe{
#'   \item{weights}{Optimal portfolio weights.}
#'   \item{metrics}{Named vector with Return, Risk, and Sharpe ratio.}
#' }
#'
#' @details
#' Exactly one of \code{risk} or \code{return} must be provided.
#' The optimization problem is solved using quadratic cone programming
#' via \code{CVXR}.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' R <- matrix(rnorm(300), ncol = 4)
#' colnames(R) <- paste0("Asset", 1:4)
#'
#' # Max return given risk constraint
#' p_estim(R, risk = 0.02)
#'
#' # Min risk given target return
#' p_estim(R, return = 0.001)
#' }
#'
#' @importFrom stats cov colMeans
#' @export
p_estim <- function(data, risk = NULL, return = NULL, rf = 0, digits = 2) {

  if (!requireNamespace("CVXR", quietly = TRUE))
    stop("CVXR kurulu değil")

  if (xor(is.null(risk), is.null(return)) == FALSE)
    stop("Sadece risk veya sadece return vermelisin")

  data <- as.matrix(data)
  mu <- colMeans(data)
  S  <- cov(data)
  L  <- chol(S)
  n  <- length(mu)

  w <- CVXR::Variable(n)

  obj <- if (is.null(return)) {
    CVXR::Maximize(t(mu - rf) %*% w)
  } else {
    CVXR::Minimize(CVXR::norm2(L %*% w))
  }

  cons <- list(
    w >= 0,
    sum(w) == 1
  )

  cons <- if (is.null(return)) {
    c(cons, list(CVXR::norm2(L %*% w) <= risk))
  } else {
    c(cons, list(t(mu) %*% w == return))
  }

  sol <- CVXR::solve(CVXR::Problem(obj, cons))

  w_hat <- as.numeric(sol$getValue(w))
  names(w_hat) <- colnames(data)

  Return  <- sum(mu * w_hat)
  Risk <- sqrt(drop(t(w_hat) %*% S %*% w_hat))
  Sharpe = (Return - rf) / Risk

  list(
    weights = round(w_hat, digits),
    metrics = round(c(Return = Return, Risk = Risk, Sharpe = Sharpe), 6)
  )
}
