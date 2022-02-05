#' Random generation for the multivariate normal distribution
#'
#' Given a dimension \eqn{n} mean vector and a dimension \eqn{n*n} covariance matrix, \code{rmnvorm} generates N realizations from the n-dimensional MVN distribution by transforming \eqn{N*n} iid standard normal samples, using the mean vector and the Cholesky decomposition of the covariance matrix.
#'
#' @param N     number of dimension \eqn{n} samples to generate.
#' @param mu    mean vector, dimension \eqn{n}.
#' @param Sigma covariance matrix, dimension \eqn{n*n}.
#'
#' @return N realizations from \eqn{MVN(μ, Σ)}, as a dimension \eqn{N*n} matrix. Each row is an observation.
#'
#' @examples
#' mu <- c(1, 5, 10, 50)
#' Sigma <- rbind(
#'     c(1.0,   0.5,  0.25, 0.125),
#'     c(0.5,   1.0,  0.5,  0.25 ),
#'     c(0.25,  0.5,  1.0,  0.5  ),
#'     c(0.125, 0.25, 0.5,  1.0  )
#' )
#' 
#' rmvnorm(100, mu, Sigma)
#'
#' @export
rmvnorm <- function(N, mu, Sigma) {
    L <- t(chol(Sigma))

    n <- length(mu)
    Z <- stats::rnorm(N * n) %>%
        matrix(nrow = n, ncol = N)

    Mu <- rep(mu, N) %>%
        matrix(nrow = n)

    X <- Mu + L %*% Z

    # Transpose so rows are observations
    t(X)
}

#' Fitting linear models using a QR decomposition
#'
#' \code{lm_qr} obtains OLS estimates of a linear regression model using the QR decomposition of the design matrix, \eqn{X = QR}.
#' The estimates solve the equation \eqn{R \beta = Q'y}.
#'
#' @param formula an object of class "formula": a symbolic description of the model to be fitted. See documentation for \code{lm} for details.
#' @param data    a data frame containing the response and covariates in the model.
#'
#' @return a named vector of intercept and coefficient estimates.
#'
#' @examples
#' lm_qr(mpg ~ wt + hp, data = mtcars)
#'
#' @export
lm_qr <- function(formula, data) {
    y <- data[[formula[[2]]]]
    X <- stats::model.matrix(formula, data)

    decomp <- qr(X)
    Q <- qr.Q(decomp)
    R <- qr.R(decomp)

    solve(R, t(Q) %*% y)[, 1]
}

#' Fitting linear models using a singular value decomposition
#'
#' \code{lm_qr} obtains OLS estimates of a linear regression model using the SVD of the design matrix, \eqn{X = UDV'}.
#' The estimates solve the equation \eqn{V' \beta = D^{-1} U'y}.
#'
#' @param formula an object of class "formula": a symbolic description of the model to be fitted. See documentation for \code{lm} for details.
#' @param data    a data frame containing the response and covariates in the model.
#'
#' @return a named vector of intercept and coefficient estimates.
#'
#' @examples
#' lm_svd(mpg ~ wt + hp, data = mtcars)
#'
#' @export
lm_svd <- function(formula, data) {
    y <- data[[formula[[2]]]]
    X <- stats::model.matrix(formula, data)

    decomp <- svd(X)
    U <- decomp$u
    D <- diag(decomp$d)
    V <- decomp$v

    # Ensure the righthand side is calculated with
    # matrix*vector multiplications for efficiency
    rhs <- solve(D) %*% (t(U) %*% y)
    solve(t(V), rhs)[, 1] %>%
        magrittr::set_names(colnames(X))
}