#' Multiply two square matrices and a vector
#' 
#' \code{multiply_matrices_and_vector} returns the product of two square matrices and a vector.
#' The product can be computed in two ways, either by first multiplying the matrices, or by first
#' multiplying the second matrix by the vector.
#' 
#' @param A,B The matrices to multiply.
#' @param x   The vector to multiply.
#' @param matrices_first If \code{TRUE}, evaluates \eqn{AB} first, otherwise evaluates \eqn{Bx} first. Useful for benchmarking the performance of the two multiplication algorithms.
#' 
#' @return The product \eqn{ABx}.
#' 
#' @examples
#' A <- matrix(1:9,   nrow = 3)
#' B <- matrix(11:19, nrow = 3)
#' x <- 1:3
#' 
#' multiply_matrices_and_vector(A, B, x)
#'
#' @export
multiply_matrices_and_vector <- function(A, B, x, matrices_first = T) {
    if (matrices_first)
        (A %*% B) %*% x
    else
        A %*% (B %*% x)
}