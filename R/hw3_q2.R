#' Blood genotype inference
#'
#' \code{blood_type_em} estimates population-level allele frequencies behind observed individuals' blood types using the EM algorithm. Specifically, we assume a Hardy-Weinberg equilibrium to relate the allele combinations (AA, BB, OO, AO, BO, AB) to allele frequencies (\eqn{p_A, p_B, p_O}).
#'
#' @param n a named integer vector containing the observed counts of each blood type. Requires names are ("a", "ab", "b", "o").
#' @param p_initial a named numerical vector containing the initial values to use when optimizing for allele frequencies. Required names are ("a", "b", "o").
#' @param tolerance used for the stopping criterion in optimization.
#'
#' @return a named vector containing the estimated allele frequencies.
#' 
#' @examples
#' blood_type_em(
#'     n = c(a = 6, ab = 4, b = 55, o = 35),
#'     p_initial = c(a = 1 / 3, b = 1 / 3, o = 1 / 3)
#' )
#' 
#' @export
blood_type_em <- function(n, p_initial, tolerance = 1e-3) {
    p <- p_initial

    repeat {
        m_aa <- n["a"] * p["a"]^2 / (p["a"]^2 + 2 * p["a"] * p["o"])
        m_ao <- n["a"] * 2 * p["a"] * p["o"] / (p["a"]^2 + 2 * p["a"] * p["o"])
        m_bb <- n["b"] * p["b"]^2 / (p["b"]^2 + 2 * p["b"] * p["o"])
        m_bo <- n["b"] * 2 * p["b"] * p["o"] / (p["b"]^2 + 2 * p["b"] * p["o"])

        p_old <- p

        p["a"] <- (2 * m_aa + m_ao + n["ab"]) / (2 * sum(n))
        p["b"] <- (2 * m_bb + m_ao + n["ab"]) / (2 * sum(n))
        p["o"] <- (2 * n["o"] + m_ao + m_bo) / (2 * sum(n))

        if (sqrt(sum((p - p_old)^2)) < tolerance)
            break

    }

    p
}