#' Fitting logistic regression models
#'
#' \code{logistic_regression} obtains maximum likelihood estimates of coefficients using either iteratively-reweighted least squares or gradient ascent.
#'
#' @param formula an object of class "formula": a symbolic description of the model to be fitted. See documentation for \code{lm} for details.
#' @param data data frame containing the response and covariates in the model.
#' @param method one of {"IRLS", "gradient"}, to obtain the MLEs using IRLS or gradient ascent.
#' @param tolerance used for the stopping criterion in optimization.
#' @param step_size used for the gradient ascent updates. Unused if \code{method = "IRLS"}.
#' @param max_iterations the maximum number of iterations to allow during optimization. A warning will be displayed if convergence is not achieved by then.
#' @param rescale logical; if \code{TRUE}, numeric columns in the data will be standardized before optimization and unstandardized after. This is for numerical stability, particularly when \code{method = "gradient"}. Coefficient and confidence interval outputs will still be in terms of the original scale.
#'
#' @return an object of class "logistic_regression". Supported S3 methods are \code{coef}, \code{confint}, and \code{plot} (for plotting log likelihoods over optimization time).
#'
#' @examples
#' logistic_regression(I(mpg > 20) ~ disp, data = mtcars)
#'
#' @export
logistic_regression <- function(formula, data,
                                method = "IRLS",
                                tolerance = 1e-4,
                                step_size = 1e-3,
                                max_iterations = 5000,
                                rescale = F) {
    # Scale data for better convergence
    if (rescale)
        data <- dplyr::mutate(data, dplyr::across(where(is.numeric), scale))

    y <- model.frame(formula, data) %>% model.response()
    x <- model.matrix(formula, data)

    # log likelihood function and its gradient/Hessian wrt. beta
    l   <- function(p) dbinom(y, size = 1, prob = p) %>% log() %>% sum()
    dl  <- function(p) t(x) %*% (y - p)
    d2l <- function(p) -t(x) %*% diag(p * (1 - p)) %*% x

    # Optimization step function (for updating beta)
    delta <- switch(method,
        IRLS     = function(p) -solve(d2l(p),  dl(p)),
        gradient = function(p) step_size * dl(p),
        stop(paste(method, "is not a supported optimization method"))
    )

    # Initialize beta by setting logit(p) = y
    beta <- 0
    p    <- 1 / (1 + exp(-y))

    # Track the log likelihoods over iterations
    log_likelihood_trace <- c()
    iteration <- 0

    # Optimize
    repeat {
        iteration <- iteration + 1

        del  <- delta(p)
        beta <- beta + del
        p    <- 1 / (1 + exp(-x %*% beta))[, 1]

        log_likelihood_trace <- c(log_likelihood_trace, l(p))

        if (sqrt(sum(del^2)) < tolerance)
            break

        if (iteration == max_iterations) {
            warning("Max iterations reached without convergence")
            break
        }
    }

    # Unscale the data to obtain the proper coefficients
    if (rescale) {
        scales <- data %>%
            dplyr::select_if(is.numeric) %>%
            lapply(\(x) attr(x, "scaled:scale")) %>%
            c(recursive = T)

        centers <- data %>%
            dplyr::select_if(is.numeric) %>%
            lapply(\(x) attr(x, "scaled:center")) %>%
            c(recursive = T)

        beta[names(scales), 1] <- beta[names(scales), 1] / scales
        beta["(Intercept)", 1] <- beta["(Intercept)", 1] -
                                  sum(beta[names(centers), 1] * centers)

        # Update x and p so d2l can be computed properly
        x[, names(scales)]  <- sweep(x[, names(scales)],  MARGIN = 2, scales,  `*`)
        x[, names(centers)] <- sweep(x[, names(centers)], MARGIN = 2, centers, `+`)
        p <- 1 / (1 + exp(-x %*% beta))[, 1]
    }

    results <- list(
        coefficients         = beta[, 1],
        log_likelihood_trace = log_likelihood_trace,
        vcov_matrix          = solve(-d2l(p))
    )

    class(results) <- "logistic_regression"

    results
}

# S3 method.
#' @export
confint.logistic_regression <- function(object) {
    se <- sqrt(diag(object$vcov_matrix))
    ci <- coef(object) + 1.96 * cbind(-se, se)
    colnames(ci) <- c("2.5 %", "97.5 %")

    ci
}

# S3 method.
#' @export
print.logistic_regression <- function(x) {
    cat("Coefficients and confidence intervals:\n")
    coefs <- as.matrix(coef(x))
    colnames(coefs) <- "Estimate"
    table <- cbind(coefs, confint(x)) %>% round(4)
    print(table)
}

# S3 method.
#' @export
plot.logistic_regression <- function(x) {
    ggplot() +
        geom_point(aes(
            x = seq_along(x$log_likelihood_trace),
            y = x$log_likelihood_trace
        )) +
        labs(
            x = "Iteration",
            y = "log likelihood",
            title = "Trace of log likelihood during optimization"
        )
}