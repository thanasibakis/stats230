#' Hidden Markov models
#'
#' \code{hmm} collects the parameters of a hidden Markov model for use with HMM-related functions in this package, including \code{sample_hmm} and \code{baum_welch}.
#'
#' @param initial_dist a vector of initial probabilities for each hidden state.
#' @param transition_probs a matrix of probabilities for transitioning between hidden states. The (i,j)th element corresponds to P(j|i).
#' @param emission_probs a matrix of probabilities for emitting obserbed values at hidden states. The (i,j)th element corresponds to E(j|i).
#' @param states an optional vector of names for the hidden states. This is used to populate the column and row names of the parameter vectors appropriately.
#' @param observed_values an optional vector of names for the possible observed values. This is used to populate the column and row names of the parameter vectors appropriately.
#'
#' @return an object of class "hmm" for use with functions in this package.
#'
#' @examples
#' model <- hmm(
#'     initial_dist     = c(1/2, 1/2),
#'     transition_probs = rbind(c(0.98, 0.02), c(0.05, 0.95)),
#'     emission_probs   = rbind(rep(1/6, 6), c(1/10, 1/10, 1/2, 1/10, 1/10, 1/10)),
#'     states           = c("A", "B")
#' )
#'
#' @export
hmm <- function(initial_dist, transition_probs, emission_probs,
                states = NULL,
                observed_values = NULL) {
    object <- list(
        # Initial distribution of hidden states
        v = initial_dist,

        # Transition probabilities between hidden states
        p = transition_probs,

        # Emission probabilities (P(observed value | hidden state))
        e = emission_probs
    )

    if (is.null(states))
        states <- seq_along(initial_dist)

    if (is.null(observed_values))
        observed_values <- seq_len(dim(emission_probs)[2])

    # The possible hidden states
    object$states <- states

    # The support of the observed variable
    object$observed_values <- observed_values

    names(object$v)    <- states
    colnames(object$p) <- states
    rownames(object$p) <- states
    rownames(object$e) <- states
    colnames(object$e) <- observed_values

    class(object) <- "hmm"

    object
}

# S3 method.
#' @export
print.hmm <- function(x) {
    cat("States: ", x$states, "\n")
    cat("Observed values: ", x$observed_values, "\n\n")

    cat("Initial distribution:\n")
    print(round(x$v, 4))

    cat("\nTransition probabilities:\n")
    print(round(x$p, 4))

    cat("\nEmission probabilities:\n")
    print(round(x$e, 4))
}

#' Hidden Markov models
#'
#' \code{sample_hmm} simulates realizations of hidden states and observed values from a hidden Markov model.
#'
#' @param n the number of realizations to generate.
#' @param model an object of type "hmm"; see \code{hmm}.
#'
#' @return a tibble with columns for generated hidden states and observed values.
#'
#' @examples
#' model <- hmm(
#'     initial_dist     = c(1/2, 1/2),
#'     transition_probs = rbind(c(0.98, 0.02), c(0.05, 0.95)),
#'     emission_probs   = rbind(rep(1/6, 6), c(1/10, 1/10, 1/2, 1/10, 1/10, 1/10)),
#'     states           = c("A", "B")
#' )
#'
#' sample_hmm(100, model)
#'
#' @export
sample_hmm <- function(n, model) {
    # Accumulate simulated values
    states       <- rep(0, n)
    observations <- rep(0, n)

    # Track the changing distribution of hidden states
    # (initialized using the model initial distribution)
    state_dist <- model$v

    # Simulate (initial sample and transitions)
    for (i in seq_len(n)) {
        states[i] <- sample(
            model$states,
            size = 1,
            prob = state_dist
        )

        observations[i] <- sample(
            model$observed_values,
            size = 1,
            prob = model$e[states[i], ]
        )

        state_dist <- model$p[states[i], ]
    }

    tibble::tibble(
        index = 1:n,
        states = states,
        observations = observations
    )
}

#' Hidden Markov models
#'
#' \code{hmm_backward} runs the HMM backward algorithm for likelihood computation.
#'
#' @param model an object of type "hmm"; see \code{hmm}.
#' @param observations a vector of observed values.
#'
#' @return the matrix \eqn{B} in the context of the algorithm, with number of rows equal to \code{length(observations)} and number of columns equal to the number of hidden states.
#'
#' @examples
#' model <- hmm(
#'     initial_dist     = c(1/2, 1/2),
#'     transition_probs = rbind(c(0.98, 0.02), c(0.05, 0.95)),
#'     emission_probs   = rbind(rep(1/6, 6), c(1/10, 1/10, 1/2, 1/10, 1/10, 1/10)),
#'     states           = c("A", "B")
#' )
#'
#' observations <- sample_hmm(100, model)$observations
#' 
#' hmm_backward(model, observations)
#'
#' @export
hmm_backward <- function(model, observations) {
    p <- model$p
    e <- model$e
    n <- length(observations)

    b <- matrix(0, nrow = n, ncol = length(model$states))
    colnames(b) <- model$states
    b[n, ] <- 1

    for (t in (n - 1):1) {
        for (state in model$states) {
            b[t, state] <- sum(
                p[state, model$states] *
                e[model$states, observations[t + 1]] *
                b[t + 1, model$states]
            )
        }
    }

    b
}

#' Hidden Markov models
#'
#' \code{hmm_forward} runs the HMM forward algorithm for likelihood computation.
#'
#' @param model an object of type "hmm"; see \code{hmm}.
#' @param observations a vector of observed values.
#'
#' @return the matrix \eqn{A} in the context of the algorithm, with number of rows equal to \code{length(observations)} and number of columns equal to the number of hidden states.
#'
#' @examples
#' model <- hmm(
#'     initial_dist     = c(1/2, 1/2),
#'     transition_probs = rbind(c(0.98, 0.02), c(0.05, 0.95)),
#'     emission_probs   = rbind(rep(1/6, 6), c(1/10, 1/10, 1/2, 1/10, 1/10, 1/10)),
#'     states           = c("A", "B")
#' )
#'
#' observations <- sample_hmm(100, model)$observations
#' 
#' hmm_forward(model, observations)
#'
#' @export
hmm_forward <- function(model, observations) {
    v <- model$v
    p <- model$p
    e <- model$e
    n <- length(observations)

    a <- matrix(0, nrow = n, ncol = length(model$states))
    colnames(a) <- model$states
    a[1, ] <- v[model$states] * e[model$states, observations[1]]

    for (t in 1:(n - 1)) {
        for (state in model$states) {
            a[t + 1, state] <- e[state, observations[t + 1]] *
                               sum(a[t, model$states] * p[model$states, state])
        }
    }

    a
}

#' Hidden Markov models
#'
#' \code{baum_welch} runs the Baum-Welch algorithm for HMM parameter inference.
#'
#' @param observations a vector of observed values.
#' @param initial_model an object of type "hmm"; see \code{hmm}. Contains the initial values of the HMM parameters to use for optimization.
#' @param tolerance used for the stopping criterion in optimization.
#' @param max_iterations the maximum number of iterations to allow during optimization. A warning will be displayed if convergence is not achieved by then.
#' @param callback an optional function. If a subset of the parameters' true values are known before optimization time, the callback can be used to specify the true values and prevent optimization of those parameters. The callback must be a function of an "hmm" object and return an "hmm" object with the known parameter values replaced.
#'
#' @return an object of type "hmm" with the locally optimal parameter values.
#'
#' @examples
#' model <- hmm(
#'     initial_dist     = c(1/2, 1/2),
#'     transition_probs = rbind(c(0.98, 0.02), c(0.05, 0.95)),
#'     emission_probs   = rbind(rep(1/6, 6), c(1/10, 1/10, 1/2, 1/10, 1/10, 1/10)),
#'     states           = c("A", "B")
#' )
#'
#' observations <- sample_hmm(100, model)$observations
#' 
#' initial_model <- hmm(
#'      initial_dist     = c(1/4, 3/4),
#'      transition_probs = rbind(c(1/2, 1/2), c(1/2, 1/2)),
#'      emission_probs   = rbind(rep(1/6, 6), rep(1/6, 6)),
#'      states           = c("A", "B")
#'  )
#' 
#' baum_welch(observations, initial_model)
#'
#' # Assume we know the true initial distribution:
#' baum_welch(observations, initial_model,
#'     callback = \(m) { m$v <- model$v; m }
#' )
#'
#' @export
baum_welch <- function(observations, initial_model,
                       tolerance = 1e-3,
                       max_iterations = 1000,
                       callback = NULL) {
    n <- length(observations)
    model <- initial_model

    iteration <- 0

    repeat {
        iteration <- iteration + 1

        # Current parameter values
        v <- model$v
        p <- model$p
        e <- model$e

        # Useful values for parameter updates
        a <- hmm_forward(model, observations)
        b <- hmm_backward(model, observations)

        gamma <- a * b / apply(a * b, 1, sum)

        g <- array(1,
            dim = c(dim(a), dim(a)[2]),
            dimnames = list(
                rownames(a), # t
                colnames(a), # i
                colnames(a)  # j
            )
        )

        for (t in 2:n)
            for (i in model$states)
                for (j in model$states)
                    g[t, i, j] <- b[t, j] *
                                  e[j, observations[t]] *
                                  p[i, j] *
                                  a[t - 1, i]

        g <- sweep(g, 1, apply(g, 1, sum), `/`)

        # Update parameter values
        new_model <- hmm(
            initial_dist     = gamma[1, ],
            transition_probs = sweep(
                apply(g[2:n, , ], 2:3, sum),
                1,
                apply(gamma[1:(n-1), ], 2, sum),
                `/`
            ),
            emission_probs   = e, # Will fill in with correct values
            states           = model$states,
            observed_values  = model$observed_values
        )

        for (state in model$states) {
            for (obs in model$observed_values) {
                new_model$e[state, obs] <-
                    sum(gamma[, state] * (observations == obs)) /
                    sum(gamma[, state])
            }
        }

        # If we're given any true parameter values, swap those in now
        if (!is.null(callback))
            new_model <- callback(new_model)

        model <- new_model

        L2 <- function(a, b) sqrt(sum((a - b)^2))

        if (L2(v, new_model$v) + L2(p, new_model$p) + L2(e, new_model$e) <= tolerance)
            break

        if (iteration == max_iterations) {
            warning("Max iterations reached without convergence")
            break
        }
    }

    model
}