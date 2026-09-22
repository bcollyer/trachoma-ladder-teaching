# SIS "ladder" model for trachoma (infection-history rungs, no age
# structure) - the odin/WODIN port of ../../teaching_app_ladder.
#
# Same mechanism as the full 300-history model: recovery gets faster and
# infectiousness gets weaker at higher rungs (nu0/nu1/nuExponent,
# lExponent), rather than a separate susceptibility reduction. Recovery
# from the top rung n_rungs loops back to S[n_rungs] (capped).
#
# NOTE: n_rungs is a user()-sized array dimension, fixed when the model is
# initialised (mod$new(n_rungs = ...)) - it cannot be changed by a plain
# parameter slider without re-creating the model instance. Requires
# n_rungs >= 2 (the array-range equations below assume at least 2 rungs).
#
# Like sis.R, this omits MDA - WODIN's "basic" app type has no UI for
# scheduling discrete treatment events.

n_rungs <- user(5, integer = TRUE)
dim(S) <- n_rungs
dim(I) <- n_rungs
dim(nu) <- n_rungs
dim(lload) <- n_rungs
dim(weighted_I) <- n_rungs
dim(incidence) <- n_rungs
dim(recovery) <- n_rungs

# Rung-dependent recovery rate and relative infectiousness.
nu[] <- (nu0 - nu1) * exp(-nuExponent * (i - 1)) + nu1
lload[] <- exp(-lExponent * (i - 1))
weighted_I[] <- lload[i] * I[i]

lambda <- beta * sum(weighted_I) / N
incidence[] <- S[i] * lambda
recovery[] <- I[i] * nu[i]

# S[1]: naive susceptibles, replenished by births.
deriv(S[1]) <- mu * N - incidence[1] - mu * S[1]
# S[2..n_rungs-1]: gain recovery from the rung below.
deriv(S[2:(n_rungs - 1)]) <- recovery[i - 1] - incidence[i] - mu * S[i]
# S[n_rungs]: top rung is capped - recovery from I[n_rungs] loops back here too.
deriv(S[n_rungs]) <- recovery[n_rungs - 1] + recovery[n_rungs] - incidence[n_rungs] - mu * S[n_rungs]

deriv(I[]) <- incidence[i] - recovery[i] - mu * I[i]

# Scarring (TS) and trichiasis (TT), driven by repeat infection (rungs >= 2).
repeat_I <- sum(I[2:n_rungs])
new_TS <- r_TS * (repeat_I / N) * (N - TS)
new_TT <- r_TT * (TS / N) * (N - TT)
deriv(TS) <- new_TS - mu * TS
deriv(TT) <- new_TT - (alpha + mu) * TT

initial(S[1]) <- N * (1 - I0_prop)
initial(S[2:n_rungs]) <- 0
initial(I[1]) <- N * I0_prop
initial(I[2:n_rungs]) <- 0
initial(TS) <- 0
initial(TT) <- 0

beta <- user(3)
nu0 <- user(0.885594)
nu1 <- user(4.208717)
nuExponent <- user(0.751195)
lExponent <- user(0.05277)
r_TS <- user(0.10)
r_TT <- user(0.02)
alpha <- user(0.04)
life_expectancy <- user(55)
mu <- 1 / life_expectancy
N <- user(1000)
I0_prop <- user(0.05)

output(prevalence_infection) <- sum(I) / N
output(prevalence_TT) <- TT / N
