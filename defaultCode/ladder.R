# SIS "ladder" model for trachoma (infection-history rungs, no age
# structure) - the odin/WODIN port of ../../teaching_app_ladder.
#
# Fully UNROLLED to plain scalar equations (S_1..S_5, I_1..I_5, TS, TT) -
# NOT the odin array (dim()/S[i]) version originally written here. That
# version compiled and ran fine with the plain `odin` R package, but WODIN
# itself rejected it with:
#   "Code error: ... deriv() and initial() must contain same set of
#   equations: in deriv() but not initial(): S, I, TS, TT"
# i.e. WODIN's own code validator does not understand odin's array
# bracket syntax at all (confirms the "arrays... not supported" line in
# the WODIN docs) - so the number of rungs is fixed at 5 here rather than
# a user()-configurable dimension. To change the rung count, this file
# needs to be regenerated/hand-edited with a different number of S_j/I_j
# pairs, following the same pattern.
#
# Same mechanism as the full 300-history model: recovery gets faster and
# infectiousness gets weaker at higher rungs (nu0/nu1/nuExponent,
# lExponent), rather than a separate susceptibility reduction. Recovery
# from the top rung (5) loops back to S_5 (capped).
#
# MDA is a CONTINUOUS approximation, not scheduled discrete events (WODIN's
# "basic" app type has no UI for those): an annual pulse clearing fraction
# p = mda_coverage*mda_efficacy once every mda_interval years clears the
# same total fraction per year as a continuous rate
# tau = -log(1-p)/mda_interval, applied I_j -> S_j at the SAME rung (MDA
# cures don't build partial immunity, unlike natural recovery).

mda_p <- min(mda_coverage * mda_efficacy, 0.995)
mda_rate <- -log(1 - mda_p) / mda_interval

# Rung-dependent recovery rate and relative infectiousness - the rung
# index (j-1) is baked in per-equation as a literal number.
nu_1 <- nu0
nu_2 <- (nu0 - nu1) * exp(-nuExponent * 1) + nu1
nu_3 <- (nu0 - nu1) * exp(-nuExponent * 2) + nu1
nu_4 <- (nu0 - nu1) * exp(-nuExponent * 3) + nu1
nu_5 <- (nu0 - nu1) * exp(-nuExponent * 4) + nu1

lload_1 <- 1
lload_2 <- exp(-lExponent * 1)
lload_3 <- exp(-lExponent * 2)
lload_4 <- exp(-lExponent * 3)
lload_5 <- exp(-lExponent * 4)

lambda <- beta * (lload_1 * I_1 + lload_2 * I_2 + lload_3 * I_3 +
                     lload_4 * I_4 + lload_5 * I_5) / N

incidence_1 <- S_1 * lambda
incidence_2 <- S_2 * lambda
incidence_3 <- S_3 * lambda
incidence_4 <- S_4 * lambda
incidence_5 <- S_5 * lambda

recovery_1 <- I_1 * nu_1
recovery_2 <- I_2 * nu_2
recovery_3 <- I_3 * nu_3
recovery_4 <- I_4 * nu_4
recovery_5 <- I_5 * nu_5

deriv(S_1) <- mu * N - incidence_1 - mu * S_1 + mda_rate * I_1
deriv(I_1) <- incidence_1 - recovery_1 - mu * I_1 - mda_rate * I_1

deriv(S_2) <- recovery_1 - incidence_2 - mu * S_2 + mda_rate * I_2
deriv(I_2) <- incidence_2 - recovery_2 - mu * I_2 - mda_rate * I_2

deriv(S_3) <- recovery_2 - incidence_3 - mu * S_3 + mda_rate * I_3
deriv(I_3) <- incidence_3 - recovery_3 - mu * I_3 - mda_rate * I_3

deriv(S_4) <- recovery_3 - incidence_4 - mu * S_4 + mda_rate * I_4
deriv(I_4) <- incidence_4 - recovery_4 - mu * I_4 - mda_rate * I_4

# Top rung is capped: recovery from I_5 loops back to S_5, not a rung 6.
deriv(S_5) <- recovery_4 + recovery_5 - incidence_5 - mu * S_5 + mda_rate * I_5
deriv(I_5) <- incidence_5 - recovery_5 - mu * I_5 - mda_rate * I_5

# Scarring (TS) and trichiasis (TT), driven by repeat infection (rungs >= 2).
repeat_I <- I_2 + I_3 + I_4 + I_5
new_TS <- r_TS * (repeat_I / N) * (N - TS)
new_TT <- r_TT * (TS / N) * (N - TT)
deriv(TS) <- new_TS - mu * TS
deriv(TT) <- new_TT - (alpha + mu) * TT

initial(S_1) <- N * (1 - I0_prop)
initial(I_1) <- N * I0_prop
initial(S_2) <- 0
initial(I_2) <- 0
initial(S_3) <- 0
initial(I_3) <- 0
initial(S_4) <- 0
initial(I_4) <- 0
initial(S_5) <- 0
initial(I_5) <- 0
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
mda_coverage <- user(0)
mda_efficacy <- user(1)
mda_interval <- user(1)

output(prevalence_infection) <- (I_1 + I_2 + I_3 + I_4 + I_5) / N
output(prevalence_TT) <- TT / N
