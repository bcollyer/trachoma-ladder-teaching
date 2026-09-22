# Simple population-level SIS model for trachoma (no age structure, no
# partial-immunity ladder) - the WODIN port of ../teaching_app.
#
# S + I = N is held constant (no births/deaths). MDA is a CONTINUOUS
# approximation, not scheduled discrete events (WODIN's "basic" app type
# has no UI for those) - see mda_rate below.

# Annual MDA as an equivalent continuous cure rate: a pulse clearing
# fraction p = coverage*efficacy once every mda_interval years clears the
# same total fraction per year as a continuous rate
# tau = -log(1-p)/mda_interval. Capped at p=0.995 so tau stays finite if
# coverage*efficacy is set to (or very near) 1.
mda_p <- min(mda_coverage * mda_efficacy, 0.995)
mda_rate <- -log(1 - mda_p) / mda_interval

deriv(S) <- -beta * S * I / N + gamma * I + mda_rate * I
deriv(I) <-  beta * S * I / N - gamma * I - mda_rate * I

initial(S) <- N - I0
initial(I) <- I0

N <- user(1000)
I0 <- user(50)
beta <- user(4)
gamma <- user(2)
mda_coverage <- user(0)
mda_efficacy <- user(1)
mda_interval <- user(1)

# Prevalence and R0 as extra traces, so both are visible on the same plot
# without needing a separate readout (R0 renders as a flat reference line).
output(prevalence) <- I / N
output(R0) <- beta / gamma

