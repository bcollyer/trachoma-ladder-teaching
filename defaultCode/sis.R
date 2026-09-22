# Simple population-level SIS model for trachoma (no age structure, no
# partial-immunity ladder, no MDA) - the WODIN port of ../teaching_app.
#
# S + I = N is held constant (no births/deaths). Note: unlike the Shiny
# prototype, this version has no MDA - WODIN's "basic" app type has no UI
# for scheduling discrete treatment events, so that's deferred until we've
# confirmed with the WODIN team how (or whether) to add it.

deriv(S) <- -beta * S * I / N + gamma * I
deriv(I) <-  beta * S * I / N - gamma * I

initial(S) <- N - I0
initial(I) <- I0

N <- user(1000)
I0 <- user(50)
beta <- user(4)
gamma <- user(2)

# Prevalence and R0 as extra traces, so both are visible on the same plot
# without needing a separate readout (R0 renders as a flat reference line).
output(prevalence) <- I / N
output(R0) <- beta / gamma
