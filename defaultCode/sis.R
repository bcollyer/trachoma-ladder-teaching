# Simple population-level SIS model for trachoma (no age structure, no
# partial-immunity ladder) - the WODIN port of ../teaching_app.
#
# S + I = N is held constant (no births/deaths). MDA is a CONTINUOUS
# approximation, not scheduled discrete events (WODIN's "basic" app type
# has no UI for those) - see mda_rate below.

# Annual MDA reproduced as a periodic pulse rather than a flat rate:
# odin exposes the integration time as `t`, so a narrow, repeating
# Gaussian centred at each multiple of mda_interval gives the sawtooth
# drop-then-rebound shape of real MDA rounds, with total "dose" per cycle
# calibrated to clear the same fraction p = coverage*efficacy as an
# instantaneous round would. Capped at p=0.995 so the dose stays finite.
mda_p <- min(mda_coverage * mda_efficacy, 0.995)
mda_target <- -log(1 - mda_p)              # total hazard to deliver per cycle
mda_sigma <- max(mda_pulse_width, 0.001) * mda_interval  # pulse width (time units)
mda_cycle_frac <- t / mda_interval - floor(t / mda_interval)
mda_dist <- min(mda_cycle_frac, 1 - mda_cycle_frac) * mda_interval  # time to nearest pulse centre, wrapped
mda_shape <- exp(-(mda_dist * mda_dist) / (2 * mda_sigma * mda_sigma))
mda_rate <- mda_target / (mda_sigma * sqrt(2 * 3.14159265358979)) * mda_shape

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
mda_pulse_width <- user(0.05)

# Prevalence and R0 as extra traces, so both are visible on the same plot
# without needing a separate readout (R0 renders as a flat reference line).
output(prevalence) <- I / N
output(R0) <- beta / gamma

