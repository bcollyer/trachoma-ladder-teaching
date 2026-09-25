# Simple population-level SIS model for trachoma (no age structure, no
# partial-immunity ladder) 
#
# S + I = N is held constant (no births/deaths). MDA is a CONTINUOUS
# approximation, not scheduled discrete events


# Model parameters
N <- user(1000)          # Total population size
I0 <- user(50)           # Initial infected individuals
beta <- user(4)          # transmission rate
gamma <- user(2)         # recovery rate
mda_coverage <- user(0)
mda_efficacy <- user(1)
mda_interval <- user(1)
mda_pulse_width <- 0.05
mda_start <- user(15)

# This block configures the MDA in a continuous time way. Normally it would be discrete.
deriv(time) <- 1           
initial(time) <- 0
mda_p <- min(mda_coverage * mda_efficacy, 0.995)
mda_target <- -log(1 - mda_p)              # total hazard to deliver per cycle
mda_sigma <- max(mda_pulse_width, 0.001) * mda_interval  # pulse width (time units)
mda_time_since_start <- time - mda_start
mda_cycle_frac <- mda_time_since_start / mda_interval - floor(mda_time_since_start / mda_interval)
mda_dist <- min(mda_cycle_frac, 1 - mda_cycle_frac) * mda_interval  # time to nearest pulse centre, wrapped
mda_shape <- exp(-(mda_dist * mda_dist) / (2 * mda_sigma * mda_sigma))#mda_rate <- mda_target / (mda_sigma * sqrt(2 * 3.14159265358979)) * mda_shape
mda_rate <- if (time < mda_start) 0 else
mda_target / (mda_sigma * sqrt(2 * 3.14159265358979)) * mda_shape


# SIS model derivatives
deriv(S) <- -beta * S * I / N + gamma * I +
         mda_rate * I                       
deriv(I) <-  beta * S * I / N - gamma * I -
         mda_rate * I

# Initials state
initial(S) <- N - I0
initial(I) <- I0



# Prevalence and R0 as extra traces, so both are visible on the same plot
# without needing a separate readout (R0 renders as a flat reference line).
output(prevalence) <- I / N
output(R0) <- beta / gamma
output(prevalence_threshold_5_percent) <- 0.05
