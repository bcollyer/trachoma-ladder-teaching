# SIS "ladder" model for trachoma (infection-history rungs, no age
# structure). Recovery gets faster and
# infectiousness gets weaker at higher rungs (nu0/nu1/nuExponent,
# lExponent). Recovery from the top rung (5) loops back to S_5 (capped).

# Model parameters

beta <- user(4.5)                  # transmission rate
nu0 <- user(0.885594)            # minimum recovery rate
nu1 <- user(4.208717)            # maximum recovery rate
nuExponent <- user(0.751195)     # exponential increase in recovery rate
lExponent <- user(0.05277)       # exponential decrease in infectiousness
r_TS <- user(0.10)               # rate of scarring
r_TT <- user(0.02)               # rate of trichiasis
alpha <- user(0.04)              # rate of surgery
life_expectancy <- user(55)      # life expectancy
mu <- 1 / life_expectancy        # death/birth rate
N <- user(1000)                  # Total population size
I0_prop <- user(0.05)            # Initial infected proportion
mda_coverage <- user(0.75)
mda_efficacy <- user(0.95)
mda_interval <- user(1)          # Time between MDA rounds (years)
mda_pulse_width <- 0.05    
mda_start <- user(50)

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

# Level-dependent recovery rate
nu_1 <- nu0
nu_2 <- (nu0 - nu1) * exp(-nuExponent * 1) + nu1
nu_3 <- (nu0 - nu1) * exp(-nuExponent * 2) + nu1
nu_4 <- (nu0 - nu1) * exp(-nuExponent * 3) + nu1
nu_5 <- (nu0 - nu1) * exp(-nuExponent * 4) + nu1

# Level-dependent infectiousness rate
lload_1 <- 1
lload_2 <- exp(-lExponent * 1)
lload_3 <- exp(-lExponent * 2)
lload_4 <- exp(-lExponent * 3)
lload_5 <- exp(-lExponent * 4)

# Force of infection
lambda <- beta * (lload_1 * I_1 + lload_2 * I_2 + lload_3 * I_3 +
                     lload_4 * I_4 + lload_5 * I_5) / N

# Incidence
incidence_1 <- S_1 * lambda
incidence_2 <- S_2 * lambda
incidence_3 <- S_3 * lambda
incidence_4 <- S_4 * lambda
incidence_5 <- S_5 * lambda

# Recoveries
recovery_1 <- I_1 * nu_1
recovery_2 <- I_2 * nu_2
recovery_3 <- I_3 * nu_3
recovery_4 <- I_4 * nu_4
recovery_5 <- I_5 * nu_5

# Model derivatives
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


# Initital state
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

# Outputs
output(prevalence_infection) <- (I_1 + I_2 + I_3 + I_4 + I_5) / N
output(prevalence_TT) <- TT / N
output(prevalence_TS) <- TS / N
output(prevalence_threshold_5_percent) <- 0.05