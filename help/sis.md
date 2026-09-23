## What is this model?

This app is a **teaching tool** for exploring the transmission dynamics of
ocular *Chlamydia trachomatis* (trachoma) using a simple **SIS
(Susceptible-Infected-Susceptible)** compartmental model.

Trachoma is modelled as SIS rather than SIR because infection does **not**
confer lasting protective immunity. Longitudinal field studies (e.g. Bailey
et al. 1999; Grassly et al. 2008) show that individuals are repeatedly
re-infected throughout life. An SIR model would predict a single epidemic
wave followed by elimination, which is not what is observed in
trachoma-endemic communities.

This is intentionally a **simplified teaching model** - no age structure,
no partial immunity - but it does include **Mass Drug Administration
(MDA)**, approximated as a continuous cure rate (see below) since WODIN's
"basic" app type has no UI for scheduling discrete treatment events.

### Compartments

- **S** - susceptible individuals (currently uninfected, can become infected)
- **I** - infected individuals (currently infected, infectious to others)

The population size $N = S + I$ is held constant (no births/deaths in this
simplified teaching version).

### Diagram

![SIS compartment diagram](sis_diagram.svg)

**S** flows to **I** via the force of infection ($\beta I/N$); **I** flows
back to **S** through natural clearance ($\gamma$) or MDA (see below).

### Equations

\begin{align}
\frac{dS}{dt} &= -\beta \frac{I}{N} S + \gamma I + \tau(t) I \\[6pt]
\frac{dI}{dt} &= \beta \frac{I}{N} S - \gamma I - \tau(t) I
\end{align}

where $\tau(t)$ is the (time-varying) MDA cure rate defined below.

| Symbol | Meaning | Units |
|---|---|---|
| $\beta$ | effective transmission rate | per year |
| $\gamma$ | recovery / spontaneous clearance rate | per year |
| $N$ | population size | people |
| mda_coverage | fraction of infections cleared per MDA round | dimensionless (0-1) |
| mda_efficacy | drug efficacy per round | dimensionless (0-1) |
| mda_interval | years between MDA rounds | years |
| mda_pulse_width | sharpness of each round, as a fraction of `mda_interval` | dimensionless (small, e.g. 0.03) |
| $R_0 = \beta/\gamma$ | basic reproduction number | dimensionless |

### Mass Drug Administration (periodic pulse approximation)

WODIN's "basic" app type has no UI for scheduling discrete treatment
events at exact times, but `odin` does expose the current simulation time
as `t`, so MDA rounds are reproduced as a **narrow, repeating pulse built
directly into the equations** rather than a flat background rate. Each
round is a Gaussian-shaped bump in the cure rate, centred at every
multiple of `mda_interval`, with total "dose" over the round calibrated
so the same fraction of infections is cleared as a single instantaneous
round with coverage $c$ and efficacy $e$ would clear:

$$
p = \min(c \times e,\ 0.995), \qquad \sigma = \text{mda\_pulse\_width} \times \text{mda\_interval}
$$

$$
\tau(t) = \frac{-\ln(1-p)}{\sigma\sqrt{2\pi}} \, \exp\!\left(-\frac{d(t)^2}{2\sigma^2}\right)
$$

where $d(t)$ is the (wrapped) time to the nearest round. This produces
the characteristic **sawtooth**: a sharp drop in prevalence at each round,
followed by a rebound towards the underlying endemic level as new
infections accumulate before the next round. Smaller `mda_pulse_width`
gives a sharper, more "instantaneous-looking" drop; too small a value can
make the pulse hard for the solver to resolve smoothly, so this is capped
at a reasonable minimum internally.

### Equilibrium prevalence

This SIS system has a stable **endemic equilibrium** whenever $R_0 > 1$, at
prevalence:

$$
I^*/N = 1 - \frac{1}{R_0} = 1 - \frac{\gamma}{\beta}
$$

If $R_0 \le 1$ the infection cannot sustain itself and dies out ($I \to 0$).

### How to use the simulator

Edit the parameters (or the code directly) and watch:

- how **prevalence** ($I/N$, plotted alongside S and I) rises to (or falls
  from) its equilibrium value;
- how increasing $\beta$ (more transmission) or decreasing $\gamma$ (slower
  clearance) raises $R_0$ (also plotted, as a flat reference line) and the
  equilibrium prevalence;
- how increasing `mda_coverage`/`mda_efficacy` or shortening `mda_interval`
  pushes the equilibrium prevalence down (set `mda_coverage` to 0 to turn
  MDA off).

See `trachoma_first_principles.R` and `trachoma_first_principles_300_history.R`
in the main modelling repository for the fuller models used for actual
fitting and projection work, and the `teaching_app_ladder` Shiny prototype
for a version with partial immunity and irreversible disease progression.
