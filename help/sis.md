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
no partial immunity, and (unlike the fuller research models) no mass drug
administration in this version.

### Compartments

- **S** - susceptible individuals (currently uninfected, can become infected)
- **I** - infected individuals (currently infected, infectious to others)

The population size $N = S + I$ is held constant (no births/deaths in this
simplified teaching version).

### Diagram

![SIS compartment diagram](sis_diagram.svg)

**S** flows to **I** via the force of infection ($\beta I/N$); **I** flows
back to **S** through natural clearance ($\gamma$).

### Equations

\begin{align}
\frac{dS}{dt} &= -\beta \frac{I}{N} S + \gamma I \\[6pt]
\frac{dI}{dt} &= \beta \frac{I}{N} S - \gamma I
\end{align}

| Symbol | Meaning | Units |
|---|---|---|
| $\beta$ | effective transmission rate | per year |
| $\gamma$ | recovery / spontaneous clearance rate | per year |
| $N$ | population size | people |
| $R_0 = \beta/\gamma$ | basic reproduction number | dimensionless |

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
  equilibrium prevalence.

See `trachoma_first_principles.R` and `trachoma_first_principles_300_history.R`
in the main modelling repository for the fuller models used for actual
fitting and projection work, and the `teaching_app_ladder` Shiny prototype
for a version with partial immunity and irreversible disease progression.
