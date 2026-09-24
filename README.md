# Trachoma Teaching Models (WODIN config)

WODIN site configuration for teaching the transmission dynamics of
trachoma, following the structure documented at
https://epimodels.dide.ic.ac.uk/demo/ ("General layout and configuration")
and the example repos https://github.com/mrc-ide/wodin-demo-config and
https://github.com/mrc-ide/wodin-shortcourse-2026.



## Apps

- `apps/sis` - simple 2-compartment (S, I) SIS model, no age structure, no
  partial immunity. Ported from the `teaching_app` Shiny prototype
- `apps/ladder` - adds the infection-history "ladder" (partial immunity via
  rung-dependent recovery/infectiousness) and irreversible disease
  progression (TS, TT), no age structure. Ported from `teaching_app_ladder`.
  **Confirmed WODIN does not support odin's array syntax** (`S[i]`,
  `dim(S) <- n`) - an array-based version compiled/ran fine with plain
  `odin` locally but was rejected by WODIN's own code validator. Rewritten
  as fully unrolled scalar equations, fixed at 5 rungs (`S_1`..`S_5`,
  `I_1`..`I_5`) - changing the rung count means hand-editing/regenerating
  the file rather than adjusting a parameter.

Both apps include MDA as a **periodic pulse approximation** - a narrow,
repeating Gaussian bump in the cure rate built from an explicit `time`
state variable (`deriv(time) <- 1`, `initial(time) <- 0` - odin doesn't
expose the integration time by default, confirmed by the WODIN team),
centred at every multiple of `mda_interval`, with total dose per round
calibrated to clear the same fraction as an instantaneous round with the
given `mda_coverage`/`mda_efficacy` would. This reproduces the
characteristic sawtooth (sharp drop at each round, rebound in between),
since WODIN's "basic" app type has no UI for scheduling true discrete
treatment events.
 

