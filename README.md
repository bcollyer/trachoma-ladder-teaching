# Trachoma Teaching Models (WODIN config)

WODIN site configuration for teaching the transmission dynamics of
trachoma, following the structure documented at
https://epimodels.dide.ic.ac.uk/demo/ ("General layout and configuration")
and the example repos https://github.com/mrc-ide/wodin-demo-config and
https://github.com/mrc-ide/wodin-shortcourse-2026.



## Apps

- `apps/sis` - simple 2-compartment (S, I) SIS model, no age structure, no
  partial immunity, no MDA. Ported from the `teaching_app` Shiny prototype
- `apps/ladder` - adds the infection-history "ladder" (partial immunity via
  rung-dependent recovery/infectiousness) and irreversible disease
  progression (TS, TT), no age structure. Ported from `teaching_app_ladder`.
  **Confirmed WODIN does not support odin's array syntax** (`S[i]`,
  `dim(S) <- n`) - an array-based version compiled/ran fine with plain
  `odin` locally but was rejected by WODIN's own code validator. Rewritten
  as fully unrolled scalar equations, fixed at 5 rungs (`S_1`..`S_5`,
  `I_1`..`I_5`) - changing the rung count means hand-editing/regenerating
  the file rather than adjusting a parameter.
 

