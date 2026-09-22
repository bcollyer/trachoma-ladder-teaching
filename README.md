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
  `n_rungs` sizes the S/I arrays via odin's user-sized-array feature -
  compiles and runs correctly with plain `odin` locally, but whether
  WODIN's UI handles an array-sizing parameter is unconfirmed.
 

