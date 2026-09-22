# Trachoma Teaching Models (WODIN config)

WODIN site configuration for teaching the transmission dynamics of
trachoma, following the structure documented at
https://epimodels.dide.ic.ac.uk/demo/ ("General layout and configuration")
and the example repos https://github.com/mrc-ide/wodin-demo-config and
https://github.com/mrc-ide/wodin-shortcourse-2026.

This repo contains config only - it is deployed by the WODIN team against
their running `wodin` + `odin.api` + `redis` stack (see the `wodin-demo`
repo for the deployment mechanics); it is not a standalone app you can run
locally without that stack.

## Apps

- `apps/sis` - simple 2-compartment (S, I) SIS model, no age structure, no
  partial immunity, no MDA. Ported from the `teaching_app` Shiny prototype
  in the main modelling repo (`SIR/teaching_app`).

## Status / open questions for the WODIN team

- `wodin.config.json`'s `baseUrl` is a placeholder
  (`https://epimodels.dide.ic.ac.uk/trachoma`) - needs confirming/assigning
  once this is deployed to the sandbox.
- MDA (mass drug administration) is deliberately left out of the `sis` app
  for now - `appType: "basic"` doesn't appear to expose scheduled discrete
  events via the UI the way our Shiny prototype's MDA rounds do. Ask
  whether/how this can be represented in odin/WODIN before adding it.
- A second app porting the "SIS ladder" model (partial immunity,
  irreversible disease progression - see `SIR/teaching_app_ladder`) is
  planned, pending confirmation of whether odin arrays are supported by
  WODIN (if not, that model needs to be code-generated as a fully unrolled
  set of scalar equations).

## Requesting a sandbox deployment

Once this looks reasonable, ask the WODIN team to deploy it to their dev
sandbox (per their standing offer) rather than production, so it can be
checked/debugged before going live.
