# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`uds8bus` is one benchmark inside the larger `benchmarks_public` collection (see
`../../README.md`). It models an 8-bus low-voltage AC/DC distribution feeder
(4 AC buses `A0`–`A4` and 3 DC buses `D2`–`D4`, plus a breaker) with a
synchronous source, two power-electronics VSCs, a grid-forming inverter (`vsg`)
and a battery (`bess_dcdc_gf`).

Unlike the other benchmarks in this repo (`ieee39`, `nts`, etc.) that follow
the `build()` / `ini()` / `run()` module pattern on `pydae-bps`, this case is
driven by **`pydae.uds`** (`UdsBuilder`) on the CasADi backend. The two stacks
live in different packages and are not interchangeable.

## Running

```bash
python run.py        # build → ini → 1 s simulation → annotated SVG
python opt.py        # parameter sweep over p_vsc_a_ref_A4
```

`run.py` and `opt.py` both **prepend a hard-coded path** to put the local
`pydae-uds` source first on `sys.path`:

```python
sys.path.insert(0, "/Users/jmmauricio/workspace/pydae/packages/pydae-uds/src")
```

This is intentional. The active conda env (`pydae_dev`) ships a different
`pydae.uds` that resolves to another repo, and the insert overrides it. If
you edit either script, preserve this line — running without it silently
imports the wrong package.

## Architecture

The pipeline in `run.py` is the canonical entry point:

1. `UdsBuilder("uds8bus.hjson", use_casadi=True).construct("uds8bus")`
   parses the HJSON network description into `grid.sys_dict`.
2. `CasadiModel(CasadiBuilder(grid.sys_dict).build())` compiles the symbolic
   DAE on the CasADi backend (no C compiler / no CFFI build step).
3. `model.ini(params, xy_0="xy_0.json")` solves the load flow seeded from
   `xy_0.json`. The seed is the previously-converged operating point; if you
   change topology or operating conditions enough that Newton-Raphson stops
   converging, regenerate it.
4. `model.run(t, {})` then `model.post()` integrates the DAE forward.
5. `model2svg(model, "uds8bus.hjson", "uds8bus.svg").set_tooltips(...)`
   annotates the topology SVG with bus voltages and line loadings at the
   current time, writing to a `*_tooltips.svg` so the source SVG template
   stays untouched.

`opt.py` is the same first three steps in a loop that sweeps
`p_vsc_a_ref_A4` from -200 kW to +200 kW in 5 kW steps and reports the slack
power `p_A0` at each point. Reuses the same `model` object across iterations.

## The HJSON network file

`uds8bus.hjson` is the single source of truth for topology, ratings, line
codes, and controller wiring. The shape of `xy_0.json` (variable names like
`V_<bus>_<node>_{r,i}`, `v_t{abc}_{r,i}_<bus>`, controller states) is derived
from this file at build time — keep them in sync.

Notable conventions in this file:

- `buses[]` entries carry `"acdc":"DC"` and an explicit `"nodes":[0,1]` for DC
  buses; AC buses default to 4-node (a, b, c, n).
- `vscs[]` contains heterogeneous converter types (`acdc_3ph_4w_vdc_q`,
  `acdc_3ph_4w_pq`, and a `vsg` grid-former). The grid-former block nests its
  control parameters under `vsg:{...}`.
- `line_codes` defines per-conductor R/X matrices used by `lines[]` via
  `code:` references.

## Files

- `uds8bus.hjson` — network description (edit this to change topology).
- `xy_0.json` — Newton-Raphson initial guess; regenerate if it stops converging.
- `run.py` — primary entry: ini + 1 s run + annotated SVG output.
- `opt.py` — VSC reference sweep example.
- `uds8bus.svg` — hand-authored topology template (read by `model2svg`).
- `uds8bus_tooltips.svg` — generated output; do not edit by hand (overwritten).
- `legacy/` — older 7-bus version built with the numba/CFFI pipeline
  (`acdc_7bus_cffi`). Kept for reference; not used by the current workflow.

## Dependencies

Python ≥ 3.10 with `pydae-core`, `pydae-uds` (local checkout at the path
above), CasADi, NumPy. No C compiler required — the model is integrated via
CasADi/SUNDIALS.
