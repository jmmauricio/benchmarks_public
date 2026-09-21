# NTS base case

Pydae pipeline for the NTS-REE base case.

- `nts_base.hjson` — network description: buses, lines, generators, loads, and the
  reference results (load-flow values and the `p_g_1` response to a `v_ref_4` step)
  used for validation.
- `nts_base_xy_0.json` — saved initial guess that seeds the Newton-Raphson solver.
- `main.py` — three functions: `build()`, `ini()`, `run()`.

## Requirements

A Python environment with `pydae-core`, `pydae-bps` and `matplotlib` installed.
No C compiler is needed — the CasADi backend folds the model into an SX graph and
integrates it with SUNDIALS/IDAS.

## Running

Run both stages (initialization, then the time-domain step response):

```bash
python main.py
```

Run a single stage:

```bash
python -c "import main; main.ini()"   # steady state + small-signal only
python -c "import main; main.run()"   # time-domain v_ref_4 step only
```

## `ini()` — steady state and small-signal analysis

1. Builds the model, applies the line 2–3 impedance adjustment, and solves the load
   flow (`model.ini`).
2. Prints the **initialization report**: bus voltages/angles, generator dispatch and
   line flows.
3. Prints the **validation table**, comparing the solved state against the reference
   values in `nts_base.hjson`.
4. Computes the reduced state matrix `A` (`model.A_eval()`) and runs **small-signal
   analysis**: `ssa.damp` (eigenvalues, damping ratios and frequencies), `ssa.eig`
   (eigenvectors and participation factors), and `ssa.get_mode` (inter-area modes in
   the 0.1–0.5 Hz band).

No output file — results are printed to the console.

## `run()` — v_ref_4 step-change simulation

1. Initializes and runs 2 s of steady operation.
2. Applies a **+0.018 pu step** to `v_ref_4` (voltage reference of generator 4) and
   integrates for 40 s.

**Output file:** `nts_base.png` — two stacked plots: generator 1 speed `omega_1`
(top) and the active power `p_g_1` (MW) compared against the NTS reference curve
(bottom).
