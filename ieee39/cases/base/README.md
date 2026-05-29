# IEEE 39-bus (New England) — base case

CasADi pipeline for the IEEE 39-bus / 10-generator New England test system.

- `ieee39.hjson` — network description: buses, lines, generators, loads, and the
  reference results (load-flow values and the 9 reference electromechanical modes)
  used for validation.
- `main.py` — three functions: `build()`, `ini()`, `run()`.

## Requirements

A Python environment with `pydae-core`, `pydae-bps` and `matplotlib` installed.
No C compiler is needed — the CasADi backend folds the model into an SX graph and
integrates it with SUNDIALS/IDAS.

## Running

Run both stages (initialization, then the time-domain fault):

```bash
python main.py
```

Run a single stage:

```bash
python -c "import main; main.ini()"   # steady state + small-signal only
python -c "import main; main.run()"   # time-domain short circuit only
```

## `ini()` — steady state and small-signal analysis

1. Builds the model and solves the load flow (`model.ini`).
2. Prints the **initialization report**: bus voltages/angles and generator dispatch.
3. Prints the **validation table**, comparing the solved state against the reference
   values in `ieee39.hjson` (each row marked OK/✓ when within tolerance).
4. Computes the reduced state matrix `A` (`model.A_eval()`) and runs **small-signal
   analysis** (`ssa.damp`): eigenvalues, damping ratios and oscillation frequencies.

**Output file:** `ieee39_eig.png` — an eigenvalue plot showing the computed modes
(blue circles) overlaid with the reference electromechanical modes (red crosses).

## `run()` — time-domain short-circuit simulation

`run(t_clear=1.2, bus_fault='16')`:

1. Initializes and runs 1 s of steady operation.
2. Applies a **short circuit** at `bus_fault` (default bus 16) by ramping the bus
   shunt conductance up, holds it through the fault, then ramps it back to clear the
   fault at `t_clear` (default 1.2 s).
3. Continues the post-fault response up to t = 10 s.

**Output file:** `ieee39_run.png` — two stacked plots: generator speeds (top) and the
voltages of buses 30–39 (bottom) versus time.
