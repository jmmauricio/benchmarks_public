# benchmarks_public

Power system benchmark cases simulated with [pydae](https://github.com/jmmauricio/pydae)
using the CasADi backend (no C compiler required).

Each case follows the same three-function pattern:

| Function | What it does |
|---|---|
| `build()` | Assembles the DAE from the HJSON network description into a CasADi SX graph |
| `ini()` | Solves the load flow, prints a bus/generator report, validates against reference results, and runs small-signal analysis |
| `run()` | Time-domain simulation with a specific disturbance; saves a results plot |

## Systems

### IEEE 39-bus (New England)

Located in `ieee39/` (separate git repository).

- **39 buses, 10 synchronous generators** — the standard New England test system.
- Reference eigenvalues from PESTR18 (9 electromechanical modes).
- **`ini()`** — load-flow report, validation, small-signal analysis, eigenvalue plot
  (`ieee39_eig.png`) comparing computed modes against PESTR18 reference.
- **`run()`** — bolted short circuit at bus 16, cleared at t = 1.2 s; generator
  speeds and bus voltages plotted (`ieee39_run.png`).

```
ieee39/cases/base/
    ieee39.hjson   # network description + reference results
    main.py        # build / ini / run
    README.md
```

### NTS (REE España)

Located in `nts/`.

- **5-bus equivalent** of the REE España system with two synchronous generators.
- Reference step-response curve from the NTS REE document.
- **`ini()`** — load-flow report, validation, small-signal analysis including
  inter-area mode identification (0.1–0.5 Hz band).
- **`run()`** — +0.018 pu step in `v_ref_4` at t = 2 s, response to t = 40 s;
  `p_g_1` (MW) plotted against the NTS reference curve (`nts_base.png`).

```
nts/cases/base/
    nts_base.hjson        # network description + reference results
    nts_base_xy_0.json    # initial guess for the Newton-Raphson solver
    main.py               # build / ini / run
    README.md
```

## Requirements

- Python ≥ 3.10
- [`pydae-core`](https://github.com/jmmauricio/pydae) and `pydae-bps` installed
- `matplotlib`, `numpy`

No C compiler is needed — CasADi integrates the model with SUNDIALS/IDAS.

## Quick start

```bash
cd ieee39/cases/base
python main.py          # load flow + small-signal analysis + short-circuit simulation

cd nts/cases/base
python main.py          # load flow + small-signal analysis + v_ref_4 step simulation
```

Run a single stage:

```bash
python -c "import main; main.ini()"   # steady state + small-signal only
python -c "import main; main.run()"   # time-domain simulation only
```
