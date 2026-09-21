# Faster time-domain solving for pydae CasADi models

This note documents a performance investigation into `pydae` time-domain
simulation (`CasadiModel.run`) on a large model, the root cause, and **two new
solving strategies** that make it 30–60× faster with no C compilation:

1. **Grid integration** — drive the existing IDAS integrator over a *grid of
   output times in a single call* instead of one call per `Dt`.
2. **A fixed-step implicit-trapezoidal integrator** built from the model's
   CasADi Jacobians, with a frozen (reused) sparse factorization.

Everything below was measured on the Spanish pre-blackout benchmark
(`esp_preblackout.hjson`): **731 differential states, 981 algebraic variables,
~130 machines**, built with `BpsBuilder(..., use_casadi=True)` +
`CasadiBuilder`/`CasadiModel`. A reference implementation lives in `main.py`
(`run_grid`, class `Trapezoidal`).

---

## 1. The symptom

- `ini()` (load flow) and small-signal analysis: **fast** (~10 ms with a good
  `xy_0` seed) — a *single* nonlinear solve.
- `model.run(...)` (time-domain): **~90 s for 1 s of simulation** at `Dt=0.01`
  (≈ 900 ms per output step). Users tried larger `xy_0`/tolerances with no luck.

## 2. What it is **not** (measured, so we stop guessing)

Profiling `run(0.2)` (20 steps) with `cProfile`:

```
run(0.2) = 20 steps: 18.0 s -> 901 ms/step
  17.911 s (99.4%)  casadi._casadi.Function_call   (40 calls, 0.448 s each)
```

- **Not dense linear algebra.** The integrator uses `linear_solver: 'csparse'`
  (sparse). No `np.linalg.solve` on the hot path. (The dense `_newton_solve` in
  `casadi_model.py` is only an `ini` fallback.)
- **Not solver accuracy.** Loosening `reltol/abstol` (1e-6→1e-3): no change
  (~900 ms/step).
- **Not `calc_ic`.** Turning off IDAS consistent-IC recomputation: only ~17%
  (900→746 ms/step).
- **Not the integrator *object* being recreated.** `run()` creates the CasADi
  integrator once (only when `Dt` changes) and reuses the object.

Single SX-function evaluation costs (for scale): **residual ≈ 0 ms**, **full
Jacobian ≈ 69 ms** (interpreted SX graph for 1712 variables).

## 3. Root cause: BDF cold-start on every `Dt` call

`CasadiModel.run` (casadi_model.py ≈ line 518) loops:

```python
if self._current_dt != self.Dt:
    self.integrator = ca.integrator('idas_int', 'idas', self.dae_dict, 0.0, self.Dt, opts)
while self.t < t_end:
    res = self.integrator(x0=self.x, z0=self.y_run, p=p_run_vec)   # ONE call per Dt
    ...
```

A CasADi integrator is a **pure, stateless function** `(x0, z0, p) → (xf, zf)`.
It does **not** carry the SUNDIALS internal state (BDF history: past states,
current order, current internal step size, the factorized iteration matrix)
between calls. So every per-`Dt` call **cold-restarts** the BDF method: it resets
to **order 1**, picks a tiny conservative initial step, and re-does the entire
startup transient.

Measured IDAS statistics (`integrator.stats()`), integrating 10 output points of
`Dt=0.01`:

| driving pattern | internal `nsteps` | Jacobian factorizations `nlinsetups` |
|---|---:|---:|
| 10× per-`Dt` calls (pydae `run`) | **110** | **110** |
| 1 call over a 10-point grid | **15** | **15** |
| *(first per-`Dt` call alone, 0.01 s)* | 11 | 11 |

Each 0.01 s call burns ~11 internal steps and ~11 sparse Jacobian factorizations
(`nlinsetups` = evaluate + factorize — the expensive op) just to ramp BDF up,
then throws it all away for the next `Dt`.

**Why CasADi doesn't reuse the state:** keeping the integrator stateless makes it
composable / differentiable / embeddable. The sanctioned way to get intermediate
outputs from one continuous solve is the **grid of output times**.

---

## 4. Solution A — grid integration (`run_grid`)

Integrate to `t_end` in **one** IDAS call over a grid of output times; SUNDIALS
initializes once, ramps BDF up, keeps/reuses the factorization across internal
steps, and simply *reports* the state at each grid point (stepping past them and
interpolating).

```python
def run_grid(model, t_end):
    model._route_dict({}, phase='run')
    Dt = model.Dt
    n  = int(round((t_end - model.t) / Dt))
    tout = model.t + Dt * np.arange(1, n + 1)
    opts = {'print_stats': False, 'max_num_steps': 100000,
            'reltol': model._integrator_reltol, 'abstol': model._integrator_abstol,
            'linear_solver': 'csparse'}
    integ = ca.integrator('grid', 'idas', model.dae_dict, model.t, tout, opts)
    p = np.concatenate((model.p_vals, model.u_run_vals))
    res = integ(x0=model.x, z0=model.y_run, p=p)
    Xf, Zf = np.array(res['xf']), np.array(res['zf'])   # (n_x, n), (n_z, n)
    # append each column to model.Time/X/Y/Z, advance model.x/y_run/t ...
```

Use **once per continuous segment**; at a discontinuity (`set_value`) call
`model.recalculate_algebraics()` first so the algebraic ICs are consistent, then
`run_grid` the next segment.

**Measured:** `run(0.5)`: per-`Dt` loop **36.3 s** → single grid call **1.14 s**
(**32×**), identical result.

### Cold-start mitigation for *many short* segments

If you must chain many short calls (≤100 ms) with changes between them, each
still cold-starts. Make the restart cheap (measured on a repeated 100 ms / 10-pt
segment, integrator object reused across segments):

| options | ms / 100 ms segment | `nsteps` |
|---|---:|---:|
| default | 1199 | 15 |
| `+ step0 = Dt` | 403 | 4 |
| `+ step0 = Dt, calc_ic = False` | **241** | 4 |

- **`step0` (initial step-size hint ≈ `Dt`)** is the big lever: BDF starts at a
  sensible step instead of ramping from tiny — `nsteps` 15→4, `nlinsetups` 15→3.
  (The valid CasADi option name is `step0`, not `first_step`/`first_time`.)
- **`calc_ic = False`** skips the in-integrator consistent-IC solve (valid when
  the previous segment left consistent ICs; after a real discontinuity call
  `model.recalculate_algebraics()` once instead).
- **Reuse the integrator object** (create once with a *relative* grid
  `[Dt, 2Dt, …, Tspan]`, feed new `x0/z0/p`, track absolute time yourself) — saves
  the ~0.1 s re-creation per call.

---

## 5. Solution B — implicit trapezoidal integrator (`Trapezoidal`)

The right tool for **many short / interactive / real-time** stepping: fixed step,
**no adaptive restart → constant per-step cost**, and the iteration matrix is
factorized once and **reused** (modified Newton).

### Formulation

DAE from `model.dae_dict`: `dx/dt = f(x, z, p)`, `0 = g(x, z, p)` (`ode`/`alg`).
Trapezoidal rule, unknown `w = [x_{n+1}; z_{n+1}]`:

```
R(w) = [ x_{n+1} - x_n - (Dt/2)(f(x_n,z_n) + f(x_{n+1},z_{n+1})) ]   (n_x rows)
       [ g(x_{n+1}, z_{n+1})                                     ]   (n_z rows)
```

Iteration matrix (from **CasADi AD**, sparse):

```
J = dR/dw = [ I - (Dt/2) f_x    -(Dt/2) f_z ]
            [ g_x                 g_z        ]
```

### Construction & stepping

```python
d = model.dae_dict; x, z, p, f, g = d['x'], d['z'], d['p'], d['ode'], d['alg']
w  = ca.vertcat(x, z);  Dt = ca.SX.sym('Dt')
xn = ca.SX.sym('xn', nx);  fn = ca.SX.sym('fn', nx)          # x_n and f_n (frozen in a step)
R  = ca.vertcat(x - xn - (Dt/2)*(fn + f), g)
R_fn = ca.Function('R', [w, xn, fn, p, Dt], [R])
J_fn = ca.Function('J', [w, xn, fn, p, Dt], [ca.jacobian(R, w)])
f_fn = ca.Function('f', [x, z, p], [f])
```

Per step: `f_n = f_fn(x,z,p)`; modified-Newton on `w` (guess = previous `w`):
factorize `J_fn(...)` once and **reuse the factorization** for `refactor_every`
steps. **Refactor** when (a) scheduled, (b) after a `set_value` discontinuity, or
(c) adaptively when a step needs more than `refactor_iters` Newton iterations;
**retry once** with a fresh factorization if a step fails, else raise.

### Linear-solver backend (the sparse solve for `J dw = -R`)

The iteration matrix `J` has **fixed sparsity for the whole simulation**, so the
symbolic factorization is done once and only the numeric factorization is
repeated. Two backends:

- **`ca.Linsol['qr']` — CasADi-native, recommended.** Stays in CasADi `DM`, no
  `DM → scipy.csc` conversion. Correct reuse pattern (getting this wrong makes
  `solve` re-factorize every call!):
  ```python
  lin = ca.Linsol('lin', 'qr', J.sparsity())   # bind fixed sparsity
  lin.sfact(J_dm)                               # symbolic factorization — ONCE
  lin.nfact(J_dm)                               # numeric factorization — on refactor only
  dw = lin.solve(J_dm, -R_dm)                   # cheap back-sub — reuses nfact
  ```
- **`scipy.sparse.linalg.splu` — fallback / default.** `J_dm → scipy.csc → splu`;
  the factor object's `.solve()` reuses. Works, but pays a `DM→CSC` conversion.

Measured on this `J` (1712×1712, nnz 7308):

| backend | factorize | reuse-solve |
|---|---:|---:|
| `ca.Linsol['qr']` | **0.043 ms** | **0.014 ms** |
| `ca.Linsol['lapacklu']` | 3.24 ms | 1.00 ms |
| `ca.Linsol['csparse']` | 16.4 ms | 4.95 ms |
| scipy `splu` | 1.21 ms *(incl. convert)* | 0.055 ms |

**Plugin choice dominates within Linsol** — `'qr'` is fast; `'csparse'`/
`'lapacklu'` are 10–60× slower in this pattern. **Caveat:** at the *full step*
level the linear solve is a small fraction — `Linsol['qr']` and scipy give
comparable ms/step because a step is dominated by **CasADi `Function`-call
overhead + `DM ↔ numpy` marshaling per Newton iteration**, not the factorization.
So the `ca.Linsol` win is mainly *native / no scipy dep / no conversion* plus a
micro-speedup; the real end-to-end lever is cutting per-iteration Python
round-trips. See `main_dev.py::Trapezoidal` (both backends) and
`prompt_add_linsol.md`.

### Validation & benchmarks (3 s transient after a −1000 MW machine step, `Dt=0.02`)

| solver | time | per step | Newton its/step | max |x − IDAS| |
|---|---:|---:|---:|---:|
| IDAS grid (reference, tol 1e-9) | 1.08 s | — | — | — |
| trapezoidal, refactor every 1 | 0.41 s | 2.7 ms | 1.8 | 6.9e-07 |
| trapezoidal, refactor every 5 | 0.24 s | 1.6 ms | 1.9 | 6.9e-07 |
| **trapezoidal, refactor every 20** | **0.20 s** | **1.4 ms** | 1.9 | **6.9e-07** |

2nd-order accurate; freezing the factorization for 20 steps costs nothing in
accuracy and ~2× in speed. Constant per-step cost, no cold-start.

### End-to-end (the 25 s inter-area `run_link` simulation)

| solver | wall time |
|---|---:|
| pydae `model.run` (per-`Dt` IDAS) | ~90 s |
| `run_grid` (single IDAS grid call) | ~3.5 s |
| **`Trapezoidal` (frozen J)** | **~1.4 s** |

All pure CasADi + SciPy, no C compilation.

---

## 6. When to use which

| scenario | recommended |
|---|---|
| one-shot load flow / SSA | existing `model.ini` / `A_eval` (CasADi is great here) |
| one long continuous run | **grid integration** (IDAS accuracy, one call) |
| many short segments / discontinuities | **grid + `step0` + reuse**, or **Trapezoidal** |
| interactive / real-time / co-simulation | **Trapezoidal** (frozen J, constant cost) |
| need adaptive high accuracy / stiff, rare outputs | IDAS (adaptive BDF) as-is |

---

## 7. Reference code

`main.py` in this folder:
- `run_grid(model, t_end)` — grid IDAS integration for one continuous segment.
- `class Trapezoidal(model, tol, max_iter, refactor_every, refactor_iters)` with
  `.run(t_end)` and `.refactor()` — the fixed-step implicit solver.
- `warm_ini(model)` — warm-started load flow (a flat start diverges on this
  100-bus / 130-machine model); needed before either integrator.
- `run_link(...)` — uses `Trapezoidal` across a `set_value` discontinuity.

---

## 8. PROMPT — implement these solvers in pydae (separate session)

> **Task.** Add two faster time-domain integration strategies to
> `pydae.core.model.casadi_model.CasadiModel` (package `pydae-core`,
> `.../src/pydae/core/model/casadi_model.py`), motivated by the analysis in this
> file. Keep the CasADi (no-compile) pipeline; do **not** require cffi.
>
> **Background (verify first).** `CasadiModel.run` (≈ line 518) currently calls a
> CasADi `idas` integrator **once per `Dt`**. Because the integrator is a
> stateless function, each call cold-restarts BDF (order 1, tiny step): measured
> ~11 internal steps + 11 sparse Jacobian factorizations per 0.01 s, giving
> ~900 ms/step on a 731-state/981-algebraic model — ~90 s for 1 s of sim.
> Reproduce this with `integrator.stats()` (`nsteps`, `nlinsetups`) comparing a
> per-`Dt` loop vs a single grid call before changing anything.
>
> **Deliverable 1 — grid integration.** Replace/augment `run` so a continuous
> `run(t_end)` builds the integrator over a **grid of output times**
> `model.t + Dt*arange(1, n+1)` and calls it **once**, then appends each output
> column to `Time/X/Y/Z` and advances `x/y_run/t`. Preserve the exact public
> behavior of the current `run` (same storage layout, `get_values`, `post`,
> `h_dict`/`_h_fn` outputs, `_route_dict(update_dict, phase='run')` semantics).
> Expected ~30× speedup, identical trajectory. Reference: `run_grid` in `main.py`.
>
> **Deliverable 2 — fixed-step implicit integrators.** Add a solver selection to
> `run` (e.g. `run(t_end, method='idas'|'idas_grid'|'trapezoidal', Dt=..., **opts)`,
> default `'idas_grid'` for backward-compatible-but-fast, keeping `'idas'` = the
> old per-`Dt` loop for exact reproducibility). Implement **implicit
> trapezoidal** (and ideally also **implicit/backward Euler** as a robust option)
> built from `self.dae_dict` (`x, z, p, ode, alg`) via `ca.jacobian`:
> - residual `R = [x - x_n - (Dt/2)(f_n + f); g]`, matrix `J = dR/d[x,z]`;
> - modified Newton: factorize `J` and **reuse the factorization** for
>   `refactor_every` steps; refactor on schedule, after discontinuities, or
>   adaptively when a step exceeds `refactor_iters`; retry once with a fresh
>   factor on failure.
> - **Pluggable linear-solver backend** — `ca.Linsol['qr']` (native: `sfact`
>   once, `nfact` on refactor, `solve` reuses) **and** `scipy.splu` (kept). See
>   the dedicated, self-contained handoff prompt **`prompt_add_linsol.md`** for
>   the exact API, plugin caveats, tests, and the "marshaling overhead is the real
>   bottleneck" finding.
> Reference: `class Trapezoidal` in `main_dev.py` (both backends; validated to
> 6.7e-07 vs IDAS, ~1.4–2.8 ms/step frozen-J).
>
> **Deliverable 3 — cold-start knobs for IDAS.** Expose `step0` (initial step),
> `calc_ic`, and optional integrator-object reuse in the IDAS options (measured
> ~5× on chained 100 ms segments). Add a `recalculate_algebraics()` + refactor/IC
> reset hook to call after `set_value` discontinuities.
>
> **API / compatibility.**
> - Do not break existing scripts: default path must produce the same results as
>   today (allow `method='idas'` to select the legacy per-`Dt` loop).
> - Keep `set_value`, `recalculate_algebraics`, `post`, `get_values(s)`,
>   `_route_dict`, `Dt`, `_integrator_reltol/abstol` working unchanged.
> - The fixed-step solvers need a consistent starting point — document that
>   `ini`/`recalculate_algebraics` must be run first (see `warm_ini`).
>
> **Tests (add to `tests/`).**
> 1. Correctness: on a small model (e.g. the milano/ieee-style fixtures) and on a
>    large one, assert `max|x_grid − x_perDt|` and `max|x_trap − x_idas|` are
>    within tolerance (~1e-4..1e-6 at a resolving `Dt`).
> 2. Discontinuity handling: `set_value` + `recalculate_algebraics` + refactor,
>    then continue — trajectory matches a reference.
> 3. Determinism/energy: an undisturbed run stays at equilibrium (0 Newton iters,
>    `omega≈1`).
> 4. Benchmark (non-gating): report ms/step and speedup vs the per-`Dt` loop.
>
> **Edge cases to handle:** `n <= 0` spans; models without `h_dict`/`_h_fn`;
> algebraic-only or differential-only systems; singular/ill-conditioned `J`
> (fail → refactor/retry → clear error); `Dt` changes between calls; keeping the
> `Time/X/Y/Z` lists vs `post()` arrays consistent.
>
> **Docs.** Update the CasadiModel docstring and any `docs/` time-domain section
> with the method options and the "why the per-`Dt` loop was slow" explanation
> from §3 of this file.
