# Prompt — add a CasADi `Linsol` linear-solver backend to the trapezoidal solver

**Task.** The fixed-step **implicit-trapezoidal** integrator already added to
`pydae` (`CasadiModel`, package `pydae-core`,
`.../src/pydae/core/model/casadi_model.py`) solves its Newton linear system with
**SciPy sparse LU** (`scipy.sparse.linalg.splu`). Add a **CasADi-native backend
using `ca.Linsol`** as a *selectable option*, and **keep the SciPy backend fully
working** (do not remove it). Make the linear solver **pluggable**.

A validated reference implementation of both the trapezoidal method and the
`ca.Linsol` usage is in `benchmarks_public/esp/pre_blackout/main_dev.py`
(class `Trapezoidal`, `linsol=` argument). Read it first.

---

## 1. Why (context)

The trapezoidal step solves `R([x,z]) = 0` with a modified Newton, reusing the
LU factorization of the iteration matrix `J = dR/d[x,z]` for `refactor_every`
steps. The **sparsity of `J` is fixed for the whole simulation**, so the symbolic
factorization can be done once.

`ca.Linsol` is the CasADi-native alternative to SciPy: it keeps everything in
CasADi `DM`, needs **no `DM → scipy.csc` conversion**, and can use several LU/QR
plugins. It supports the same factorize-once / solve-many pattern SciPy does.

## 2. The correct `ca.Linsol` reuse pattern (get this exactly right)

`ca.Linsol.solve(A, b)` **re-factorizes** unless you use the three-tier API. The
verified pattern is:

```python
lin = ca.Linsol('name', 'qr', J.sparsity())   # bind the FIXED sparsity pattern
lin.sfact(J_dm)                                #   symbolic factorization  -> ONCE (structure never changes)
# ...each time J is refreshed (refactor_every, or after a set_value discontinuity):
lin.nfact(J_dm)                                #   numeric factorization
# ...each Newton iteration:
dw = lin.solve(J_dm, -R_dm)                    #   cheap back-substitution, REUSES nfact
```

Measured on the 1712x1712 (nnz 7308) esp benchmark: `nfact ≈ 0.043 ms`,
`solve ≈ 0.014 ms` (solve/nfact ratio ~0.33 → genuinely reusing, not
re-factorizing). Compare SciPy `splu`: factorize+convert ≈ 1.21 ms, solve ≈
0.055 ms.

**Plugin choice matters a lot.** Use **`'qr'`** by default. In this pattern
`'csparse'` (~165 ms/step) and `'lapacklu'` (~36 ms/step) were **10–60x slower**
than `'qr'` (~2.8 ms/step). Detect availability at runtime and fall back
gracefully (some plugins like `ma27`/`ma57` need HSL and may be missing).

## 3. What to build

**A small pluggable linear-solver interface** with two implementations, e.g.:

```python
class _LinBackend:                      # scipy
    def refresh(self, J_dm): ...        # (re)factorize from a CasADi DM Jacobian
    def solve(self, b): ...             # solve using the stored factorization, return ndarray/DM

class _ScipyLU(_LinBackend):            # existing behaviour, keep as default or option
    # J_dm.sparsity().get_triplet() + nonzeros() -> scipy.csc -> splu; solve reuses factor
class _CasadiLinsol(_LinBackend):       # NEW
    # ca.Linsol(name, plugin, sparsity); sfact ONCE; nfact on refresh; solve reuses
```

Wire a **`linsol` option** into the trapezoidal solver / `CasadiModel.run(...)`,
e.g. `linsol='scipy'` (default, unchanged behaviour) or `linsol='qr'` (or a
generic `linsol='casadi:<plugin>'`, default plugin `'qr'`). Keep
`refactor_every`, `refactor_iters`, `tol`, `max_iter` working for both backends.

The Newton loop is backend-agnostic:
- factorize (`backend.refresh(J_dm)`) when scheduled / forced;
- iterate: eval `R`, check `‖R‖`, `dw = backend.solve(-R)`, `w += dw`.

## 4. Constraints / compatibility

- **Do NOT remove or change the SciPy path's numerics.** Existing scripts that
  use the current trapezoidal solver must produce identical results. SciPy stays
  the default unless you make `qr` default *and* prove bit-parity in tests.
- Preserve the public surface: `set_value`, `recalculate_algebraics`
  (must still force a refactor), `post`, `get_value(s)`, `Time/X/Y/Z` storage,
  `Dt`, `_route_dict(update_dict, phase='run')`.
- `sfact` **once** at solver construction (sparsity is fixed); `nfact` on refresh
  only. Do not call `sfact` per step.
- Handle plugin-not-available: try the requested `ca.Linsol` plugin at
  construction; if it raises, fall back to SciPy (with a warning) so runs never
  hard-fail on a missing HSL/plugin.
- Keep the CasADi `DM ↔ numpy` boundary minimal (see §6).

## 5. Tests (add under `tests/`)

1. **Parity**: on a small fixture and on a larger model, trapezoidal with
   `linsol='scipy'` vs `linsol='qr'` — assert `max|x_scipy - x_qr|` ≤ ~1e-9 over
   a transient (they solve the same equations).
2. **Reference accuracy**: both vs an IDAS grid solve — `max|x - x_idas|` within
   ~1e-6 at a resolving `Dt` (reference: 6.7e-7 on the esp benchmark).
3. **Reuse check**: assert `Linsol.solve` time ≪ `nfact` time (factorization is
   actually reused), and that `sfact` is called once.
4. **Discontinuity**: `set_value` + `recalculate_algebraics` + forced refactor,
   continue — trajectory matches a reference.
5. **Fallback**: requesting an unavailable plugin falls back to SciPy without
   error.
6. **Benchmark (non-gating)**: report ms/step per backend/plugin.

## 6. Set expectations (measured)

The linear solve is a **small fraction of a trapezoidal step**: at the step level
`Linsol['qr']` and SciPy `splu` are comparable (~1.4–2.8 ms/step) because the step
is dominated by **CasADi `Function`-call overhead + `DM ↔ numpy` marshaling per
Newton iteration**, not the factorization. So the `ca.Linsol` win is mainly
*cleanliness / no SciPy dependency / no format conversion / stays in DM*, plus a
small micro-speedup — not a large end-to-end gain. If you want a real end-to-end
speedup, the follow-up lever is **reducing per-iteration Python round-trips**
(evaluate `R`, `J`, and the solve in as few native CasADi calls as possible, e.g.
a purpose-built Newton/step `Function`), which is out of scope here but worth a
note in the docstring/docs.

## 7. Deliverable

- New `ca.Linsol` backend + pluggable selection, SciPy retained and default-safe.
- `linsol=` option documented in the `CasadiModel.run` / solver docstring, with
  the `sfact/nfact/solve` note and the `'qr'` recommendation.
- Tests above. Reference implementation: `main_dev.py::Trapezoidal`.
