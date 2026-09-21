# Spanish Power-System pydae Equivalent — End-to-End Modeling Guide

**Goal.** Build a dynamic (RMS phasor) model of the Spanish transmission system in
[`pydae`](https://github.com/pydae) / `pydae-bps` (`.hjson` format), representing the
operating conditions just **before the 28 April 2025 Iberian blackout**, and use it for
small-signal analysis (SSA), inter-area oscillation studies, and overvoltage studies.

The pipeline is:

```
PyPSA-Eur (Europe)  ─►  PyPSA-Spain (fork)  ─►  single pre-blackout snapshot (PF)
        ─►  pypsa_to_pydae.py  ─►  esp_preblackout.hjson (pydae-bps)
        ─►  dynamic refinements (inertia spread, PSS tuning, WECS wind)
        ─►  analyses: ini/SSA, run_link (inter-area), run_overvoltage
```

This document explains **every step**, the **decisions** taken and **why**, the
**modifications** made to both PyPSA and pydae, the **current status**, and how to
**reproduce it from scratch**.

---

## 0. Glossary

| Term | Meaning |
|---|---|
| **VRE** | Variable Renewable Energy — solar PV + wind (non-synchronous, inverter-based). |
| **NUTS3** | Eurostat statistical region level 3 (≈ province). Used as the bus granularity. |
| **PF** | Power flow (load flow). |
| **SSA** | Small-signal analysis (eigenvalues of the linearized DAE). |
| **Inter-area mode** | Low-frequency (~0.2 Hz) electromechanical oscillation of Iberia swinging against continental Europe. |
| **GENROU / GENSAL** | Round-rotor / salient-pole synchronous machine models. |
| **AVR / gov / PSS** | Automatic Voltage Regulator / turbine governor / Power System Stabilizer. |
| **AGC** | Automatic Generation Control (secondary frequency control / slack). |
| **Ferranti effect** | Voltage rise on a lightly-loaded long transmission line (shunt capacitance dominates). |
| **RMS phasor model** | Positive-sequence, fundamental-frequency dynamic model (captures electromechanical + control dynamics and reactive-balance voltage behaviour; not sub-cycle EM transients). |

---

## 1. Environments

Two independent toolchains are used; keep them separate.

### 1.1 PyPSA-Spain (the grid-building side)
- Repo: `/Users/jmmauricio/workspace/pypsa/pypsa-spain` — a **fork of PyPSA-Eur**.
- Package manager: **pixi** (preferred) or conda-lock; conda env name **`pypsa-eur`**.
- Install (macOS ARM): `make install-lock-macos-arm` (targets exist for `linux64`, `windows`, `macos64`).
- Run tasks: `pixi run <task>` or `conda run -n pypsa-eur <cmd>`.
- Snakemake workflow; **always pass the Spain config**:
  ```bash
  snakemake all --configfile config/config_ES.yaml --cores 4
  ```

### 1.2 pydae (the dynamic-modeling side)
- Repo: `/Users/jmmauricio/workspace/pydae` (mono-repo: `pydae-core`, `pydae-bps`).
- conda env **`pydae_dev`** has `pydae` + `hjson` + `casadi`.
- The conda **`base`** env has `geopandas` + `plotly` (used for the maps) but **not** `hjson`
  → the plotting scripts load the model with plain `json` (the `.hjson` we emit is valid JSON).
- Set `KMP_DUPLICATE_LIB_OK=TRUE` when running the RAG in `base`.

### 1.3 Working folder for the pydae model
`/Users/jmmauricio/workspace/benchmarks_public/esp/pre_blackout/` — holds the `.hjson`,
the analysis scripts (`main.py`, `export_pf.py`, `plot_sld*.py`), the init seed
(`xy_0.json`), and all outputs.

---

## 2. PyPSA-Eur → PyPSA-Spain

**PyPSA-Eur** is the open European sector-coupled energy-system model (a Snakemake workflow
that retrieves data, builds a network, and solves an optimisation). **PyPSA-Spain** is a fork
that specialises it for the Spanish system and adds Spain-specific data/features while
tracking upstream (so most code is unmodified upstream and easy to re-merge).

Key conventions of the fork (see `pypsa-spain/CLAUDE.md`):
- **Config layering** (each overrides the previous): `config.default.yaml` →
  `plotting.default.yaml` → `config.yaml` (local) → **`config_ES.yaml`** (loaded last as a
  directive → always wins; holds the `pypsa_spain:` block).
- Spain-specific edits inside upstream files are flagged with `#####` markers and
  `[PyPSA-Spain]` log prefixes so they survive re-merges.
- The `pypsa_spain:` block toggles Spain features: `q2q_transform`, `electricity_demand`
  (NUTS3 profiles), `interconnections`, `update_elec_capacities` (esios real capacities),
  `ISA_class`, `H2_valley_demands`, `H2_imports_exports`, `regional_network_focus`,
  `pop_layouts_HR`, `industry_scenario`.
- We use the **sector-coupled** network (not electricity-only).

---

## 3. Spatial resolution — the NUTS decision

The dynamic model needs **buses**. PyPSA-Eur normally clusters the grid to an arbitrary
number of nodes; PyPSA-Spain can cluster to **administrative (NUTS) regions**.

**Decision: use NUTS3 (`adm`) regions as buses.** Rationale:
- One bus per province (≈50 peninsular NUTS3 regions) gives a **geographically meaningful,
  reproducible** equivalent — generation, load and interconnections map to real places.
- It is a good compromise between fidelity and a tractable dynamic model (~50 electrical
  buses, before adding machine terminal buses).
- The NUTS3 geometry (`resources/ES_test/nuts3_shapes.geojson`, index = NUTS3 code such as
  `ES512`) is reused for all maps; a copy lives next to the pydae model
  (`nuts3_shapes.geojson`) so the plotting scripts are self-contained.

Buses are named by their NUTS3 code (`ES111`, `ES512`, …); France is a single external node
`FR`.

**Single-snapshot power flow, not the full year.** Solving the whole 8760-hour optimisation
is unnecessary (and slow). For a dynamic snapshot we take a **pre-solved network**, apply a
**single-snapshot economic dispatch**, and run `n.pf()` (AC power flow) for that instant.
This avoids the ~20-min full-year LP.

---

## 4. The power-system RAG (blackout knowledge base)

The official incident report is indexed in a retrieval-augmented store:
- Location: `/Users/jmmauricio/workspace/mdbib_bps/` (source markdown in `mds/…Final Report
  on the Grid Incident in Spain and Portugal on 28 April 2025…`).
- Query (in `base`, needs `KMP_DUPLICATE_LIB_OK=TRUE`):
  ```bash
  cd /Users/jmmauricio/workspace/mdbib_bps/rag
  KMP_DUPLICATE_LIB_OK=TRUE conda run -n base python rag.py "your question"
  ```
- It returns reranked passages with sources. Used throughout to **ground modeling choices in
  the real event** (inter-area frequency, overvoltage mechanism, protection behaviour).

---

## 5. The 28 April 2025 Iberian blackout (what the report says)

Facts used to calibrate/validate the model (from the RAG):
- A **~0.2 Hz inter-area oscillation** (the "East-Centre-West continental mode", Iberia vs
  continental Europe) appeared ~12:19–12:22. Operators damped it by **cutting Spain→France
  exports, coupling southern lines, and changing the FR-ES HVDC mode**, which **raised
  Iberian voltage**.
- A cascade of **overvoltage-driven generation trips** followed: units disconnected →
  **loss of reactive-power absorption → steep voltage rise** (Finding #7). Example: Granada
  400 kV rose **400.6 → 417.2 kV in ~44 s**; its reactive went from **−165 → −96 Mvar**
  (absorbing less as it rose). Many **overvoltage protections were mis-set** (below the
  1.1–1.2 pu ride-through requirement, some with **zero time delay**), tripping prematurely.
- **Shunt reactors** were being switched by operators; **Palmela's reactor tripped at 12:19
  on _under_voltage** (379.8 kV, 380 kV/2 s setting). Crucially, "a substantial reactive
  power capacity from shunt reactors was **available but not activated during the voltage
  rise**" (manual action, lead time).
- **Loss of synchronism at 12:33:19.6**; **AC FR-ES lines tripped at 12:33:21.5**; system
  collapse ~12:33:30.

**Modeling implications:** the event is fundamentally a **reactive-balance / voltage
phenomenon** on a **low-inertia, high-VRE** grid, plus a **weakly-damped 0.2 Hz inter-area
mode** — both are within reach of an RMS phasor model.

---

## 6. Building the pre-blackout PyPSA snapshot

Modifications applied to the PyPSA-Spain network to represent the pre-blackout instant
(Monday 28 April 2025, around midday — high solar, moderate wind, exporting to France):

1. **Pick the snapshot** (midday, pre-blackout conditions).
2. **Model the interconnection as an external slack.** Decision: represent **France as an
   external slack node** (`FR`) rather than a PQ tie — it provides the continental voltage/
   frequency reference the Iberian system swings against. Portugal and France
   interconnections considered when matching flows.
3. **Nuclear must-run floor** set to the pre-blackout nuclear output.
4. **Match the real generation mix** (per carrier) and the **loads**; fine-tune (relaxing the
   nuclear limit where needed) so the dispatch resembles the reported mix.
5. Solve `n.pf()` for that snapshot → bus voltages/angles + per-generator P/Q.

Helper: `preblackout_fr_slack.py` (`build()` returns the France-as-slack solved network).

---

## 7. PyPSA → pydae conversion (`pypsa_to_pydae.py`)

The converter (in the PyPSA-Spain repo) turns the solved snapshot into a pydae-bps `.hjson`.
Key mappings and decisions:

- **Voltage base:** the electrical model is collapsed to a **single 400 kV layer**
  (PyPSA-Spain simplifies to 380 kV internally; line R/X are scaled by `(400/380)²` and the
  nominal set to **400 kV** to match the real EHV grid).
- **Synchronous machines:**
  - **GENROU** for thermal + nuclear (`CCGT_*`, `nuclear_*`, `biomass_*`).
  - **GENSAL** for hydro and run-of-river (`hydro_*`, `ror_*`).
  - Each machine sits on its **own terminal bus** behind a **step-up transformer** to the
    region bus (this keeps the DAE square when several machines share a region and lets each
    have its own AVR).
- **Non-synchronous (VRE):** `pv_pq_ss` grid-following inverters at **unity power factor**
  (`q_s_ppc = 0`), named `solar_*` and `onwind_*`. (Wind later migrated to a WECS model — §12.)
- **Controls (as in the ieee39 sample):**
  - `avr: sexs` (K_a=100, E_min=0, E_max=5) — note the **E_min=0** floor: the key
    reactive-absorption limit.
  - `gov: tgov1`; `pss: pss_kundur_2`.
  - PV/PQ machines get a **load control** `lc: {K_i, p_c_lc}` (active-power set-point).
  - The **slack / AGC reference** is set in the `agc` block on **`SLACK_FR`** (France), a
    large `genrou` (S_n=50 GVA, H=40) that motors to absorb the Spanish export.
- **Shunt compensation:** the inverters are modelled unity-PF, so the reactive they carried in
  the PyPSA solution is re-injected as a **shunt at each bus** (`X_pu = −100/q`), letting the
  synchronous machines relax to physical reactive. This is why the model has ~44 shunt
  **reactors** (absorbing the line charging) + a couple of capacitors.
- **`results` block:** the converter stores the raw PyPSA per-device `P_MW`/`Q_Mvar` and bus
  `V_pu/theta_deg` **for reference/warm-start**. ⚠️ The raw PyPSA reactive on small devices is
  unphysical (e.g. a 16 MVA biomass unit at −897 Mvar) — see §11.

Result: **`esp_preblackout.hjson`** — 101 buses (≈50 region + machine-terminal buses),
89 lines, 54 transformers, 55 synchronous machines, 77 VRE plants, 46 loads, 46 shunts.

---

## 8. The pydae model (`.hjson`) and how it is built

Top-level sections: `system`, `agc`, `buses`, `loads`, `lines`, `transformers`, `shunts`,
`syns`, `pvs` (+ `wecs` after §12), `sources`, `results`.

Build with the **CasADi backend** (no C compilation):
```python
from pydae.bps import BpsBuilder
from pydae.core.builder.casadi_builder import CasadiBuilder
from pydae.core.model.casadi_model import CasadiModel
grid = BpsBuilder("esp_preblackout.hjson", use_casadi=True); grid.uz_jacs = False
grid.construct("esp_preblackout")
model = CasadiModel(CasadiBuilder(grid.sys_dict).build())
```
- `use_casadi=True` folds the symbolic DAE into a CasADi SX graph.
- `uz_jacs=False` skips the u/z Jacobians (only the A matrix is needed for SSA).

---

## 9. Initialization — the `xy_0.json` seed

Flat start diverges on a ~100-bus / 130-machine network, so init is **seeded**:
```python
model.ini({}, xy_0="xy_0.json")   # simple, no hand-rolled Newton
```
- **`grid.construct()` auto-writes `xy_0.json` and `<name>_xy_0.json`** from the model's
  `xy_0_dict` init hints (in `pydae-bps/.../bmapu_builder_line_exp.py`) **on every build**.
  → **Never hand-edit `xy_0.json`** — it is overwritten by `build()` before `ini()` reads it.
  If `ini` fails, **fix the model's `xy_0_dict` init hints**, not the file (this is exactly
  what was needed for the WECS wind model — §12).
- The seed stays valid after **dynamic-only** parameter changes (H, reactances, PSS gains) —
  those don't move the load flow.

Analyses that follow `ini()`:
- **`report_buses` / `report_gens`** print the solved operating point.
- **SSA:** `model.A_eval()` builds `A = Fx − Fy·inv(Gy)·Gx`; `ssa.damp(A, model=model)`
  gives eigenvalues / damping / frequencies; `ssa.plot_eig(...)` plots the s-plane.

---

## 10. Solvers added to pydae (time-domain performance)

`CasadiModel.run` originally called the CasADi `idas` integrator **once per output step**,
cold-restarting BDF each time → ~90 s for 1 s of simulation. Fixes (all pure CasADi/SciPy):
1. **Grid integration** — one `ca.integrator(..., tout_grid)` over all output times (~32×).
2. **Implicit trapezoidal fixed-step solver** — Jacobian from `ca.jacobian`, sparse Newton
   with a **frozen/reused factorization** (modified Newton), `refactor_every` / `refactor_iters`
   guards. ~1.4 ms/step, matches IDAS to ~1e-6, no cold-start.
3. **`ca.Linsol` linear backend** — factorize-once (`sfact`) / refresh (`nfact`) / reuse
   (`solve`); the `'qr'` plugin is fastest.
Reference implementation: `main_dev.py` (`class Trapezoidal`). Write-ups: `new_solvers.md`,
`prompt_add_linsol.md`. The trapezoidal solver is what `run`, `run_link`, and `run_overvoltage`
use (`method='trapezoidal'`), with discrete events applied via `set_value` +
`recalculate_algebraics` between segments.

---

## 11. Dynamic-model refinements

### 11.1 Break the biomass inertia degeneracy
The converter emitted **10 biomass GENROU machines with identical parameters** (H=3.5, generic
reactances). Identical machines form a **degenerate cluster** (~1.75 Hz) that masqueraded as
the "critical mode" and hid the true inter-area mode. Fix — **differentiate dynamic parameters
only** (load flow unchanged):
- `H` by carrier + unit MVA size (nuclear 5.5, CCGT 5.0, biomass 3.0, hydro 3.5 base), scaled
  0.9–1.3× by size, plus ±15 % per-machine jitter.
- `X_d, X_q, X1d, T1d0` given ±12 % per-machine spread.
- Every machine stays at its real bus → **geography preserved**.

After this the true **Spain-vs-Europe inter-area mode surfaces at ~0.2 Hz** and, at the base
tie, is critical (matching the report) — 100 % of Spanish machines coherent, in anti-phase
with France.

### 11.2 PSS tuning (stabilizing the 0.2 Hz mode)
All 55 machines originally carried an **identical generic `pss_kundur_2`** (K_stab=1,
T1/T2=1/0.05, T3/T4=3/0.5) providing ~91° of phase lead at 0.2 Hz — **over-compensating** and
*feeding* the instability. Fix — retune the **10 largest units** (geographically spread:
`CCGT_ES620`, `nuclear_ES514`, `CCGT_ES213`, `nuclear_ES424`, and the big hydros):
- **~30° lead centered at 0.2 Hz** (T1=T3≈1.01, T2=T4≈0.59), washout T_w=10, **K_stab = 10**.
- Result at full-load base: **inter-area mode 0.202 Hz, ζ = +15.5 %**; system-wide minimum
  damping **+6.8 %** (a ~1 Hz local mode); **0 unstable modes**.
- ⚠️ `K_stab` must be ≥ 8 at full load; at K_stab=2 the mode is unstable (−24 %). Keep K=10.

### 11.3 Reactive-power corrections (display/reference)
The `results` block held **raw PyPSA reactive** (line-charging dumped onto whatever device sat
at a bus → 16 MVA biomass at −897 Mvar, etc.). The **physical** values come from the **solved
pydae model**: `p_g/q_g` for machines, `p_s/q_s` for inverters. `export_pf.py` writes those to
`esp_preblackout_pf.json` (+ a markdown report). In the solved model PV/wind are unity-PF
(Q≈0) and machines relax to physical reactive (biomass ≈ +1 Mvar, CCGT_ES620 ≈ +318 Mvar).

---

## 12. Wind: `pv_pq_ss` → `pmsm_pq_ss` (WECS)

The 35 wind units were migrated from the algebraic `pv_pq_ss` (in the `pvs` section) to the
**PMSM full-converter WECS model `pmsm_pq_ss`** in a new **`wecs` section** (dispatched by
`add_wecs`; `add_pvs` does **not** handle it, so wind must live in `wecs`, solar stays in
`pvs`). Each turbine adds **7 dynamic states** (wind-speed ramp, drivetrain torsion, turbine &
rotor speed, pitch integrator, pitch, MPPT power) → model went **731 → 976 states**.

Decisions / fixes required:
- **Setpoint semantics:** `p_s_ppc` behaves like the PV set-point — **output = min(p_s_ppc,
  available wind power)**. With `nu_w_0 = 16 m/s` (over nominal) the available power
  (~K_mppt3·1.2³ ≈ 0.69 pu) exceeds every `p_s_ppc` (≤ 0.34), so **output = p_s_ppc exactly**
  and the load flow is preserved (verified 0.00 MW error, Q=0).
- **CasADi build:** the model initially hard-coded `import sympy`; it was ported to the
  builder's `grid.backend` API so it emits CasADi SX (like `pv_pq_ss`).
- **Init fix (in the model's `xy_0_dict`):** curtailed operation needs a **large pitch angle**
  (~22–26°) to shed the over-nominal wind down to the low `p_s_ppc`. The original linear guess
  `beta = 3·(nu_w_0−11) = 15°` was too low → rotor runaway → singular Newton. Fixed by
  **inverting the c_p aero curve** for `p_w(beta) == p_seed` in step 20 of `pmsm_pq_ss.py`.
- **Post-migration SSA (base, K_stab=10):** still **0 unstable, inter-area 0.202 Hz / +15.5 %,
  min damping +6.8 %**, plus new well-damped (**+44 %**) wind drivetrain/pitch modes at ~0.1 Hz.

Tooling updated for the `wecs` section: `export_pf.py`, `plot_sld.py`, `plot_sld_interactive.py`.

---

## 12b. VRE Q–V droop — the faithful generation-trip overvoltage actuator

The base VRE inverters are unity-PF (`q_s_ppc = 0`), so in the original model **tripping
generation caused _under_voltage** (removing P), the *opposite* of the report. Real grid-code
VRE instead **absorbs reactive at high voltage** (Granada PV at −165 Mvar) — so tripping it
removes *absorption* → **overvoltage**. To capture this an optional **Q–V droop** was added to
`pv_pq_ss` (pydae-bps `pvs/pv_pq_ss.py`):

$$q_{s,ppc}^{\text{eff}} = q_{s,ppc} + K_{qv}\,(V_{qref} - V_m)$$

- New parameters **`K_qv`** (droop gain, default **0** → identical to the old model) and
  **`V_qref`** (reference voltage, default 1.0). Both are runtime-settable per unit
  (`K_qv_<name>`, `V_qref_<name>`) and a `q_s_ppc_droop_<name>` output is exposed.
- The existing soft current-saturation (`I_max`) naturally bounds the reactive the inverter
  can absorb, so no extra limiter is needed.
- **Load-flow preservation:** seed `V_qref` to each inverter's *equilibrium* terminal voltage
  (helper `enable_qv_droop` in `main.py`). Then the droop term is 0 at the operating point, so
  the seed / `xy_0.json` init is unchanged — verified **max Δq_s ≈ 7e-9 pu, max ΔV ≈ 4e-9 pu**
  at K_qv=2; SSA still **0 unstable** (min damping 6.84 % → 5.75 %, mildly reduced by the added
  voltage-feedback loop). No re-init needed (same trick as PSS tuning).

This makes VRE a **reactive absorber** at elevated voltage, so tripping it now removes
absorption → overvoltage, matching the report's actual actuator (§13 `run_gentrip_cascade`).

---

## 13. Analyses / tests available (`main.py`)

- **`ini()`** — steady-state init + SSA; writes `esp_preblackout_eig.{png,csv,html}`.
- **`run()`** — short trapezoidal time-domain run from the saved seed.
- **`run_link(dP_pu, machine, …)`** — steps the Spain→France interconnection power and shows
  the **0.2 Hz inter-area response**; overlays tuned vs generic-PSS (the latter runs away,
  reproducing the growing 0.2 Hz oscillation → breakdown, like the real event).
- **`run_overvoltage(trip_frac, …)`** — single-event **overvoltage** study: trips shunt
  **reactors** (`b_shunt → 0`) → reactive surplus → voltages rise; AVRs cut excitation toward
  `E_min=0` and saturate.
  - 2 reactors → settles at ~1.03 pu (stable rise).
  - ≥4 reactors → **growing overvoltage oscillation, AVR rail-to-rail saturation, voltage
    collapse** (~9 s) — the divergence is genuine instability, caught and plotted.
- **`run_gentrip_cascade(...)`** — the **faithful generation-trip overvoltage cascade**
  (Finding #7). Enables the VRE **Q–V droop** (§12b) on the 42 solar units, seeds an overvoltage
  by tripping a few shunt reactors, then applies **overvoltage protection with a definite-time
  delay** `t_delay` (a unit trips only after its bus V stays above the setting for `t_delay` s —
  the report faulted the *zero-delay* settings). Each trip removes reactive **absorption** →
  voltage rises → more units time out → **cascade**. Overlays two settings (`compare=True`):
  - **compliant 1.10 pu** → the droop holds the disturbance, **no trips**, settles.
  - **mis-set 1.045 pu** (below the requirement) → **11 solar trip at ~1.6 s**, then a slow
    growing voltage swing to **collapse at ~11 s**.
  Defaults are softened for a readable trace (`K_qv=1`, seed 4 reactors, `t_delay=0.6 s`); a
  gentler mis-set (`V_trip=1.05`) gives a **partial cascade that arrests** (7 trips, settles
  ~1.084 pu). Plot: `esp_preblackout_gentrip.png`. The §16-#1/#2 milestone, done.
- **`analyze_syn_q.py`** — records the **synchronous-machine reactive power** (`q_g·S_n`) and AVR
  field voltage through the cascade. Shows the reactive-balance shift: the 55 machines start net
  **supplying +1493 Mvar** (Spain) / France slack +544 Mvar; as solar trip and voltage rises the
  sexs AVRs cut field voltage toward the **E_min=0 floor** and the machines swing to net
  **absorbing** (Σ Spain ≈ −2100 Mvar in the arresting case; big units like hydro_ES432 +37 →
  −729, nuclear_ES514 +351 → −244, CCGT_ES620 +318 → +4 Mvar). Plot:
  `esp_preblackout_gentrip_synq.png`.

**Important reactive-posture finding (for triggers):** in this model the synchronous machines
are net reactive **suppliers** (+1493 Mvar) and the **shunt reactors are the absorbers**;
VRE is unity-PF. Therefore **tripping generation causes _under_voltage** here (removing supply),
while **tripping reactors causes overvoltage** (removing absorption). The report's "generation
trip → loss of absorption → overvoltage" maps to this model **once the VRE carry a Q–V droop**
(§12b) — then they too are reactive absorbers (like Granada at −165 Mvar) and tripping them
removes absorption → overvoltage. This is now implemented and driven by `run_gentrip_cascade`.

---

## 14. File inventory (working folder)

| File | Purpose |
|---|---|
| `esp_preblackout.hjson` | The pydae-bps model (base). |
| `esp_preblackout_dev.hjson` | Active dev model (**`main.py` uses this**): differentiated inertia, K_stab=10 PSS, `results` reactive cleaned, **wind as `wecs`/`pmsm_pq_ss`**. |
| `xy_0.json` / `esp_preblackout_xy_0.json` | Init seed (auto-written by `construct()`). |
| `main.py` | `build`, `ini`, `run`, `run_link`, `run_overvoltage`, `enable_qv_droop`, `run_gentrip_cascade`. |
| `verify_qvdroop.py` | Checks the Q–V droop preserves the equilibrium + SSA (§12b). |
| `analyze_syn_q.py` | Records synchronous-machine reactive power + AVR field voltage through the cascade. |
| `main_dev.py` | Reference trapezoidal + `ca.Linsol` solver harness. |
| `export_pf.py` | Solves PF, writes `esp_preblackout_pf.json` + `esp_preblackout_pf.md`. |
| `plot_sld.py` | Static geographic single-line diagram (SVG+PNG, hover `<title>` tooltips). |
| `plot_sld_interactive.py` | Interactive Plotly SLD (`esp_preblackout_sld.html`): size-by S_n/P/\|Q\|/\|I\|, label toggle, legend filtering, region-voltage choropleth. |
| `nuts3_shapes.geojson` | NUTS3 geometry (local copy, for the maps). |
| `new_solvers.md`, `prompt_add_linsol.md` | Solver write-ups / porting prompts. |
| `esp_preblackout_eig.{png,csv,html}` | SSA outputs. |
| `esp_preblackout_link.png`, `esp_preblackout_overvoltage.png`, `esp_preblackout_gentrip.png` | Time-domain outputs. |

In the PyPSA-Spain repo: `pypsa_to_pydae.py` (converter), `preblackout_fr_slack.py`
(France-slack snapshot), `pydae_init.py` / `pydae_ssa.py` (older init/SSA helpers,
superseded by the saved-seed workflow).

---

## 15. Current status

- ✅ Pre-blackout snapshot converted to a **976-state** pydae-bps model, **stable at full
  load**: inter-area mode **0.202 Hz / +15.5 %**, min damping **+6.8 %**, 0 unstable modes.
- ✅ **Wind as PMSM WECS** (`pmsm_pq_ss`) delivering exactly `p_s_ppc` (curtailed), init robust.
- ✅ **Inter-area study** (`run_link`) reproduces the growing/damped 0.2 Hz mode.
- ✅ **Overvoltage study** (`run_overvoltage`) reproduces reactive-surplus voltage rise and
  collapse via reactor tripping; AVR `E_min` saturation captured.
- ✅ **Generation-trip overvoltage cascade** (`run_gentrip_cascade`) — with the VRE **Q–V droop**
  (§12b) the solar act as reactive absorbers, so overvoltage-protection trips remove absorption
  → **positive-feedback cascade → collapse** with mis-set (1.04 pu) protection, while the
  compliant (1.10 pu) setting rides through. Faithful to Finding #7.
- ⚠️ `main.py` `DATA = esp_preblackout_dev.hjson`; the **base `esp_preblackout.hjson` still has
  wind in `pvs`** (not migrated). Wind (`pmsm_pq_ss`) does **not** yet carry the Q–V droop —
  only solar (`pv_pq_ss`) does.

---

## 16. Future steps

1. ✅ **DONE — VRE Q–V droop** on `pv_pq_ss` (§12b) + ✅ **overvoltage-trip cascade**
   (`run_gentrip_cascade`, §13): gen tripping now removes reactive absorption → overvoltage
   → documented positive-feedback runaway, resolved as a cascade.
2. **Q–V droop on the WECS wind** (`pmsm_pq_ss`) as well, so wind also contributes to the
   reactive-absorption cascade (currently only solar carries the droop).
3. **Islanding scenario:** open the FR AC tie (loss of synchronism / 12:33:21.5) with the
   Spanish AGC taking the slack — the sudden loss of the export path + charging surplus.
4. Migrate the **base `esp_preblackout.hjson`** to the WECS wind (keep base and dev consistent).
5. Calibrate solar **overvoltage protection settings**, the droop gain `K_qv`, and inverter
   current limits to the report's requirements (1.1 pu indefinite, 1.2 pu transient, no
   zero-delay). Present defaults (K_qv=3, mis-set 1.04 pu) are illustrative, not calibrated.
6. Optionally port the trapezoidal + `ca.Linsol` solvers into pydae core (see `prompt_add_linsol.md`).

---

## 17. Reproduce from scratch — checklist

1. **Install PyPSA-Spain** (`pypsa-spain`, `make install-lock-macos-arm`; env `pypsa-eur`).
2. **Configure** for Spain (`config/config_ES.yaml`, `pypsa_spain:` block) with **NUTS3 (adm)**
   clustering and the desired features (electricity_demand, interconnections, esios capacities).
3. **Solve/prepare a snapshot** near the pre-blackout instant; build the **France-as-slack**
   single-snapshot power flow (`preblackout_fr_slack.py`); apply the **nuclear must-run** and
   **generation-mix / load** matching.
4. **Convert** to pydae-bps: run `pypsa_to_pydae.py` → `esp_preblackout.hjson` (GENROU/GENSAL,
   `pv_pq_ss` VRE, terminal buses + transformers, shunt compensation, `SLACK_FR` AGC slack,
   400 kV base, `results` block).
5. **Set up the RAG** (`mdbib_bps`) with the incident report to ground calibration.
6. **Refine dynamics:** differentiate machine inertia/reactances (§11.1); tune PSS on the 10
   largest units to K_stab=10 (§11.2); clean the `results` reactive via `export_pf.py` (§11.3).
7. **Migrate wind to `pmsm_pq_ss`** (`wecs` section, `nu_w_0=16`, aero-curve pitch seed in the
   model's `xy_0_dict`) (§12).
8. **Init & analyze** (env `pydae_dev`):
   ```bash
   conda run -n pydae_dev python -c "import main; main.ini()"          # SSA
   conda run -n pydae_dev python -c "import main; main.run_link()"     # inter-area
   conda run -n pydae_dev python -c "import main; main.run_overvoltage()"  # overvoltage
   conda run -n pydae_dev python export_pf.py                          # solved PF + tables
   conda run -n base       python plot_sld_interactive.py             # interactive map
   ```

---

*RMS phasor model — captures electromechanical (0.2 Hz inter-area), control (AVR/gov/PSS),
and reactive-balance voltage dynamics. It does not represent sub-cycle electromagnetic
transients, harmonics, or unbalanced faults.*
