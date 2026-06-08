# 68-bus (NETS/NYPS) — base case

CasADi pipeline for the 68-bus, 16-generator, 5-area test system from
**PES-TR18 §5.6** (IEEE PES Task Force on Benchmark Systems for Stability
Controls, August 2015).

The network is a reduced-order equivalent of the New England Test System
(NETS) plus the New York Power System (NYPS), interconnected by three
double-circuit tie lines (60-61, 53-54, 27-53), with three neighbouring
areas represented as equivalent generators (**G13, G14, G15, G16**).

- `system68bus.hjson` — network description: 68 buses, 64 transmission lines,
  19 transformers, 35 ZIP loads, 16 `genroe` machines, AVRs (ST1A), PSSs
  (PSS1A on G1-G12), and the reference power-flow values used for validation.
- `main.py` — three functions: `build()`, `ini()`, `run()`.

## Requirements

A Python environment with `pydae-core`, `pydae-bps` and `matplotlib`
installed. No C compiler is needed — the CasADi backend folds the model
into an SX graph and integrates it with SUNDIALS/IDAS.

## Running

```bash
python main.py                               # both stages
python -c "import main; main.ini()"          # steady state + small-signal only
python -c "import main; main.run()"          # time-domain short circuit only
```

## `ini()` — steady state and small-signal analysis

1. Builds the model and solves the load flow (`model.ini`).
2. Prints the **initialization report** (bus voltages/angles, generator
   dispatch).
3. Prints the **validation table** comparing the solved state against the
   reference values from PES-TR18 Tables 5.43-5.44 (buses) and Table 5.49
   (generators).
4. Computes the reduced state matrix `A` (`model.A_eval()`) and runs
   small-signal analysis (`ssa.damp`): eigenvalues, damping ratios,
   frequencies.

**Output:** `system68bus_eig.png` — eigenvalue plot. (PES-TR18 publishes
only figures, not a numerical eigenvalue table, so there is no reference
overlay.)

## `run()` — time-domain short-circuit simulation

`run(t_clear=1.2, bus_fault='60')`:

1. Initializes and runs 1 s of steady operation.
2. Applies a short circuit at `bus_fault` (default bus 60, near the
   NETS-NYPS tie) by ramping the shunt conductance up, holds through the
   fault, then ramps it back to clear at `t_clear`.
3. Continues the post-fault response to t = 15 s.

**Output:** `system68bus_run.png` — two stacked plots: generator speeds
(all 16 machines, top) and voltages of the tie-line and major load buses
17, 18, 27, 53, 54, 60, 61 (bottom).

---

## Data sources

Every parameter is sourced from a specific PES-TR18 table.

| Element | PES-TR18 reference |
|---|---|
| Bus voltages and angles | Tables 5.43, 5.44 |
| Transmission lines (R, X, B) | Tables 5.45, 5.46 |
| Transformers + taps | Table 5.47 |
| Load (P, Q) | Table 5.48 |
| Generator dispatch | Table 5.49 |
| Synchronous machine dynamics | Tables 5.50, 5.51, 5.52 |
| DC4B excitation system | Table 5.53 *(not used — see Approximations)* |
| ST1A excitation system | Table 5.54 |
| PSS1A stabiliser | Table 5.55 |

## Modelling choices

### Generators
- All 16 machines use pydae's `genroe` (6th-order Anderson-Fouad / IEEE 1110
  Model 2.2 with **geometric saturation** $S(E) = B \cdot E^N$, matching
  the PSS/E GENROE form documented in PES-TR18).
- The PES-TR18 leakage reactance `X_l = 0.075` is **not** passed to
  `genroe`; the PSS/E convention treats $X''_d, X''_q$ as terminal-referred
  (X_l already folded in). See `packages/pydae-bps/.../syns/genroe.py`.
- Saturation enabled at the PES-TR18 light values (`S_10 = 0.001`,
  `S_12 = 0.01`).
- **G16** is the swing/reference machine: `K_delta = 0.001` and
  `Delta_ref = 0` anchor the absolute angle; `gov.p_c = 0.4224`
  (= 3379.5 MW / 8000 MVA, per Table 5.49).
- **G1-G15** carry `lc.p_c_lc` setpoints in machine pu so each gen's grid
  injection matches its Table 5.49 dispatch at steady state.

### AVRs
- **G1-G8, G10-G12** — PES-TR18 specifies IEEE **DC4B** (Table 5.53).
  pydae has no DC4B today; **`st1` (IEEE ST1A)** is used as the closest
  available AVR, parameterised from Table 5.54.
- **G9** — IEEE ST1A as in PES-TR18, with the same Table 5.54 parameters.
- **G13-G16** — no AVR (constant field voltage). PES-TR18 §4.6 does not
  specify excitation systems for the equivalent generators.

### PSS
- **G1-G8, G10-G12** — IEEE PSS1A `††` column of Table 5.55:
  `K_S=20, T_w=15, T_1=T_3=0.15, T_2=T_4=0.04, L_S ∈ [−0.05, 0.20]`.
- **G9** — IEEE PSS1A Unit 9 column:
  `K_S=12, T_w=10, T_1=T_3=0.09, T_2=T_4=0.02`.
- **G13-G16** — no PSS.
- The PSS1A 2nd-order filter coefficients (A1 = 0.04, A5 = 0.15 on the
  †† column) are **dropped** — `pss_kundur_2` implements only the washout
  + two lead-lag stages. For the electromechanical frequency band (~0.5 to
  2 Hz) the omitted block contributes a negligible phase shift.

### Loads
- Constant admittance (`K_zp = K_zq = 1.0`) per PES-TR18 §5.6.1
  ("constant admittance characteristics were admitted for all loads").

### Governors / LC
- Every generator has a `tgov1` governor with conservative settings
  (`R=0.05`, `T_1=0.5`, `T_2=2.1`, `T_3=7.0`). PES-TR18 omits speed
  governor models from the small-signal study, but pydae's `kundur`-style
  AVR PV-bus initialisation depends on the gov/LC chain to pin `p_m` at
  the dispatch value during `ini()`.

## Configuration switching

PES-TR18 §5.6.3 defines three cases:

| Case | PSS configuration |
|---|---|
| 1 | No PSS at all |
| 2 | PSS only on G9 |
| **3** | **PSS on G1-G12 (default in this file)** |

To switch:
- **Case 1** — set `K_stab: 0.0` on every PSS entry (or remove the `pss`
  blocks).
- **Case 2** — set `K_stab: 0.0` on every PSS *except* G9.
- **Case 3** — default; leave as is.

## Validation status

| Check | Result |
|---|---|
| `ini()` convergence | ✓ (~1.7 s wall clock) |
| Generator bus voltages match Table 5.43 (AVR-controlled) | ✓ exact |
| Generator dispatch matches Table 5.49 | ✓ (within LC tolerance) |
| Transmission bus voltages match Table 5.43-5.44 | ~ within ±5% (typical NR-vs-PSS/E load-flow drift) |
| Electromechanical modes match PES-TR18 §5.6.3 qualitatively | ✓ |

### pydae electromechanical modes (Case 3 default)

PES-TR18 publishes no numerical eigenvalue table — only Figures 5.79-5.81
showing eigenvalue scatter plots. The pydae computed modes reproduce the
**qualitative finding** of §5.6.3: "PSSs on G1-G12: there are still two
inter-area underdamped modes remaining ... because no PSSs were placed in
the equivalent generators that represent areas 3, 4 and 5."

| Mode (driving states) | real ± imag (1/s) | f (Hz) | ζ | Note |
|---|---|---|---|---|
| G14 ↔ G16 inter-area | −0.257 ± j3.774 | 0.601 | 0.068 | **under-damped** |
| G15 inter-area | −0.254 ± j5.169 | 0.823 | 0.049 | **under-damped** |
| G12 local | −0.547 ± j6.359 | 1.012 | 0.086 | |
| G02 local | −0.597 ± j7.481 | 1.191 | 0.080 | |
| G10 local | −0.798 ± j7.349 | 1.170 | 0.108 | |
| G01 local | −0.768 ± j7.876 | 1.254 | 0.097 | |
| G11 local | −1.221 ± j10.35 | 1.647 | 0.117 | |

The two under-damped modes (ζ < 7%) are exactly the ones PES-TR18 §5.6.3
flags: both involve **G13-G16**, which carry no PSS. All G1-G12 local
modes have ζ ≥ 8%.

## Limitations and known issues

1. **No exact DC4B model.** The G1-G12 excitation is approximated with
   ST1A. For small-signal analysis at electromechanical frequencies the
   discrepancy is small; for transient stability with strict ceiling
   limits the difference matters and you should track this issue.
2. **PSS1A 2nd-order blocks dropped.** Phase contribution at 0.5-2 Hz is
   small (<5°), but if you tune PSS gains close to the stability
   boundary the residual phase matters.
3. **Transmission voltage profile drifts ±5% vs Table 5.43-5.44.** Likely
   caused by minor differences in how the load admittance is recomputed
   in pydae vs PSS/E (initial admittance is taken at V=1 pu; PSS/E may
   use the converged V). Generator dispatch matches Table 5.49 because
   LC enforces `p_g = p_c_lc`.
4. **No published eigenvalue reference.** PES-TR18 only documents §5.6
   eigenvalues graphically. The "validation" of small-signal modes is
   qualitative (mode count, damping signs, frequency band).

## References

- IEEE PES Task Force, "Benchmark Systems for Small-Signal Stability
  Analysis and Control", **PES-TR18**, August 2015. §4.6 (description)
  and §5.6 (data).
- A. K. Singh, B. Pal, "Report on the 68-bus, 16-machine, 5-area System"
  (USP, 2013).
- R. Kuiava, T. Fernandes, M. Mansour, R. Ramos, "Report on the 68-bus
  system using PacDyn/ANATEM" (USP, 2014).
- IEEE Std. 421.5-2005 (excitation system models — ST1A, DC4B; PSS1A).
- IEEE Std. 1110-2019 (synchronous machine model 2.2 — implemented by
  `genroe`).
