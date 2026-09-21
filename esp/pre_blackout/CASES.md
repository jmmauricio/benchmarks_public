# The three cases, and how to obtain each

This directory carries one network and three distinct dynamic cases built from it.
They differ in **dispatch**, **interconnection topology** and **stabiliser tuning**,
and nothing else. This document says exactly how each is produced, what it is for,
what it produces, and what may be claimed about it.

| | case | file | for |
|---|---|---|---|
| 1 | Pre-blackout, 28 April 2025 ≈ 12:20 | `esp_case1_preblackout.hjson` | reproducing the observed 0.2 Hz inter-area oscillation |
| 2 | Stressed N-1, at the 5% criterion | `esp_case2_pod_bench.hjson` | developing POD controllers on PV plant |
| 3 | Well-damped normal operation | `esp_case3_normal.hjson` | a healthy reference / control case |

Everything derives from `esp_preblackout_dev.hjson`, the base network.

**The three case files are self-contained.** Each carries its own dispatch,
topology and stabiliser settings, and reproduces the modes stated in its own header
comment with nothing passed in — no arguments, no environment variables. Build one
and initialise it and you get the documented numbers.

```bash
cd esp/pre_blackout
python tools/make_cases.py      # regenerate all three from the base network
python tools/verify_cases.py    # check each reproduces its header
```

The sections below explain how each is constructed and why, and the intermediate
files (`esp_preblackout_ce_chain.hjson`, `_pod_bench`, `_pss_tuned`) are the working
artefacts those steps produce.

---

## Before anything: two traps that will cost you a day

**`model.ini()` is not reproducible on a model instance that has already been
initialised.** A second call on the same object converges to a *different*
operating point — `‖x‖` was observed jumping from 455 to 5004, with an entirely
different spectrum. Every script here therefore builds a fresh model per variant.
If you need an open-loop baseline to compare against, take it from
`TuningReport.floor_before` / `.modes_before`, never from a second `ini()`.

Two exceptions, both used deliberately: **stabiliser parameters and governor droop
do not move the equilibrium** (`v_pss` is zero at steady state, and droop acts on a
speed deviation that is zero there). So gains and time constants may be changed and
`eig()` re-run without re-initialising. Network reactances and dispatch may not.

**Bus `FR` is declared `U_kV: 1`.** `pydae/bps/lines/lines.py` converts a line's
ohms to per-unit using **`bus_j`'s** base voltage. A line written as
`bus_j: "FR"` therefore gets `Z_base = 0.01 Ω`, turning 22 Ω into 2200 pu — a
silent open circuit that builds, solves, and yields plausible-looking modes. Always
put the 400 kV bus first: `bus_j: "CE_C", bus_k: "FR"`.

---

## Case 1 — Pre-blackout, 28 April 2025 around 12:20

### What it is for

Reproducing the **East-Centre-West inter-area oscillation at 0.2 Hz** that ran from
12:19 to 12:22, at up to 200 mHz frequency amplitude, and subsided only after
operator countertrading. Use it to study the phenomenon. Do **not** use it as a POD
test bench — it contains an unstable mode (see below).

### Why the base network cannot produce it

The 0.2 Hz mode is a mode of the **whole Continental Europe synchronous area**
(incident report §4.2.3.1), slow because it spans CE, with Iberia merely at its West
end. The base case represents everything past the border as one stiff machine
(`SLACK_FR`, 50 GVA, H = 40, i.e. 2,000,000 MW·s), which yields an *Iberia versus
France* mode at ~0.33 Hz and nothing slower.

Two measurements settle the point:

- **The tie flow moves the wrong way.** Read off Figures 2-22 and 2-23 (RE SCADA),
  the ES–FR physical flow ran **1,300–1,450 MW through 12:15–12:22**, peaking ~1,450
  at 12:21.5. The 866 MW (12:24) and 190 MW (12:28) quoted in the report text are the
  countertrading *response*, after the event. Correcting the model's 2,025 MW down to
  ~1,450 MW moves the mode from 0.327 Hz / 9.28% to **0.343 Hz / 10.64%** — higher
  frequency and *more* damping.
- **The inertia is already right.** Model Iberian stored energy is 106,754 MW·s
  against the 119,474 MW·s reported at 12:30 (Table 2-4) — 11% low. Inertia is not
  the missing ingredient, and raising it lowers frequency only slowly.

So the missing element is structural: Continental Europe has to be able to oscillate.

### How it is built

`SLACK_FR` is replaced by a three-group chain, Iberia at the West end:

```
Iberia ══(existing ES-FR ties)══ West ══X_wc══ Centre ══X_ce══ East
```

`SLACK_FR`'s 2,000,000 MW·s is *redistributed*, not changed: West 600,000, Centre
800,000, East 600,000 MW·s (H = 40 throughout, sized via `S_n`). Total CE inertia is
identical; only its distribution differs, which is the entire point.

The three equivalents also receive **area-equivalent machine data** — rotor flux time
constants long enough that the fluxes are effectively constant over the mode period,
i.e. the classical machine form that `gencls` exists for. This is not cosmetic:
`SLACK_FR`'s `T1q0 = 0.8 s` puts the q-axis corner at 1.25 rad/s, directly on top of a
0.2 Hz mode, and that one parameter alone contributed **12 points of damping**
(17.27% → 5.07% when frozen). An aggregate representing a whole synchronous area
must not inherit one turbo-generator's q-axis damping.

### Commands

```bash
cd esp/pre_blackout

# 22 ohm on both CE corridors, from the ORIGINAL stabiliser settings
CHAIN_SRC=$PWD/esp_preblackout_dev.hjson \
python tools/make_chain.py 22 22

# modes, with the four-group shape; the 625 argument trims Iberian generation
# to bring the ES-FR export to the observed ~1,450 MW
python tools/chain_modes.py 625
```

`make_chain.py` takes `X_wc X_ce [KE_W KE_C KE_E]`, all in ohms and MW·s.

### What it produces

```
ES-FR export = 1,451 MW

  f [Hz]  zeta%       Iberia          West         Centre           East
  0.202    3.53   1.000@175.8   0.804@179.9    0.083@ -8.0    0.872@  0.0
  0.301    1.54   0.994@  0.7   0.414@ -0.9    1.000@179.8    0.743@  0.0
  0.384   -2.94   1.000@162.5   0.261@ -8.4    0.074@-171.1   0.009@  0.0
```

The **0.202 Hz mode is the East-Centre-West mode**: Iberia and West swing together in
antiphase with East, Centre near a node (amplitude 0.083). An independent four-mass
analytic chain predicts 0.183 Hz with the same shape, which is a reasonable check.

**Damping is 3.53%, and is reported as obtained — not tuned to a target.** For
reference the incident report's own sensitivity study finds ~3% even with Spanish
inertia tripled, so this is the right order.

### Read this before using it

- **The 0.384 Hz mode is unstable at −2.94%.** This case is *not* a stable system.
  That is a consequence of running the original stabiliser settings (43 of 55 at zero
  gain, 7 more at a token 0.01) on this network. It is fine for studying the 0.2 Hz
  oscillation and wrong for anything needing a stable base point.
- **The export trim is a runtime argument, not stored in the file.** The hjson gives
  2,025 MW; the 1,451 MW figure comes from the `625` passed to `chain_modes.py`. Any
  other script must apply the same trim or it is analysing a different operating
  point.
- The CE group sizes and corridor reactances are **plausible, not sourced**. They
  are chosen so the West-end mode lands at 0.2 Hz. Say so if you publish this.

### Tuning the chain

Frequency responds to the corridor reactance roughly as `f ∝ X^-0.42`:

| X_wc = X_ce | f | ζ |
|---|---|---|
| 27 Ω | 0.184 Hz | 6.61% |
| 22 Ω | **0.202 Hz** | 3.53% |
| 15 Ω | 0.192 Hz* | 16.20%* |

\* measured before the equivalents were given area-equivalent data; the frequency
relation holds, the damping does not.

---

## Case 2 — Stressed N-1 at the 5% criterion (the POD bench)

### What it is for

Developing and evaluating **power oscillation damping controllers on PV plant**. It
is stable, nothing is unstable anywhere, and the inter-area mode is lightly damped so
a controller has something to act on and something to improve.

Full framing, provenance and suggested paper wording: **`README_pod_bench.md`**.
Read it before publishing anything based on this case.

### How it is built

Three ingredients, none of which is a statement about any operator's practice:

1. **Forward-looking converter penetration.** Synchronous plant committed at 52% of
   base by stored energy, PV taking up the released 4,072 MW. The machines here are
   plant aggregates, so scaling `S_n` and `p_c_lc_mw` together is exactly *fewer
   units of that plant committed* — per-unit machine data untouched.
2. **N-1 on one ES–FR circuit.** One of the two circuits is commented out. The
   exchange is reduced to ~1,000 MW accordingly: the remaining circuit has
   `P_max = 400²/60.4 = 2,647 MW`, so the pre-contingency 1,950 MW would sit at 74%
   of it and ~48° across the tie — past where the load flow will solve at all. The
   trim is part of the contingency scenario, not a numerical convenience.
3. **Stabilisers tuned to the NTS-SEPE criterion of 5%**, not beyond it, by the
   residue method (`pydae.bps.utils.pss_tuner.tune_psss`).

### Commands

```bash
cd esp/pre_blackout

BENCH_N1=1 BENCH_EXPORT_TRIM=1000 BENCH_ZETA=0.05 \
python tools/finalize_bench.py 0.52
```

That regenerates the dispatch, applies N-1, tunes the stabilisers **and writes them
into the hjson**, then prints the modes, the sensitivity row and the residues.

To explore without writing the file, use `tools/bench_eval.py` with the same
environment variables — it tunes in memory only.

| knob | meaning |
|---|---|
| positional arg | committed synchronous capacity, as a fraction of base |
| `BENCH_N1=1` | open one ES–FR circuit |
| `BENCH_EXPORT_TRIM` | MW of PV backed off, to reduce the exchange |
| `BENCH_ZETA` | stabiliser tuning target |
| `BENCH_BASELOAD=1` | hold nuclear and CCGT at full commitment (see caveat) |

### What it produces

| | |
|---|---|
| Iberian stored energy | 55,512 MW·s (46% of the 119,474 MW·s at 12:30) |
| Synchronous share | 15.9% — 4,411 MW vs 19,970 MW PV + 3,398 MW wind |
| ES–FR | one circuit, 988 MW export |
| **Inter-area mode** | **0.262 Hz, ζ = 5.14%** |
| Worst mode in the system | 1.113 Hz, ζ = 4.24% |
| Committed stabilisers | `hydro_ES432_hydro` (K = 14.69), `hydro_ES113_hydro` (K = 0.38) |

**Sensitivity — the table the paper rests on:**

| stabiliser state | worst mode | inter-area mode |
|---|---|---|
| As shipped, tuned to 5% | 4.24% | 0.262 Hz, **5.14%** |
| Removed | **−16.22%** | 0.296 Hz, −16.22% |

The margin is *not* obtained by weakening stabilisers. Without them the same case is
unstable by 16 points.

### Read this before using it

- **The obvious stronger claim is false.** Synchronous machines *can* damp this mode:
  re-tuning against a 10% target reaches ~14%, at several commitment levels, with and
  without N-1. The supportable claim is narrower — tuned to the *criterion* rather
  than beyond it, the system sits at the criterion, and the inter-area mode is what
  sits on it. Damping here is set by the tuning target, by construction.
- **5% is the floor compatible with this framing.** Reaching 1–3% requires either
  tuning below the grid code, or commitment below ~50% where the system collapses to
  −15% and worse with no usable window between. Both reintroduce the problem the
  framing exists to avoid.
- **`BENCH_BASELOAD=1` cannot reach the marginal regime.** Nuclear and the large
  CCGTs carry 67,082 of the 106,754 MW·s in this case — 63% of the inertia — so
  holding them floors stored energy at 63% of base whatever else is decommitted,
  close to the uniform-60% case that came out comfortably damped. Merit order alone
  is not enough; getting there means taking some baseload off too.

---

## Case 3 — Well-damped normal operation

### What it is for

A healthy reference. Use it as the control case against which cases 1 and 2 are
compared, and as the "before" in any study of what stress does.

### How it is built

Base network, base dispatch, no contingency — the only change from
`esp_preblackout_dev.hjson` is that the stabilisers are tuned by the residue method
against a **10% target**, which is beyond the 5% criterion.

The tuner keeps 3 of 55 stabilisers and comments the rest out rather than leaving
them at zero gain, so the model carries no inert stabiliser states — 976 states
becomes 820.

### Commands

```bash
cd esp/pre_blackout
python - <<'PY'
from pydae.bps import BpsBuilder
from pydae.core.builder.casadi_builder import CasadiBuilder
from pydae.core.model.casadi_model import CasadiModel
from pydae.bps.utils.pss_tuner import tune_psss

grid = BpsBuilder('esp_preblackout_dev.hjson', use_casadi=True)
grid.uz_jacs = False
grid.construct('esp_case3')
model = CasadiModel(CasadiBuilder(grid.sys_dict).build())

report = tune_psss(model, 'esp_preblackout_dev.hjson', zeta_target=0.10,
                   xy_0='xy_0.json', write='esp_preblackout_pss_tuned.hjson')
print(report.markdown())
PY
```

### What it produces

| machine | K_stab | T₁ = T₃ | T₂ = T₄ |
|---|---|---|---|
| `CCGT_ES620` | 6.136 | 0.1534 | 0.0100 |
| `hydro_ES113_hydro` | 12.218 | 0.1521 | 0.0100 |
| `hydro_ES432_hydro` | 2.125 | 0.1119 | 0.0100 |

Worst electromechanical damping **−3.15% → +6.33%**; every mode stable.

| mode | f | ζ |
|---|---|---|
| worst | 1.242 Hz | 6.33% |
| | 0.996 Hz | 6.50% |
| inter-area | 0.327 Hz | 9.28% |

### Read this before using it

- **It stops at 6.33%, short of the 10% asked for**, and the tuner's skip table says
  why: machine after machine needs 145–178° of lead at the 1.0–1.24 Hz modes, i.e.
  73–89° per stage against a 70° limit. That is the excitation system's phase
  roll-off, not a tuner limitation — at 1.74 Hz *every* machine in this case,
  including the 3.8 GVA CCGT, has `arg(GEP) ≈ −167°`. Three lead-lag stages or faster
  excitation is the lever, not more tuning.
- The three surviving designs carry 24–36° of residual phase error and all pin `T₂`
  at its 10 ms floor. Two stages cannot properly cover 0.34 Hz and 1 Hz at once. The
  designs work, but they are stretched.

---

## Cases 4a / 4b / 4c — the robustness family

### What they are for

Case 2 gives one inter-area mode at one frequency. A POD controller is phase
compensated exactly like a stabiliser, so its design is frequency sensitive in
exactly the same way — and in this system that bites hard: the residue-tuned
networks already carry 24–36° of residual phase error and pin `T₂` at its floor,
because two lead-lag stages cannot cover a wide band at once. A controller tuned at
0.262 Hz has no demonstrated validity at 0.21 or 0.40 Hz.

These three hold the dispatch of case 2 fixed and vary **only** the interconnection,
so the controller can be designed against a locus rather than a point. This is the
X_L sweep NTS-SEPE prescribes for robustness assessment, which is what makes 4c's
raised reactance a sanctioned sensitivity rather than an invented operating point.

### Commands

```bash
cd esp/pre_blackout
python tools/make_robustness.py          # all of 4a, 4b, 4c and 5
python tools/make_robustness.py 4c 5     # or just some
```

### What they produce

| | interconnection | inter-area mode | worst mode | unstable |
|---|---|---|---|---|
| 4a `esp_case4a_tie_strong.hjson` | both circuits | **0.396 Hz, 6.25%** | 1.514 Hz, 6.15% | 0 |
| 4b `esp_case4b_tie_base.hjson` | one circuit (N-1) | **0.262 Hz, 5.14%** | 1.113 Hz, 4.24% | 0 |
| 4c `esp_case4c_tie_weak.hjson` | one circuit, 73 Ω | **0.212 Hz, 12.07%** | 1.085 Hz, 5.39% | 0 |

The modal locus spans **0.21 – 0.40 Hz**, which is the band a POD design has to
cover. 4b is identical to case 2 and is included so the family is self-contained.

Stabilisers are re-tuned to the 5% criterion on **each** variant. Carrying one
tuning across a changed network would make any margin an artefact of stale tuning
rather than a property of the case.

### Read this before using them

**Weakening the tie is not a monotone route to a lightly damped mode.** Damping
ratio is ζ = −σ/|λ|. Weakening the interconnection collapses the synchronising
torque, so ω falls sharply while σ barely moves, and the *ratio* inflates even
though the mode is no better damped in absolute terms. At 104 Ω the inter-area mode
came out at 0.122 Hz with ζ = 43%. 4c is set at 73 Ω, past which the family stops
being useful. If you extend the sweep, watch σ as well as ζ.

That is also why damping is not constant across the family (5.1% → 12.1%): the
sweep moves frequency and damping together, and cannot be used to isolate one.

## Case 5 — low irradiance

### What it is for

The first question anyone will ask about PV-based damping is what happens when the
sun is not shining. In case 2 the PV units run at about 0.45 pu, so there is
modulation headroom in both directions — a favourable condition, chosen. This case
documents where the approach loses authority.

Synchronous commitment is raised to 120% of base with PV backed off to cover it, so
a PV-based controller has much less to modulate.

```bash
python tools/make_robustness.py 5
```

| | |
|---|---|
| Inter-area mode | **0.261 Hz, ζ = 10.80%** |
| Worst mode | 1.022 Hz, ζ = 5.26% |
| Unstable | none |

The compensating point, and it is worth making explicitly: the system is also **less
stressed** in this condition. More synchronous plant is committed, so the inter-area
mode is better damped (10.8% against 5.1%) precisely when the PV controller has
least authority. Reduced authority is therefore not as costly as it first appears —
which is a stronger result for a paper than quietly omitting the case.

An earlier attempt at 150% commitment produced a 0.936 Hz mode at 0.42% damping —
nearly unstable, and the tuner could not recover it. 150% is an aggressive
extrapolation of machines that are already plant aggregates. 120% is the setting
that behaves.

## Summary

| | dispatch | interconnection | stabilisers | inter-area mode | anything unstable? |
|---|---|---|---|---|---|
| 1 Pre-blackout | base, export 1,451 MW | intact + CE chain | as-found (mostly off) | 0.202 Hz, 3.53% | **yes**, 0.384 Hz at −2.94% |
| 2 POD bench | 52% sync, export 988 MW | **N-1** | tuned to 5% | 0.262 Hz, 5.14% | no |
| 3 Normal | base, export 2,025 MW | intact | tuned to 10% | 0.327 Hz, 9.28% | no |
| 4a Tie strong | 52% sync, 1,036 MW | both circuits | tuned to 5% | 0.396 Hz, 6.25% | no |
| 4b Tie base | 52% sync, 988 MW | N-1 | tuned to 5% | 0.262 Hz, 5.14% | no |
| 4c Tie weak | 52% sync, 988 MW | N-1 at 73 Ω | tuned to 5% | 0.212 Hz, 12.07% | no |
| 5 Low irradiance | 120% sync, 367 MW | N-1 | tuned to 5% | 0.261 Hz, 10.80% | no |

Case 1 answers *what happened*. Case 2 is what you develop controllers against, and
4a/4b/4c are the locus you demonstrate robustness over. Case 3 is what a healthy
system looks like, and case 5 is where the approach runs out of resource.

Every one of these is checked against its own header by `tools/verify_cases.py`,
which exits non-zero on a mismatch.
