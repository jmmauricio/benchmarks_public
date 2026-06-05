# Hydro SMIB benchmark — gensal + HYGOV + SEXS

Single 250 MVA salient-pole hydro unit, step-up transformer to a 230 kV
infinite bus. Demonstrates the pydae `gensal` machine model
(IEEE 1110 Model 2.1 / PSS/E `GENSAL`) wired to a HYGOV turbine governor
and a SEXS AVR — the canonical hydro stack.

## Parameters (Kundur §13.6 style)

- 250 MVA, 13.8 kV, 50 Hz, $H = 4$ s
- Salient-pole: $X_d = 1.05$, $X_q = 0.65$ pu
- $X'_d = 0.30$, $T'_{d0} = 6.0$ s
- $X''_d = X''_q = 0.22$ pu (terminal-referred); $T''_{d0} = 0.038$ s,
  $T''_{q0} = 0.099$ s
- $S_{1.0} = 0.10$, $S_{1.2} = 0.40$ — d-axis saturation
- HYGOV: $R = 0.05$ (permanent droop), $R_r = 0.5$, $T_r = 5$ s
  (dashpot), $T_w = 1$ s (water column), $A_t = 1.2$,
  $Q_{nl} = 0.08$
- SEXS: $K_a = 100$, fast lead-lag

## Run

```bash
PYTHONPATH=/path/to/pydae/packages/pydae-uds/src \
  python main.py         # gate-opening step  (0.80 → 0.95 → 0.80)
PYTHONPATH=/path/to/pydae/packages/pydae-uds/src \
  python load_trip.py    # load trip          (0.80 → 0.20)
```

`main.py` generates `hydro_smib_step.png` — load-reference step
$p_{c,lc} = 0.80 \to 0.95 \to 0.80$ pu, exercising the
**gate-opening** inverse response.

`load_trip.py` generates `hydro_smib_load_trip.png` — a sudden
load-trip step $p_{c,lc} = 0.80 \to 0.20$ pu at $t = 5$ s, exercising
the **gate-closing** inverse response (output briefly rises while the
water column keeps flowing through the partly-closed gate).

## Expected behaviour

- **Gate-opening inverse response (`main.py`).** When the gate opens,
  the water column has inertia and the immediate effect is a
  *reduction* in turbine power before the flow accelerates and power
  climbs to the new setpoint. `p_g` dips below the pre-step value for
  the first few seconds of the gate move.
- **Gate-closing inverse response (`load_trip.py`).** Symmetrically,
  when the gate closes (a sudden drop in dispatch — load was shed or
  another unit picked up demand), the water that was already in motion
  has nowhere to go and momentarily *raises* the turbine output before
  the column decelerates and power falls to the new setpoint. The
  observed peak is ≈ 210 MW (0.84 pu, ~5 % above the pre-trip 200 MW)
  at $t \approx 17$ s with $T_w = 1$ s; with a longer penstock the
  overshoot would be larger.
- **Slow recovery in both cases.** The transient-droop dashpot
  ($R_r$, $T_r$) rolls off over ~5 s; the LC integrator
  ($K_i = 0.01$, $\tau \approx 100$ s) eventually pins $p_g$ to the
  setpoint as the dispatch catches up.
- **Tight voltage regulation.** SEXS with $K_a = 100$ keeps $V_1$
  near 1.0 pu through both scenarios.
