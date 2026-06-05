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
  python main.py
```

Generates `hydro_smib_step.png` (4-panel: power, frequency, gate/flow,
terminal voltage) showing the response to a load-reference step
$p_{c,lc} = 0.80 \to 0.95 \to 0.80$ pu.

## Expected behaviour

- **Non-minimum-phase response.** When the gate opens, the water
  column has inertia and the immediate effect is a *reduction* in
  turbine power before the flow accelerates and power climbs to the
  new setpoint. This is the canonical "inverse response" of a
  Francis-style hydro turbine — `p_m` and `p_g` dip below their
  pre-step value for the first few seconds of the gate move.
- **Slow recovery.** The transient-droop dashpot ($R_r$, $T_r$)
  rolls off over ~5 s; the LC integrator ($K_i = 0.01$, $\tau \approx
  100$ s) eventually pins $p_g$ to the setpoint as the dispatch
  catches up.
- **Tight voltage regulation.** SEXS with $K_a = 100$ keeps $V_1$
  near 1.0 pu through both steps.
