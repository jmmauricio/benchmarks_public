"""Track synchronous-machine reactive power through the gen-trip overvoltage cascade.

Softened settings (K_qv=1, gentle seed, definite-time protection delay) so the
cascade develops over seconds and the reactive-balance trace is clean.
"""
import numpy as np
from matplotlib import pyplot as plt

import main
from pydae.core.model.casadi_model import CasadiModel

data = main.data
# arresting partial cascade -> clean reactive-balance trace that settles
K_qv, n_seed, V_trip, t_delay = 1.0, 4, 1.05, 0.6
t_settle, t_end, Dt, dt_check = 1.0, 12.0, 0.02, 0.1

syns = list(data['syns'])
S_n = {s['name']: s['S_n'] for s in syns}

model = CasadiModel(main.build())
model.Dt = Dt
model.ini({}, xy_0='xy_0.json')
solar = main.enable_qv_droop(model, K_qv=K_qv)
solar_names = [p['name'] for p in solar]
solar_bus = {p['name']: p['bus'] for p in solar}

reactors = []
for b in data['buses']:
    try:
        val = float(model.get_value(f"b_shunt_{b['name']}"))
    except Exception:
        continue
    if val < -1e-9:
        reactors.append((b['name'], val))
reactors.sort(key=lambda x: x[1])
seed = reactors[:n_seed]

model.run(t_settle, method='trapezoidal')
for bus, _ in seed:
    model.set_value(f"b_shunt_{bus}", 0.0)
model.recalculate_algebraics(tol=1e-8)
print(f"seed {len(seed)} reactors at t={t_settle}s, K_qv={K_qv}, "
      f"OV protection {V_trip} pu / {t_delay}s delay")

tripped, above, t = set(), {nm: 0.0 for nm in solar_names}, t_settle
while t < t_end - 1e-9:
    t_next = min(t + dt_check, t_end)
    try:
        model.run(t_next, method='trapezoidal')
    except RuntimeError as e:
        print(f"collapse at t~{t_next:.2f}s -- {e}")
        break
    t = t_next
    newly = []
    for nm in solar_names:
        if nm in tripped:
            continue
        V = float(model.get_value(f"V_{solar_bus[nm]}"))
        if V > V_trip:
            above[nm] += dt_check
            if above[nm] >= t_delay - 1e-9:
                newly.append(nm)
        else:
            above[nm] = 0.0
    for nm in newly:
        model.set_value(f"p_s_ppc_{nm}", 0.0)
        model.set_value(f"K_qv_{nm}", 0.0)
        tripped.add(nm)
    if newly:
        model.recalculate_algebraics(tol=1e-8)
        print(f"  t={t:5.2f}s  tripped {len(newly)} solar -> {len(tripped)}/{len(solar)} off")

model.post()
tt = np.array(model.Time)[:-1]        # drop last (non-converged) step

Qsyn = {s['name']: np.array(model.get_values(f"q_g_{s['name']}"))[:len(tt)] * S_n[s['name']] / 1e6
        for s in syns}
Qtot_es = np.sum([Qsyn[s['name']] for s in syns if s['name'] != 'SLACK_FR'], axis=0)
Qslack = Qsyn['SLACK_FR']
Qsolar = np.zeros(len(tt))
for p in solar:
    Qsolar += np.array(model.get_values(f"q_s_{p['name']}"))[:len(tt)] * p['S_n'] / 1e6

for label, idx in [('operating point', 0), ('final (pre-collapse)', -1)]:
    print(f"{label:22s} t={tt[idx]:5.2f}s  Q_syn(ES)={Qtot_es[idx]:+6.0f}  "
          f"Q_slack(FR)={Qslack[idx]:+6.0f}  Q_solar={Qsolar[idx]:+6.0f} Mvar")
dQ = {nm: Qsyn[nm][-1] - Qsyn[nm][0] for nm in Qsyn if nm != 'SLACK_FR'}
print("biggest ES machine swings [Mvar]:")
for nm in sorted(dQ, key=lambda n: dQ[n])[:6]:
    print(f"  {nm:20s} {Qsyn[nm][0]:+7.0f} -> {Qsyn[nm][-1]:+7.0f}  (dQ={dQ[nm]:+.0f})")

fig, ax = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
big = sorted([s for s in syns if s['name'] != 'SLACK_FR'], key=lambda s: -s['S_n'])[:6]
for s in big:
    ax[0].plot(tt, Qsyn[s['name']], lw=1.3, label=s['name'])
ax[0].set_ylabel('Machine Q [Mvar]'); ax[0].grid(True); ax[0].legend(loc='best', fontsize=7)
ax[0].set_title(f'Synchronous-machine reactive power during gen-trip cascade '
                f'(K_qv={K_qv}, OV {V_trip} pu / {t_delay}s delay)')
ax[1].plot(tt, Qtot_es, 'b-', lw=1.8, label='Σ Spain syn Q')
ax[1].plot(tt, Qslack, 'navy', lw=1.2, ls='--', label='France slack Q')
ax[1].plot(tt, Qsolar, 'r-', lw=1.6, label='Σ solar Q (droop)')
ax[1].axhline(0, color='k', lw=0.6)
ax[1].set_ylabel('Total Q [Mvar]'); ax[1].grid(True); ax[1].legend(loc='best', fontsize=8)
for s in big:
    ax[2].plot(tt, np.array(model.get_values(f"v_f_{s['name']}"))[:len(tt)], lw=1.2, label=s['name'])
ax[2].axhline(0.0, color='k', ls=':', lw=1, label='E_min=0')
ax[2].set_ylabel('AVR field v_f [pu]'); ax[2].set_xlabel('Time [s]'); ax[2].grid(True)
ax[2].legend(loc='best', fontsize=7)
fig.tight_layout(); fig.savefig('esp_preblackout_gentrip_synq.png')
print('saved esp_preblackout_gentrip_synq.png')
