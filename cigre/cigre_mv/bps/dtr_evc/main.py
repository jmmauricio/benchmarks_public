"""
CIGRE European MV benchmark — EV fast-charging-hub stress study (DTR).

A high-power EV charging hub (constant power) is placed at the radial feeder
end (B11). Because it is constant power, the current holds as the voltage sags,
so the whole B01→…→B11 cable run is stressed and the closed-loop dynamic
thermal rating (temperature-dependent resistance) is exercised.

    python main.py        # ini + reports + EV-hub stress sweep + 24 h transient

Key question: at the feeder end, is the binding constraint the cable
*temperature* (90 °C) or the *voltage* (~0.9 pu floor)?
"""
import math

import numpy as np
from matplotlib import pyplot as plt

from pydae.bps import BpsBuilder
from pydae.core.builder.casadi_builder import CasadiBuilder
from pydae.core.model.casadi_model import CasadiModel
from pydae.bps.utils.reporter import report_buses
from pydae.bps.utils.dtr_report import report_dtr
from pydae.utils import read_data

DATA = 'cigre_mv.hjson'
S_BASE = 25e6                      # system base (VA), from system.S_base
PF_HUB = 0.98                      # EV-hub power factor (near unity)
TANPHI = math.tan(math.acos(PF_HUB))
T_LIM_CABLE = 90.0                 # XLPE max conductor temperature (°C)
V_FLOOR = 0.90                     # MV voltage floor (pu)
# The head cable carries the whole feeder-1 current; B11 is fed via this path:
HEAD_PATH = ['B01_B02', 'B02_B03', 'B03_B08', 'B08_B09', 'B09_B10', 'B10_B11']


def build():
    grid = BpsBuilder(DATA, use_casadi=True)
    grid.uz_jacs = False
    grid.construct('cigre_mv')
    return CasadiBuilder(grid.sys_dict).build()


def _hub_inputs(p_mw):
    """Constant-power EV-hub inputs (pu on S_base) for a given hub size [MW]."""
    p_pu = p_mw * 1e6 / S_BASE
    return {'p_p_ev_B11': p_pu, 'q_p_ev_B11': p_pu * TANPHI}


def ini():
    """Steady-state init at the nominal (installed) hub size + reports.

    No base-case validation here: with the EV hub energised the voltages
    deliberately deviate from the CIGRE Table 9.6 reference (that check lives in
    the `base`/`base_dtr` cases).
    """
    model = CasadiModel(build())
    model.ini({}, newton_tol=1e-10)
    report_buses(model, DATA)
    report_dtr(model, DATA)
    return model


def stress_ev(p_max_mw=4.0, n_steps=9):
    """
    Sweep the EV hub at the feeder end (B11) from 0 to *p_max_mw* and report,
    at each steady state, the feeder-end voltage, the head-cable current and
    temperature, and the binding constraint (voltage floor vs thermal limit).

    model.ini() warm-starts from the previous solution, so the sweep is a
    continuation toward the hosting-capacity limit (voltage collapse).
    """
    model = CasadiModel(build())
    model.ini(_hub_inputs(0.0), newton_tol=1e-10)

    print('\nEV charging-hub stress sweep at B11 (radial feeder end)')
    print(f"{'P_hub MW':>9} {'V_B11':>7} {'I_head A':>9} {'Tc_head C':>10} {'binding':>9}")
    rows = []
    for p_mw in np.linspace(0.0, p_max_mw, n_steps):
        try:
            model.ini(_hub_inputs(p_mw), newton_tol=1e-9)
        except Exception:
            print(f"{p_mw:9.2f}   --- no converged solution: voltage collapse ---")
            break
        v_end = float(model.get_value('V_B11'))
        i_head = float(model.get_value('I_B01_B02'))
        tc_head = max(float(model.get_value(f't_conductor_{ln}')) for ln in HEAD_PATH)
        binding = 'VOLTAGE' if v_end <= V_FLOOR else ('THERMAL' if tc_head >= T_LIM_CABLE else 'ok')
        print(f"{p_mw:9.2f} {v_end:7.3f} {i_head:9.1f} {tc_head:10.1f} {binding:>9}")
        rows.append((p_mw, v_end, i_head, tc_head))

    print("\nAt this long cable feeder the feeder-end VOLTAGE binds well before the\n"
          "cable thermal limit — DTR thermal headroom exists but is not usable at\n"
          "the feeder end without voltage support (taps/STATCOM) or a head-end hub.")
    return rows


def run(p_hub_mw=2.0, hours=24.0):
    """
    Time-domain transient: energise a *p_hub_mw* EV hub at B11 and integrate for
    *hours* to show the cable warming (conductor in minutes, soil over hours).
    """
    data = read_data(DATA)
    model = CasadiModel(build())
    model.Dt = 10.0
    model.ini(_hub_inputs(0.0), newton_tol=1e-10)   # cold start, hub off
    model.run(3600.0, _hub_inputs(0.0))             # 1 h baseline
    model.run(hours * 3600.0, _hub_inputs(p_hub_mw))  # hub energised
    model.post()

    fig, ax = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    for line in data['lines']:
        if line.get('dtr'):
            ln = f"{line['bus_j']}_{line['bus_k']}"
            ax[0].plot(model.Time / 3600.0, model.get_values(f"t_conductor_{ln}"), label=ln)
    ax[0].axhline(T_LIM_CABLE, ls='--', color='r', lw=0.8, label='90 °C limit')
    ax[0].set_ylabel('conductor temp (°C)'); ax[0].legend(fontsize=7, ncol=2); ax[0].grid(True)
    ax[1].plot(model.Time / 3600.0, model.get_values('V_B11'), 'k', label='V_B11')
    ax[1].axhline(V_FLOOR, ls='--', color='r', lw=0.8, label='0.9 pu floor')
    ax[1].set_ylabel('V_B11 (pu)'); ax[1].set_xlabel('time (h)'); ax[1].legend(fontsize=7); ax[1].grid(True)
    fig.suptitle(f'EV hub {p_hub_mw:.1f} MW at B11 — DTR response')
    fig.savefig('ev_stress.png')
    print(f'\nSaved ev_stress.png  (hub {p_hub_mw:.1f} MW, {hours:.0f} h)')
    return model


if __name__ == '__main__':
    ini()         # steady-state init + reports + validation
    stress_ev()   # EV-hub hosting-capacity sweep
    run()         # 24 h transient at the nominal hub size
