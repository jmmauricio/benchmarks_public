"""Hydro SMIB load-trip scenario — gensal + HYGOV + SEXS.

Drops the dispatch setpoint p_c_lc from 0.80 → 0.20 pu at t = 5 s,
simulating the moment when the load this unit was supplying is shed
(or another generator picks it up). The gate closes against a still-
flowing water column, producing the **closing inverse response**:
the rotor sees a brief power *rise* before output decays to the new
operating point.

This is the symmetric counterpart of the gate-opening response shown
in main.py: opening dips before rising, closing rises before falling.
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from pydae.bps import BpsBuilder
from pydae.core import Builder, Model


def main():
    here = Path(__file__).resolve().parent
    grid = BpsBuilder(str(here / "hydro_smib.hjson"))
    grid.checker()
    grid.uz_jacs = False
    grid.construct("hydro_smib")

    Builder(grid.sys_dict, target="cffi", sparse=False).build()
    m = Model("hydro_smib")
    m.Dt = 0.005
    m.decimation = 4

    ok = m.ini({}, xy_0=dict(grid.dae["xy_0_dict"]))
    assert ok, "ini failed"
    print(f"initial: p_g={m.get_value('p_g_1'):.4f} pu, "
          f"x_g={m.get_value('x_g_gov_1'):.4f}, "
          f"V_1={m.get_value('V_1'):.4f}")

    m.run( 5.0, {})                          # settle 5 s at 0.8 pu
    m.run(60.0, {"p_c_lc_1": 0.20})          # LOAD TRIP — drop to 0.2 pu
    m.run(60.0, {})                          # let it settle
    m.post()

    t = m.Time
    S_n_MW = 250.0

    fig, axes = plt.subplots(4, 1, figsize=(9, 9), sharex=True)
    axes[0].plot(t, m.get_values("p_g_1") * S_n_MW, label="p_g", color="C0")
    axes[0].plot(t, m.get_values("p_m_1") * S_n_MW, label="p_m",
                 color="C1", ls="--")
    axes[0].plot(t, m.get_values("p_c_lc_1") * S_n_MW, label="p_c_lc",
                 color="k", ls=":")
    axes[0].axvline(5.0, color="red", lw=0.7, alpha=0.5)
    axes[0].set_ylabel("Power (MW)")
    axes[0].grid(); axes[0].legend()

    axes[1].plot(t, (m.get_values("omega_1") - 1.0) * 100.0, color="C2")
    axes[1].axvline(5.0, color="red", lw=0.7, alpha=0.5)
    axes[1].set_ylabel(r"$(\omega-1)\times 100$ (%)")
    axes[1].grid()

    axes[2].plot(t, m.get_values("x_g_gov_1"), label=r"gate $x_g$", color="C3")
    axes[2].plot(t, m.get_values("q_gov_1"),   label=r"flow $q$",
                 color="C4", ls="--")
    axes[2].axvline(5.0, color="red", lw=0.7, alpha=0.5)
    axes[2].set_ylabel("Gate / flow")
    axes[2].grid(); axes[2].legend()

    axes[3].plot(t, m.get_values("V_1"), color="C5")
    axes[3].axvline(5.0, color="red", lw=0.7, alpha=0.5)
    axes[3].set_ylabel(r"$V_1$ (pu)")
    axes[3].set_xlabel("Time (s)")
    axes[3].grid()

    fig.suptitle("Hydro SMIB — load trip at t=5 s, p_c_lc: 0.80 → 0.20 pu")
    fig.tight_layout()
    out = here / "hydro_smib_load_trip.png"
    fig.savefig(out, dpi=130)
    print(f"saved {out}")

    # Report the closing inverse response peak
    i_trip = int(np.searchsorted(t, 5.05))
    win = m.get_values("p_g_1")[i_trip:i_trip + int(20.0 / (m.Dt * m.decimation))]
    i_peak_rel = int(np.argmax(win))
    t_peak = t[i_trip + i_peak_rel]
    p_peak = win[i_peak_rel] * S_n_MW
    print(f"closing inverse-response peak: p_g={p_peak:.2f} MW at t={t_peak:.2f} s")
    print(f"final  (t={t[-1]:.1f} s): p_g={m.get_value('p_g_1') * S_n_MW:.2f} MW")


if __name__ == "__main__":
    main()
