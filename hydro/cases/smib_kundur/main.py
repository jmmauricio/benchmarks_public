"""Hydro SMIB step-response demo — gensal + HYGOV + SEXS.

Loads the case, inits, applies a p_c_lc setpoint step 0.80 → 0.95 → 0.80,
plots p_g, p_m, omega, gate, water flow, terminal voltage. The plot
shows the characteristic non-minimum-phase (inverse) response of a
hydro turbine: opening the gate momentarily reduces output before the
water column accelerates.
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
    assert ok, "ini failed — check seeds"

    print(f"initial: p_g={m.get_value('p_g_1'):.4f}, "
          f"x_g={m.get_value('x_g_gov_1'):.4f}, "
          f"V_1={m.get_value('V_1'):.4f}")

    m.run( 5.0, {})                          # 5 s settle
    m.run(60.0, {"p_c_lc_1": 0.95})          # step up
    m.run(60.0, {"p_c_lc_1": 0.80})          # step back
    m.post()

    t = m.Time
    S_n_MW = 250.0

    fig, axes = plt.subplots(4, 1, figsize=(9, 9), sharex=True)
    axes[0].plot(t, m.get_values("p_g_1") * S_n_MW, label="p_g", color="C0")
    axes[0].plot(t, m.get_values("p_m_1") * S_n_MW, label="p_m", color="C1", ls="--")
    axes[0].plot(t, m.get_values("p_c_lc_1") * S_n_MW, label="p_c_lc",
                 color="k", ls=":")
    axes[0].set_ylabel("Power (MW)")
    axes[0].grid(); axes[0].legend()

    axes[1].plot(t, (m.get_values("omega_1") - 1.0) * 100.0, color="C2")
    axes[1].set_ylabel(r"$(\omega-1)\times 100$ (%)")
    axes[1].grid()

    axes[2].plot(t, m.get_values("x_g_gov_1"), label=r"gate $x_g$", color="C3")
    axes[2].plot(t, m.get_values("q_gov_1"),   label=r"flow $q$",
                 color="C4", ls="--")
    axes[2].set_ylabel("Gate / flow")
    axes[2].grid(); axes[2].legend()

    axes[3].plot(t, m.get_values("V_1"), color="C5")
    axes[3].set_ylabel(r"$V_1$ (pu)")
    axes[3].set_xlabel("Time (s)")
    axes[3].grid()

    fig.suptitle("Hydro SMIB — gensal + HYGOV + SEXS, "
                 "p_c_lc step 0.80 → 0.95 → 0.80")
    fig.tight_layout()
    out = here / "hydro_smib_step.png"
    fig.savefig(out, dpi=130)
    print(f"saved {out}")


if __name__ == "__main__":
    main()
