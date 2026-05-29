"""
NTS base case — CasADi pipeline.

    python main.py     # ini: load flow + report + validation + small-signal
                       # run: v_ref_4 step-change response
"""
import time

import numpy as np
from matplotlib import pyplot as plt

from pydae import ssa
from pydae.bps import BpsBuilder
from pydae.core.builder import CasadiBuilder
from pydae.core.model import CasadiModel
from pydae.bps.utils.reporter import report_all
from pydae.bps.utils.validator import validate_all
from pydae.bps.lines import change_line
from pydae.utils import read_data

DATA = 'nts_base.hjson'     # network description with buses, generators, lines and reference results
XY_0 = 'nts_base_xy_0.json'  # saved initial guess for the Newton-Raphson solver


def build():
    # Assemble the CasADi DAE for the NTS base case from the HJSON description.
    grid = BpsBuilder(DATA, use_casadi=True)
    grid.construct('nts_base')          # concatenate all component equations into grid.sys_dict
    return CasadiBuilder(grid.sys_dict).build()  # fold sys_dict into a CasADi SX graph (no C compile)


def ini():
    model = CasadiModel(build())        # runtime model wrapping the CasADi graph
    model.decimation = 10               # store every 10th integration step
    # Adjust line 2-3 impedance before initialization.
    change_line(model, {"bus_j": "2", "bus_k": "3",
                        "X_pu": 0.6, "R_pu": 0.0, "Bs_pu": 0.0, "S_mva": 100})

    model.ini({}, XY_0)                 # Newton-Raphson load-flow initialization

    # Initialization report: bus voltages/angles, generator dispatch and line flows.
    report_all(model, DATA)
    # Validation: compare the solved state against the reference values in the HJSON.
    validate_all(model, DATA)

    # Small-signal analysis.
    model.A_eval()                                    # reduced state matrix A = Fx - Fy*inv(Gy)*Gx
    ssa.damp(model.A, model=model, sort='damp')       # eigenvalues, damping ratios and frequencies
    ssa.eig(model)                                    # eigenvectors and participation factors
    ssa.get_mode(model, f_min=0.1, f_max=0.5)        # identify inter-area modes in 0.1-0.5 Hz band

    return model


def run():
    model = CasadiModel(build())        # runtime model wrapping the CasADi graph
    model.decimation = 10               # store every 10th integration step
    # Adjust line 2-3 impedance before initialization.
    change_line(model, {"bus_j": "2", "bus_k": "3",
                        "X_pu": 0.6, "R_pu": 0.0, "Bs_pu": 0.0, "S_mva": 100})

    model.ini({}, XY_0)                 # load-flow initialization
    model.run(2.0, {})                  # 2 s of steady operation before the disturbance

    # Step change in the voltage reference of generator 4 (+0.018 pu).
    model.run(40.0, {"v_ref_4": model.get_value('v_ref_4') + 0.018})
    model.post()                        # copy solver buffers into the public Time/X/Y/Z arrays

    # Reference data from the HJSON: columns [t (s), p_g_1 (MW)].
    nts_results = np.array(read_data(DATA)['results']['step_vref4']['data'])

    # Plot generator speed (top) and active power vs. NTS reference (bottom).
    fig, axes = plt.subplots(2, 1, figsize=(8, 8), sharex=True)

    axes[0].plot(model.Time, model.get_values("omega_1"), label="omega_1", color="b")
    axes[0].set_ylabel("Speed (pu)"); axes[0].legend(); axes[0].grid(True)

    axes[1].plot(model.Time, model.get_values("p_g_1") * model.get_value('S_n_1') / 1e6,
                 label="p_g_1", color="b")
    axes[1].plot(nts_results[:, 0], nts_results[:, 1], label="p_g_1 (NTS ref)", color="r")
    axes[1].set_ylabel("Power (MW)"); axes[1].legend(); axes[1].grid(True)
    axes[1].set_ylim((1330, 1375))
    axes[1].yaxis.set_major_locator(plt.MultipleLocator(5))
    axes[1].set_xlabel("Time (s)")

    fig.savefig("nts_base.png", dpi=300)

    return model


if __name__ == "__main__":
    ini()   # steady-state initialization + small-signal analysis
    run()   # time-domain v_ref_4 step-change simulation
