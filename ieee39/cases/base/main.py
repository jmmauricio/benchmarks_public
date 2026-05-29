"""
IEEE 39-bus New England test system — CasADi pipeline.

    python main.py        # ini: load flow + report + validation + small-signal
"""
import os
import time

import numpy as np
from matplotlib import pyplot as plt

from pydae import ssa
from pydae.bps import BpsBuilder
from pydae.core.builder.casadi_builder import CasadiBuilder
from pydae.core.model.casadi_model import CasadiModel
from pydae.bps.utils.reporter import report_buses, report_gens
from pydae.bps.utils.validator import validate_all
from pydae.utils import read_data

DATA = 'ieee39.hjson'   # network description (buses, lines, generators, loads, reference results)


def build():
    # Assemble the symbolic DAE for the IEEE 39-bus system from the HJSON description.
    grid = BpsBuilder(DATA, use_casadi=True)
    grid.uz_jacs = False            # skip the u/z Jacobians: only the A matrix is needed here
    grid.construct('ieee39')        # concatenate every component's equations into grid.sys_dict
    # Fold sys_dict into a CasADi SX graph (no C compilation required).
    return CasadiBuilder(grid.sys_dict).build()


def ini():
    t0 = time.perf_counter()
    model = CasadiModel(build())            # runtime model wrapping the CasADi graph
    model.ini({}, newton_tol=1e-13)         # Newton-Raphson load-flow initialization

    # Initialization report: bus voltages/angles and generator dispatch.
    report_buses(model, DATA)
    report_gens(model, DATA)
    # Validation: compare the solved steady state against the reference results in the HJSON.
    validate_all(model, DATA)

    # Small-signal analysis.
    model.A_eval()                                   # reduced state matrix A = Fx - Fy*inv(Gy)*Gx
    ss = ssa.damp(model.A, model=model, sort='damp')  # eigenvalues, damping ratios and frequencies

    # Overlay the computed eigenvalues (blue) with the reference electromechanical modes (red).
    eigs_ref = np.array([complex(e['real'], e['imag'])
                         for e in read_data(DATA)['results']['eigenvalues']])
    fig = ssa.plot_eig(ss['eigvalues'], x_min=-3.5, x_max=-1, y_min=0, y_max=2,
                       fig='', mark='o', color='blue', label='')
    ssa.plot_eig(eigs_ref, x_min=-3.5, x_max=-1, y_min=0, y_max=2,
                 fig=fig, mark='x', color='red', label='Selected modes')
    fig.savefig('ieee39_eig.png')

    return model


def run(t_clear=1.2, bus_fault='16'):
    model = CasadiModel(build())
    model.Dt = 0.01                     # integration step size (s)
    model.ini({}, newton_tol=1e-8)      # load-flow initialization
    model.run(1.0, {})                  # 1 s of steady operation before the disturbance

    # Short circuit at bus_fault: ramp the shunt conductance up to a large magnitude so the
    # algebraic recalculation keeps converging, hold through the fault, then ramp it back to clear.
    for i in range(1, 50):
        model.set_value(f"g_shunt_{bus_fault}", -i**2)
        model.recalculate_algebraics(tol=1e-8)
    model.run(t_clear, {})              # integrate through the fault until the clearing time
    for i in range(1, 50):
        model.set_value(f"g_shunt_{bus_fault}", -(100 - i)**2)
        model.recalculate_algebraics(tol=1e-8)

    model.run(10, {})                   # post-fault response up to t = 10 s
    model.post()                        # copy the solver buffers into the public Time/X/Y/Z arrays

    print('Run complet!')

    # Plot generator speeds (top) and bus voltages 30-39 (bottom).
    fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
    for g in range(1, 11):
        axes[0].plot(model.Time, model.get_values(f"omega_G{g:02}"), label=f"G{g:02}")
    for bus in range(30, 40):
        axes[1].plot(model.Time, model.get_values(f"V_{bus:02}"), label=f"V_{bus:02}")
    axes[0].set_ylabel("Speed (pu)"); axes[0].legend(); axes[0].grid(True)
    axes[1].set_ylabel("Voltage (pu)"); axes[1].legend(); axes[1].grid(True)
    axes[1].set_xlabel("Time (s)")
    fig.savefig('ieee39_run.png')

    return model


if __name__ == '__main__':
    ini()   # steady-state initialization + small-signal analysis
    run()   # time-domain short-circuit simulation
