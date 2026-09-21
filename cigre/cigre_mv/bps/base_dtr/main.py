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
from pydae.bps.utils.dtr_report import report_dtr
from pydae.bps.utils.validator import validate_all
from pydae.utils import read_data

DATA = 'cigre_mv.hjson'   # network description (buses, lines, generators, loads, reference results)


def build():
    # Assemble the symbolic DAE for the IEEE 39-bus system from the HJSON description.
    grid = BpsBuilder(DATA, use_casadi=True)
    grid.uz_jacs = False            # skip the u/z Jacobians: only the A matrix is needed here
    grid.construct('cigre_mv')        # concatenate every component's equations into grid.sys_dict
    # Fold sys_dict into a CasADi SX graph (no C compilation required).
    return CasadiBuilder(grid.sys_dict).build()


def ini():
    t0 = time.perf_counter()
    model = CasadiModel(build())            # runtime model wrapping the CasADi graph
    model.ini({}, newton_tol=1e-13)         # Newton-Raphson load-flow initialization

    # Initialization report: bus voltages/angles and generator dispatch.
    report_buses(model, DATA)
    report_gens(model, DATA)
    report_dtr(model, DATA)              # DTR conductor temperature & heat balance
    # Validation: compare the solved steady state against the reference results in the HJSON.
    validate_all(model, DATA)

    # Small-signal analysis.
    model.A_eval()                                   # reduced state matrix A = Fx - Fy*inv(Gy)*Gx
    ss = ssa.damp(model.A, model=model, sort='damp')  # eigenvalues, damping ratios and frequencies

    # Plot the computed eigenvalues; overlay reference modes if the HJSON ships a
    # results.eigenvalues section (cigre_mv has none, so the overlay is skipped).
    fig = ssa.plot_eig(ss['eigvalues'], fig='', mark='o', color='blue', label='computed')
    results = read_data(DATA).get('results', {})
    if 'eigenvalues' in results:
        eigs_ref = np.array([complex(e['real'], e['imag']) for e in results['eigenvalues']])
        ssa.plot_eig(eigs_ref, fig=fig, mark='x', color='red', label='reference')
    fig.savefig('cigre_mv_eig.png')

    model.report_u()

    model.report_z()

    return model


def run():
    data = read_data(DATA)
    model = CasadiModel(build())
    model.Dt = 10
    model.ini({}, xy_0='xy_0.json')      # load-flow initialization
    model.run(3600.0, {})                  # 1 s of steady operation before the disturbance
    model.run(10*3600.0, {'p_p_B03': 5/25})                  # 1 s of steady operation before the disturbance
    model.run(12*3600.0, {'p_p_B03': 4.9/25})                  # 1 s of steady operation before the disturbance

    # Short circuit at bus_fault: ramp the shunt conductance up to a large magnitude so t

    model.post()                        # copy the solver buffers into the public Time/X/Y/Z arrays

    # Plot generator speeds (top) and bus voltages 30-39 (bottom).
    fig, axes = plt.subplots(1, 1, figsize=(8, 7), sharex=True)
    for line in data['lines']:
        if 'dtr' in line:
            line_name = f'{line['bus_j']}_{line['bus_k']}'
            axes.plot(model.Time/60, model.get_values(f"t_conductor_{line_name}"), label=f"{line_name}")
    
    fig.savefig('temperatures.png')

    return model


if __name__ == '__main__':
    ini()   # steady-state initialization + small-signal analysis
    run()   # time-domain short-circuit simulation
