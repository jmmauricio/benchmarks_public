"""
Kundur Two-Area, 4-generator test system — PES-TR18 variant (Case 6, with PSS).

Reference: PES-TR18 §5.3 and Appendix C (Lima & Silva, July 2015).

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

DATA = 'kundur2area.hjson'   # network description (PESTR18 §5.3, Case 6)


def build():
    grid = BpsBuilder(DATA, use_casadi=True)
    grid.uz_jacs = False
    grid.construct('kundur2area_pestr18')
    return CasadiBuilder(grid.sys_dict).build()


def ini():
    t0 = time.perf_counter()
    model = CasadiModel(build())
    model.ini({}, newton_tol=1e-13)

    report_buses(model, DATA)
    report_gens(model, DATA)
    validate_all(model, DATA)

    model.A_eval()
    ss = ssa.damp(model.A, model=model, sort='damp')

    eigs_ref = np.array([complex(e['real'], e['imag'])
                         for e in read_data(DATA)['results']['eigenvalues']])
    fig = ssa.plot_eig(ss['eigvalues'], x_min=-3.5, x_max=0.5, y_min=0, y_max=8,
                       fig='', mark='o', color='blue', label='')
    ssa.plot_eig(eigs_ref, x_min=-3.5, x_max=0.5, y_min=0, y_max=8,
                 fig=fig, mark='x', color='red', label='PESTR18 Table 5.30')
    fig.savefig('kundur2area_pestr18_eig.png')

    return model


def run(t_clear=1.2, bus_fault='08'):
    model = CasadiModel(build())
    model.Dt = 0.01
    model.ini({}, newton_tol=1e-8)
    model.run(1.0, {})

    # Short circuit at bus_fault: ramp shunt conductance, hold, then ramp clear.
    for i in range(1, 50):
        model.set_value(f"g_shunt_{bus_fault}", -i**2)
        model.recalculate_algebraics(tol=1e-8)
    model.run(t_clear, {})
    for i in range(1, 50):
        model.set_value(f"g_shunt_{bus_fault}", -(100 - i)**2)
        model.recalculate_algebraics(tol=1e-8)

    model.run(15, {})
    model.post()

    print('Run complete!')

    fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
    for g in range(1, 5):
        axes[0].plot(model.Time, model.get_values(f"omega_G{g:02}"), label=f"G{g:02}")
    for bus in range(5, 12):
        axes[1].plot(model.Time, model.get_values(f"V_{bus:02}"), label=f"V_{bus:02}")
    axes[0].set_ylabel("Speed (pu)"); axes[0].legend(); axes[0].grid(True)
    axes[1].set_ylabel("Voltage (pu)"); axes[1].legend(); axes[1].grid(True)
    axes[1].set_xlabel("Time (s)")
    fig.savefig('kundur2area_pestr18_run.png')

    return model


if __name__ == '__main__':
    ini()
    run()
