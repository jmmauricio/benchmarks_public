"""
68-bus (NETS/NYPS), 16-generator, 5-area benchmark — CasADi pipeline.

Reference: PES-TR18 §5.6 (IEEE PES Task Force on Benchmark Systems for
Stability Controls, August 2015).

    python main.py        # ini: load flow + report + validation + small-signal
"""
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

DATA = 'system68bus.hjson'


def build():
    grid = BpsBuilder(DATA, use_casadi=True)
    grid.uz_jacs = False
    grid.construct('system68bus')
    return CasadiBuilder(grid.sys_dict).build()


def ini():
    t0 = time.perf_counter()
    model = CasadiModel(build())
    model.ini({}, newton_tol=1e-10)

    report_buses(model, DATA)
    report_gens(model, DATA)
    validate_all(model, DATA)

    model.A_eval()
    ss = ssa.damp(model.A, model=model, sort='damp')

    # PESTR18 §5.6.3 publishes no numerical eigenvalue reference for this case
    # (only figures).  Plot the computed eigenvalues; user can overlay their own
    # reference by editing read_data(DATA)['results']['eigenvalues'].
    eigs_ref = np.array([complex(e['real'], e['imag'])
                         for e in read_data(DATA)['results']['eigenvalues']])
    fig = ssa.plot_eig(ss['eigvalues'], x_min=-6.0, x_max=1.0, y_min=0, y_max=12,
                       fig='', mark='o', color='blue', label='')
    if eigs_ref.size:
        ssa.plot_eig(eigs_ref, x_min=-6.0, x_max=1.0, y_min=0, y_max=12,
                     fig=fig, mark='x', color='red', label='reference')
    fig.savefig('system68bus_eig.png')

    print(f"ini() done in {time.perf_counter() - t0:.1f} s")
    return model


def run(t_clear=1.2, bus_fault='60'):
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

    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    for g in range(1, 17):
        axes[0].plot(model.Time, model.get_values(f"omega_G{g:02}"), label=f"G{g:02}")
    for bus in [27, 53, 54, 60, 61, 17, 18]:    # tie-line buses + load buses
        axes[1].plot(model.Time, model.get_values(f"V_{bus:02}"), label=f"V_{bus:02}")
    axes[0].set_ylabel("Speed (pu)"); axes[0].legend(ncol=2, fontsize=7); axes[0].grid(True)
    axes[1].set_ylabel("Voltage (pu)"); axes[1].legend(ncol=2, fontsize=7); axes[1].grid(True)
    axes[1].set_xlabel("Time (s)")
    fig.savefig('system68bus_run.png')

    return model


if __name__ == '__main__':
    ini()
    run()
