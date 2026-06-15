"""
NTS base case — CasADi pipeline.

    python main.py     # ini:   load flow + report + validation + small-signal
                       # run:   v_ref_4 step-change response
                       # sweep: bus 2-3 X_L sweep, overlay NTS Figura 18
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
from pydae.control.continuos.stab_1lpf_2wo_3ll import stab_1lpf_2wo_3ll
from pydae.control.continuos.pi import pi

DATA = 'nts_mpe_wecc.hjson'     # network description with buses, generators, lines and reference results
XY_0 = 'nts_mpe_wecc_xy_0.json'  # saved initial guess for the Newton-Raphson solver


def build():
    # Assemble the CasADi DAE for the NTS base case from the HJSON description.
    grid = BpsBuilder(DATA, use_casadi=True)
    grid.construct('nts_mpe_wecc')          # concatenate all component equations into grid.sys_dict

    # # POD: 2-washout/3-lead-lag stabilizer filter chain (pydae.control.continuos),
    # # input = omega_1 (generator 1 speed), output exposed as P_5 (power-modulation signal).
    # pod_data = {'T_lpf': 0.02, 'T_wo1': 10.0, 'T_wo2': 10.0,
    #              'T_1': 0.1, 'T_2': 0.02, 'T_3': 0.1, 'T_4': 0.02, 'T_5': 0.1, 'T_6': 0.02,
    #              'K_stab': 10, 'Limit': 0.1}
    # stab_1lpf_2wo_3ll(grid.sys_dict, 'pod_5', grid.backend, pod_data, input='omega_1', output='P_5')

    # POI active-/reactive-power regulators driving the MPE WECC stack.
    #   • Inputs are in W / VAR (S_base = 1e9, so 1350 MW = 1.35e9, etc.).
    #   • Input order is [reference, measurement] so the PI error
    #     u = p_ref − p_meas is positive when more output is needed; paired
    #     with positive K_p / K_i this is a negative-feedback regulator.
    #     (The reverse order [meas, ref] saturates the integrator and
    #     gives a non-regulating fixed point at ini.)
    #   • P PI output is Pref_5 (REEC_B Pref, pu on 1500 MVA): natural
    #     range [0, 1], so z_max ≥ 0.9 is required to deliver 1350 MW.
    #   • Q PI output is Qext_5 (pu on 1500 MVA): bounded to [−0.5, 1.1].
    pi(grid.sys_dict, 'poi_p', grid.backend,
       {'K_p': 1e-9, 'K_i': 1e-10, 'K_aw': 2.0, 'z_min': 0.0,  'z_max': 1.0, 'u': 0.0, 'xi_0': 0.0},
       input=['p_poi_ref', 'p_line_5_2'], output='Pref_5')

    pi(grid.sys_dict, 'poi_q', grid.backend,
       {'K_p': 1e-9, 'K_i': 1e-10, 'K_aw': 2.0, 'z_min': -0.5, 'z_max': 1.1, 'u': 0.0, 'xi_0': 0.0},
       input=['q_poi_ref', 'q_line_5_2'], output='Qext_5')

    return CasadiBuilder(grid.sys_dict).build()  # fold sys_dict into a CasADi SX graph (no C compile)


def ini():
    model = CasadiModel(build())        # runtime model wrapping the CasADi graph
    model.decimation = 10               # store every 10th integration step
    # Adjust line 2-3 impedance before initialization.
    change_line(model, {"bus_j": "2", "bus_k": "3",
                        "X_pu": 0.6, "R_pu": 0.0, "Bs_pu": 0.0, "S_mva": 100})

    model.ini({'p_poi_ref':1250e6, 'q_poi_ref':0.1e9}, XY_0)                                 # Newton-Raphson load-flow initialization

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
                        "X_pu": 0.01, "R_pu": 0.0, "Bs_pu": 0.0, "S_mva": 100})

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


def sweep():
    # X_L sweep on the bus 2-3 line. For each value we re-initialize the model
    # and pick the dominant electromechanical mode (max participation on
    # delta_1, narrowed to 0.1-1.2 Hz). The locus is then overlaid with the
    # NTS Figura 18 read-off stored in the HJSON.
    model = CasadiModel(build())
    model.decimation = 10

    nts_ref = read_data(DATA)['results']['eigenvalues_XL_sweep_nts_fig18']
    X_L_values = [entry['X_L'] for entry in nts_ref]

    pydae_eigs = []
    for X_L in X_L_values:
        change_line(model, {"bus_j": "2", "bus_k": "3",
                            "X_pu": X_L, "R_pu": 0.0, "Bs_pu": 0.0, "S_mva": 100})
        model.ini({}, XY_0)
        ssa.eig(model)                                            # populates model.eigenvalues
        _, lam = ssa.get_mode(model, f_min=0.1, f_max=1.2, report=False)
        pydae_eigs.append(lam)
        print(f"X_L = {X_L:.2f}  →  λ = {lam.real:+.4f} {lam.imag:+.4f}j  "
              f"(f = {abs(lam.imag)/(2*np.pi):.3f} Hz, "
              f"ζ = {-lam.real/abs(lam)*100:.2f}%)")

    pydae_eigs = np.array(pydae_eigs)
    nts_eigs   = np.array([complex(e['real'], e['imag']) for e in nts_ref])

    # Plot the locus on the complex S-plane: damping zones (green/orange/red)
    # are added by ssa.plot_eig; pydae as blue 'o', NTS Figura 18 as red 'x'.
    fig = ssa.plot_eig(pydae_eigs, x_min=-2.2, x_max=0.1, y_min=0.0, y_max=1.1,
                       fig='', mark='o', color='blue', label='')
    ssa.plot_eig(nts_eigs, x_min=-2.2, x_max=0.1, y_min=0.0, y_max=1.1,
                 fig=fig, mark='x', color='red', label='')

    ax = fig.axes[0]
    # Connect the points along the sweep direction so the trajectory is clear.
    ax.plot(pydae_eigs.real, pydae_eigs.imag / (2 * np.pi), '-',
            color='blue', alpha=0.4)
    ax.plot(nts_eigs.real, nts_eigs.imag / (2 * np.pi), '--',
            color='red', alpha=0.4)
    # Proxy artists so the legend shows the combined marker + line style.
    ax.plot([], [], 'o-',  color='blue', label='pydae (genrou)')
    ax.plot([], [], 'x--', color='red',  label='NTS Figura 18')
    ax.legend(loc='upper left')
    ax.set_title('Bus 2-3 line $X_L$ sweep — dominant electromechanical mode')

    fig.savefig('nts_base_eig_sweep.png', dpi=300)

    return model


if __name__ == "__main__":
    model = ini()    # steady-state initialization + small-signal analysis
    model.report_u()
    model.report_z()
    # run()    # time-domain v_ref_4 step-change simulation
    # sweep()  # X_L sweep, overlay with NTS Figura 18 read-off
