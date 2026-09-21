"""
Spanish equivalent system — residue-based PSS tuning.

    python pss_design.py                  # design for the default machine
    python pss_design.py nuclear_ES424    # ... or for any machine with a pss block

Method (see bloque_emec/residuos): the PSS closes a loop from its input signal
(rotor speed, for pss_kundur_2) back into the AVR summing junction.  The residue
of that loop for the target mode is

    R_k = (C v_k)(w_k^T B)

and the closed-loop eigenvalue shift for a small gain is dl_k = R_k H(l_k).  To
push the mode straight to the left the total phase must be 180 deg, so the
lead-lag network has to supply

    phi = 180 - arg(R_k) - arg(washout)
"""
import sys

import numpy as np

from pydae import ssa
from pydae.bps import BpsBuilder
from pydae.core.builder.casadi_builder import CasadiBuilder
from pydae.core.model.casadi_model import CasadiModel
from pydae.ssa import eig
from pydae.utils import read_data

DATA = 'esp_preblackout_dev.hjson'   # network description (buses, lines, generators, loads)

TARGET = 'CCGT_ES620'   # machine whose PSS is being tuned (override on the command line)
T_W = 10.0              # washout time constant, kept as in the data file
N_COMPS = 2             # lead-lag stages in series (pss_kundur_2 has two)
F_MIN, F_MAX = 0.1, 2.5 # band searched for the electromechanical mode [Hz]
ZETA_MAX = 0.10         # only modes below this damping ratio are worth a stabiliser
GAINS = (0.0, 0.01, 0.02, 0.03, 0.05, 0.08, 0.12, 0.2, 0.4, 0.8, 1.5, 3.0)

data = read_data(DATA)


def build():
    # Assemble the symbolic DAE for the Spanish pre-blackout system and fold it
    # into a CasADi SX graph (no C compilation required).
    grid = BpsBuilder(DATA, use_casadi=True)
    grid.uz_jacs = False            # skip the u/z Jacobians: only the A matrix is needed here
    grid.construct('esp_preblackout')
    return CasadiBuilder(grid.sys_dict).build()


def loop_signals(model):
    """Actuator column and sensor row of the loop the PSS actually closes.

    The stabiliser output v_pss enters the sexs summing junction as
    v_2 = v_ref - V + v_pss, i.e. through exactly the same path as v_ref, so the
    B column of v_ref is the actuator.  The sensor is whatever the stabiliser
    measures — for pss_kundur_2 that is the rotor speed, *not* the terminal
    voltage.
    """
    idx_u = model.u_run_list.index(f'v_ref_{TARGET}')
    B = np.array(model.B)[:, idx_u].reshape(-1, 1)

    C = np.zeros((1, model.N_x))
    C[0, model.x_list.index(f'omega_{TARGET}')] = 1.0
    return B, C


def pick_mode(model, B, C):
    """Poorly damped electromechanical mode this machine has most authority over.

    Two filters, both needed.  Damping, because a mode that is already well
    damped needs no stabiliser — and because the centring formula in design()
    places the maximum phase on the jw axis, which is only where the eigenvalue
    sits when it is lightly damped.  Then residue magnitude, because a badly
    damped mode the target machine can neither see nor move is not a candidate
    for its stabiliser either.
    """
    ev = model.eigenvalues
    V, W = model.right_eigenvectors, model.left_eigenvectors
    residues = {}
    for k in range(len(ev)):
        if ev[k].imag <= 0.05:
            continue                                  # skip real modes and conjugates
        if not (F_MIN <= ev[k].imag / (2 * np.pi) <= F_MAX):
            continue
        if -ev[k].real / abs(ev[k]) > ZETA_MAX:
            continue
        residues[k] = ((C @ V[:, [k]]) * (W[[k], :] @ B))[0, 0]

    if not residues:
        raise ValueError(f'no mode under zeta = {ZETA_MAX} found in '
                         f'{F_MIN}-{F_MAX} Hz: nothing here needs a stabiliser')

    print(f'\nmodes in {F_MIN}-{F_MAX} Hz under zeta = {100*ZETA_MAX:.0f} %, ranked by '
          f'residue magnitude (v_ref_{TARGET} -> omega_{TARGET})')
    print(f"{'idx':>5} {'lambda':>22} {'f [Hz]':>8} {'zeta [%]':>9} {'|R|':>11} {'arg R':>9}")
    ranked = sorted(residues, key=lambda k: -abs(residues[k]))
    for k in ranked[:8]:
        lam, R = ev[k], residues[k]
        print(f'{k:5d} {lam.real:9.4f}{lam.imag:+9.4f}j {lam.imag/(2*np.pi):8.3f} '
              f'{-100*lam.real/abs(lam):9.2f} {abs(R):11.4e} {np.angle(R, deg=True):9.2f}')

    k_t = ranked[0]
    return k_t, ev[k_t], residues[k_t]


def design(lam, R):
    """Lead-lag time constants and gain sign that make arg(R H) = 180 deg at s = lam."""
    G_w = (lam * T_W) / (1 + lam * T_W)               # washout evaluated at s = lam
    wrap = lambda a: (a + 180) % 360 - 180

    # This is the 180-degree rule.  Compensating arg(R) itself (instead of
    # 180 - arg(R)) rotates the eigenvalue shift the wrong way and can turn the
    # lead network into a lag network.
    #
    # A negative K_stab contributes 180 deg on its own, so the lead-lag chain
    # only has to cancel arg(R).  Both branches land on 180 deg and are equally
    # valid to first order, so the conventional one wins: positive gain with a
    # lead network.  Inverted polarity is the fallback for when that asks for
    # more phase than a lead-lag cascade can produce, which is what happens when
    # the residue already sits near 0 deg.
    base = np.angle(R, deg=True) + np.angle(G_w, deg=True)
    options = {+1: wrap(180.0 - base), -1: wrap(-base)}
    reachable = [s for s in (+1, -1) if abs(options[s] / N_COMPS) < 85.0]

    if not reachable:
        raise ValueError(f'{options[+1]/N_COMPS:.1f} deg/stage (or '
                         f'{options[-1]/N_COMPS:.1f} inverted) is out of reach for a '
                         f'lead-lag; use more stages than N_COMPS={N_COMPS}')

    sign = reachable[0]                               # +1 preferred, -1 as fallback
    phi = options[sign]
    phi_stage = phi / N_COMPS

    alpha = (1 + np.sin(np.deg2rad(phi_stage))) / (1 - np.sin(np.deg2rad(phi_stage)))
    T_2 = 1 / (lam.imag * np.sqrt(alpha))             # max phase centred on the mode
    T_1 = alpha * T_2

    G_c = ((1 + T_1 * lam) / (1 + T_2 * lam)) ** N_COMPS

    print(f'\ntarget mode        : {lam:.4f}  ({lam.imag/(2*np.pi):.3f} Hz, '
          f'zeta = {-100*lam.real/abs(lam):.2f} %)')
    print(f'residue            : |R| = {abs(R):.4f}, arg R = {np.angle(R, deg=True):.2f} deg')
    print(f'washout at s = lam : {np.angle(G_w, deg=True):.2f} deg')
    print(f'gain sign          : {"+" if sign > 0 else "-"} '
          f'(positive-gain branch would need {options[+1]:.2f} deg, '
          f'inverted branch {options[-1]:.2f} deg)')
    print(f'phase to supply    : {phi:.2f} deg  ({phi_stage:.2f} deg per stage), '
          f'alpha = {alpha:.3f}')
    print(f'G_c = ((1 + {T_1:.4f} s)/(1 + {T_2:.4f} s))^{N_COMPS}')
    achieved = np.angle(sign * R * G_w * G_c, deg=True)
    print(f'check              : arg(R H) = {achieved:.2f} deg (180 = pure damping)')
    if abs(abs(achieved) - 180) > 10:
        print('  WARNING: the compensator does not land on 180 deg at s = lam. The '
              'centring formula assumes lam sits near the jw axis.')

    return T_1, T_2, sign * G_w * G_c, sign


def sweep(model, lam, R, H1, k_t, gains, sign):
    """Move K_stab and watch where the mode really goes.

    No ini() here: the stabiliser output is zero at steady state, so the
    equilibrium does not depend on K_stab or on the lead-lag constants — only
    the Jacobian does.  (Calling ini() a second time on the same model is also
    not reproducible: it converges to a different operating point.)
    """
    ref = model.right_eigenvectors[:, k_t]
    ref = ref / np.linalg.norm(ref)

    print(f'\nclosed loop (first-order prediction: lambda + K dl/dK, '
          f'dl/dK = {abs(R*H1):.4f} < {np.angle(R*H1, deg=True):.2f} deg)')
    print(f"{'K_stab':>8} {'lambda':>22} {'f [Hz]':>8} {'zeta [%]':>9} "
          f"{'predicted':>22} {'corr':>6}")
    trace = []
    for K in gains:
        model.set_value(f'K_stab_pss_{TARGET}', sign * K)
        eig(model)
        ev, V = model.eigenvalues, model.right_eigenvectors
        osc = [k for k in range(len(ev)) if ev[k].imag > 0.05]
        corr = [abs(np.vdot(ref, V[:, k] / np.linalg.norm(V[:, k]))) for k in osc]
        lam_cl = ev[osc[int(np.argmax(corr))]]
        pred = lam + K * R * H1
        print(f'{sign*K:8.3f} {lam_cl.real:9.4f}{lam_cl.imag:+9.4f}j '
              f'{lam_cl.imag/(2*np.pi):8.3f} {-100*lam_cl.real/abs(lam_cl):9.2f} '
              f'{pred.real:11.4f}{pred.imag:+9.4f}j {max(corr):6.3f}')
        # Only trust the row while the tracker is still on the same mode: at high
        # gain the correlation decays and the eigenvalue can jump to a different
        # frequency altogether, which would otherwise look like a great result.
        held = max(corr) >= 0.8 and abs(lam_cl.imag - lam.imag) <= 0.3*lam.imag
        trace.append((K, -lam_cl.real/abs(lam_cl), held))

    # Smallest gain within 0.5 % damping of the best: past that the root locus
    # has bent over and extra gain only buys noise.
    kept = [(K, z) for K, z, held in trace if held]
    best = max(z for _, z in kept)
    K_knee = min(K for K, z in kept if z >= best - 0.005)
    z_knee = dict(kept)[K_knee]
    K_best = sign * K_knee
    lost = [K for K, _, held in trace if not held]
    if lost:
        print(f'\nmode tracking lost at |K_stab| >= {min(lost)}; those rows ignored')
    print(f'\nbest damping {100*best:.2f} %; knee of the locus at K_stab = {K_best} '
          f'({100*z_knee:.2f} %), within half a point of the best')
    return K_best


def ini():
    """Steady-state initialisation + small-signal analysis."""
    model = CasadiModel(build())

    # The residue has to be measured with this machine's loop open, so only the
    # stabiliser under design is switched off — the ones already tuned on other
    # machines stay in and count as part of the plant.
    model.set_value(f'K_stab_pss_{TARGET}', 0.0)
    model.ini({}, xy_0='xy_0.json')       # load-flow initialisation from a saved seed
    eig(model)                            # A, B, C, eigenvalues and eigenvectors

    ssa.damp(model.A, model=model, sort='damp',
             csv='esp_preblackout_eig.csv',
             html='esp_preblackout_eig.html')

    B, C = loop_signals(model)
    k_t, lam, R = pick_mode(model, B, C)
    T_1, T_2, H1, sign = design(lam, R)

    for name, value in (('T_1', T_1), ('T_3', T_1), ('T_2', T_2), ('T_4', T_2)):
        model.set_value(f'{name}_pss_{TARGET}', value)

    K_best = sweep(model, lam, R, H1, k_t, GAINS, sign)

    print(f'\nfor the "pss" block of {TARGET} in {DATA}:')
    print(f'"K_stab": {K_best}, "T_w": {T_W}, "T_1": {T_1:.4f}, "T_2": {T_2:.4f}, '
          f'"T_3": {T_1:.4f}, "T_4": {T_2:.4f},')

    return model


if __name__ == '__main__':
    if len(sys.argv) > 1:
        TARGET = sys.argv[1]
    print(f'tuning the PSS of {TARGET}')
    ini()       # steady-state initialization + small-signal analysis
