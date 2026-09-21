"""
Spanish equivalent system — CasADi pipeline.

    python main.py        # step the interconnection power and observe the inter-area response
"""
import time

import hjson
import numpy as np
from matplotlib import pyplot as plt

from pydae import ssa
from pydae.bps import BpsBuilder
from pydae.core.builder.casadi_builder import CasadiBuilder
from pydae.core.model.casadi_model import CasadiModel
from pydae.bps.utils.reporter import report_buses, report_gens
from pydae.utils import read_data

DATA = 'esp_preblackout_dev.hjson'   # network description (buses, lines, generators, loads, reference results)


data = read_data(DATA)

def build():
    # Assemble the symbolic DAE for the Spanish pre-blackout system and fold it
    # into a CasADi SX graph (no C compilation required).
    grid = BpsBuilder(DATA, use_casadi=True)
    grid.uz_jacs = False            # skip the u/z Jacobians: only the A matrix is needed here
    grid.construct('esp_preblackout')
    return CasadiBuilder(grid.sys_dict).build()


def ini():
    """Steady-state initialisation + small-signal analysis."""
    model = CasadiModel(build())

    CCGTs = [ 'CCGT_ES213', 'CCGT_ES512', 'CCGT_ES620' ]
    CCGTs = [ 'CCGT_ES620' ]

    params = {}
    #for ccgt in CCGTs:
    #    params.update({f'p_c_lc_{ccgt}': 0.5})  #  
    model.ini(params, xy_0='xy_0.json')  # load-flow initialisation from a saved seed

    report_buses(model, DATA)
    report_gens(model, DATA)

    model.A_eval()                                    # reduced state matrix A = Fx - Fy*inv(Gy)*Gx
    ss = ssa.damp(model.A, model=model, sort='damp', 
                  csv='esp_preblackout_eig.csv',
                  html='esp_preblackout_eig.html')  # eigenvalues, damping ratios and frequencies
    fig = ssa.plot_eig(ss['eigvalues'], x_min=-6, x_max=0.2, y_min=0, y_max=2.5,
                       fig='', mark='o', color='blue', label='Spain pre-blackout')
    fig.savefig('esp_preblackout_eig.png')

    # model.report_u()

    return model


def run():
    """Short time-domain run from a saved operating point, comparing solvers.

    For each method in ``methods`` the model is re-initialised from ``xy_0.json``
    and integrated to ``t_end``; wall time and the final state are reported so the
    fixed-step ``trapezoidal`` solver can be checked against ``idas_grid``.
    """
    t0 = time.perf_counter()
    builder = build()
    print(f"instantiation done in {time.perf_counter() - t0:.2f} s")

    model = CasadiModel(builder)
    model.Dt = 0.1                         # output step size (s)
    model.ini({}, xy_0='xy_0.json')        # load-flow initialisation from a saved seed
    ta = time.perf_counter()
    model.run(1.0, method='trapezoidal')        # time-domain integration

    CCGTs = [ 'CCGT_ES620' ]
    params = {}
    for ccgt in CCGTs:
        params.update({f'dp_lc_{ccgt}': -0.5})  #  

    model.run(30.0, params, method='trapezoidal')        # time-domain integration

    model.post()

    return model


def run_link(dP_pu=-0.20, machine='CCGT_ES620', t_settle=1.0, t_end=30.0, Dt=0.05,
             compare=True):
    """Step the Spain->France interconnection power and watch the 0.2 Hz inter-area
    response, integrated with the built-in fixed-step ``trapezoidal`` solver.

    The step is applied as a load-control command change ``dp_lc_<machine>`` on the
    largest CCGT; the France slack takes up the mismatch, exciting the inter-area
    mode. With ``compare=True`` the same disturbance is also run with the ten
    largest units reset to the generic (untuned) PSS, so the tuned response can be
    overlaid on the untuned/unstable one.
    """
    spain = [s for s in data['syns'] if s['name'] != 'SLACK_FR']
    S_FR = next(s['S_n'] for s in data['syns'] if s['name'] == 'SLACK_FR')
    H = {s['name']: s['H'] for s in spain}
    targets = [s['name'] for s in sorted(spain, key=lambda s: -s['S_n'])[:10]]

    def scenario(generic):
        model = CasadiModel(build())
        model.Dt = Dt
        model.ini({}, xy_0='xy_0.json')                 # simple seeded init, no warm-up
        if generic:                                     # revert the ten tuned units to the default PSS
            for nm in targets:
                model.set_value(f'K_stab_pss_{nm}', 1.0)
                model.set_value(f'T_1_pss_{nm}', 1.0); model.set_value(f'T_2_pss_{nm}', 0.05)
                model.set_value(f'T_3_pss_{nm}', 3.0); model.set_value(f'T_4_pss_{nm}', 0.5)
        model.run(t_settle, method='trapezoidal')
        try:
            model.run(t_end, {f'dp_lc_{machine}': dP_pu}, method='trapezoidal')
        except RuntimeError as e:
            print(f"  [{'generic' if generic else 'tuned'}] integration stopped early: {e}")
        model.post()
        t = np.array(model.Time)
        p_tie = -np.array(model.get_values('p_g_SLACK_FR')) * S_FR / 1e6
        w = np.array([model.get_values(f"omega_{s['name']}") for s in spain])   # (n_spain, nt)
        f_es = np.average(w, axis=0, weights=[H[s['name']] for s in spain])     # inertia-weighted COI
        f_fr = np.array(model.get_values('omega_SLACK_FR'))
        return t, p_tie, f_es, f_fr

    print(f"run_link: {dP_pu:+.2f} pu step at {machine}, Dt={Dt}s, t_end={t_end}s")
    t0 = time.perf_counter()
    tt, pt, fe, ff = scenario(generic=False)
    print(f"  tuned run: {len(tt)} points in {time.perf_counter() - t0:.1f} s")

    fig, ax = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    if compare:
        tg, pg, ge, gf = scenario(generic=True)
        ax[0].plot(tg, pg, 'r--', lw=1.1, label='generic PSS (untuned)')
        ax[1].plot(tg, (ge - 1) * 50e3, 'r--', lw=1.1, label='Spain COI (generic)')
    ax[0].plot(tt, pt, 'g-', lw=1.6, label='tuned PSS')
    ax[0].set_ylabel('Spain $\\to$ France tie [MW]'); ax[0].grid(True); ax[0].legend(loc='best')
    ax[0].set_title(f'Interconnection power step {dP_pu:+.2f} pu at {machine}  (trapezoidal, Dt={Dt}s)')
    ax[1].plot(tt, (fe - 1) * 50e3, 'g-', lw=1.6, label='Spain COI (tuned)')
    ax[1].plot(tt, (ff - 1) * 50e3, 'navy', lw=1.0, label='France (slack)')
    ax[1].set_ylabel('Freq. deviation [mHz]'); ax[1].set_xlabel('Time [s]')
    ax[1].grid(True); ax[1].legend(loc='best')
    fig.tight_layout(); fig.savefig('esp_preblackout_link.png')
    print('saved esp_preblackout_link.png')
    return fig


def run_overvoltage(trip_frac=0.05, t_settle=1.0, t_end=10.0, Dt=0.01):
    """Single-event overvoltage study: trip shunt reactors and watch bus voltages rise.

    The 400 kV lines inject ~390 Mvar of charging; the shunt reactors (negative
    susceptance ``b_shunt`` < 0) absorb it. Disconnecting them (``b_shunt`` -> 0)
    dumps that reactive onto the network. The sexs AVRs cut excitation toward
    E_min = 0 and then saturate, so the voltage settles at an elevated level --
    the fundamental-frequency overvoltage mechanism of the 28-Apr-2025 event.

    A pure fixed-parameter step (no protection): shows the *initial* voltage rise.
    """
    model = CasadiModel(build())
    model.Dt = Dt
    model.ini({}, xy_0='xy_0.json')                     # simple seeded init

    region_buses = [b['name'] for b in data['buses']
                    if not b['name'].startswith('GT_') and b['name'] != 'FR']
    # shunt reactors: inductive (b_shunt < 0), i.e. reactive absorbers
    reactors = []
    for b in data['buses']:
        try:
            val = float(model.get_value(f"b_shunt_{b['name']}"))
        except Exception:
            continue
        if val < -1e-9:
            reactors.append((b['name'], val))
    reactors.sort(key=lambda x: x[1])                   # most absorbing first
    n = max(1, int(round(len(reactors) * trip_frac)))
    tripped = reactors[:n]
    print(f"tripping {n}/{len(reactors)} shunt reactors at t={t_settle}s")

    model.run(t_settle, method='trapezoidal')
    for bus, _ in tripped:
        model.set_value(f"b_shunt_{bus}", 0.0)          # disconnect the reactors
    model.recalculate_algebraics(tol=1e-8)              # algebraic voltage jump at the event
    collapsed = False
    try:
        model.run(t_end, method='trapezoidal')          # dynamic response (AVRs relax to E_min)
    except RuntimeError as e:                            # voltage runaway -> solver can't continue
        collapsed = True
        print(f"  VOLTAGE INSTABILITY: integration stopped early -- {e}")
    model.post()

    t = np.array(model.Time)
    Vmat = np.array([model.get_values(f'V_{b}') for b in region_buses])
    vmax, vmean = Vmat.max(axis=0), Vmat.mean(axis=0)
    vend = Vmat[:, -1]
    order = np.argsort(vend)[::-1]
    print("highest final bus voltages [pu]:",
          [(region_buses[i], round(float(vend[i]), 3)) for i in order[:6]])
    print(f"pre-event max V = {vmax[t < t_settle].max():.3f} pu -> post-event max V = {vmax[-1]:.3f} pu")

    fig, ax = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    for i in order[:5]:
        ax[0].plot(t, Vmat[i], lw=1, alpha=0.8, label=region_buses[i])
    ax[0].plot(t, vmax, 'r-', lw=2, label='max bus V')
    ax[0].plot(t, vmean, 'b-', lw=1.2, label='mean bus V')
    ax[0].axhline(1.10, color='orange', ls='--', lw=1, label='1.10 pu')
    ax[0].axhline(1.15, color='red', ls='--', lw=1, label='1.15 pu')
    ax[0].set_ylabel('Bus voltage [pu]'); ax[0].grid(True); ax[0].legend(loc='best', fontsize=8)
    if collapsed:
        ax[0].axvline(t[-1], color='k', ls='-', lw=1.2)
        ax[0].annotate('voltage collapse\n(solver stops)', (t[-1], vmax.max()),
                       textcoords='offset points', xytext=(-70, -10), fontsize=9, color='k')
    status = f'VOLTAGE COLLAPSE at t={t[-1]:.2f}s' if collapsed else 'settles (stable)'
    ax[0].set_title(f'Overvoltage: {n} shunt reactors tripped at t={t_settle}s  ->  {status}')
    # AVR field voltages driven toward the E_min = 0 floor
    for s in data['syns'][:8]:
        try:
            ax[1].plot(t, model.get_values(f"v_f_{s['name']}"), lw=1, label=s['name'])
        except Exception:
            pass
    ax[1].axhline(0.0, color='k', ls=':', lw=1, label='E_min = 0')
    ax[1].set_ylabel('AVR field voltage v_f [pu]'); ax[1].set_xlabel('Time [s]')
    ax[1].grid(True); ax[1].legend(loc='best', fontsize=7)
    fig.tight_layout(); fig.savefig('esp_preblackout_overvoltage.png')
    print('saved esp_preblackout_overvoltage.png')
    return model


def enable_qv_droop(model, K_qv=3.0):
    """Turn on the grid-code Q-V droop on every solar (pv_pq_ss) inverter.

    Each inverter's droop reference ``V_qref`` is seeded to its present terminal
    voltage, so the droop term ``K_qv*(V_qref - V_m)`` is zero at the operating
    point (the load flow / init seed is preserved). Above ``V_qref`` the inverter
    absorbs reactive power (voltage support); when it is tripped that absorption
    disappears -- the reactive-balance actuator of the 28-Apr-2025 event.

    Returns the list of solar unit dicts that received the droop.
    """
    solar = [p for p in data['pvs'] if p.get('type') == 'pv_pq_ss']
    for p in solar:
        V = float(model.get_value(f"V_{p['bus']}"))
        model.set_value(f"V_qref_{p['name']}", V)
        model.set_value(f"K_qv_{p['name']}", K_qv)
    return solar


def run_gentrip_cascade(K_qv=1.0, n_seed_reactors=4, V_trip=1.045, V_trip_ok=1.10,
                        t_delay=0.6, t_settle=1.0, t_end=18.0, Dt=0.02, dt_check=0.1,
                        compare=True):
    """Faithful generation-trip overvoltage cascade (Finding #7 of the report).

    Mechanism reproduced:
      1. Solar inverters carry a grid-code **Q-V droop** (``enable_qv_droop``): as
         voltage rises they **absorb** reactive power, like the real VRE (Granada
         at -165 Mvar). At the operating point they are unity-PF (droop = 0).
      2. A seed disturbance (tripping ``n_seed_reactors`` shunt reactors) lifts the
         Iberian voltage -- as the operators' export cut / line coupling did.
      3. **Overvoltage protection** with a **definite-time delay** ``t_delay``: a
         solar unit trips only after its terminal voltage has stayed above the
         setting continuously for ``t_delay`` seconds. Each trip removes that
         unit's reactive **absorption** -> voltage rises -> more units time out ->
         a **slow positive-feedback cascade** developing over seconds (rather than
         the violent zero-delay runaway).

    With ``compare=True`` two protection settings are overlaid:
      * **mis-set** ``V_trip`` (~1.03 pu, below the 1.1-1.2 pu ride-through
        requirement) -> cascade / collapse (the real event);
      * **compliant** ``V_trip_ok`` (1.10 pu) -> the droop rides the disturbance
        through, no trips.

    Tripping is applied between fixed-step integration segments (``dt_check``) via
    ``set_value`` (p and droop -> 0) + ``recalculate_algebraics``.
    """
    region_buses = [b['name'] for b in data['buses']
                    if not b['name'].startswith('GT_') and b['name'] != 'FR']

    def scenario(v_trip):
        model = CasadiModel(build())
        model.Dt = Dt
        model.ini({}, xy_0='xy_0.json')
        solar = enable_qv_droop(model, K_qv=K_qv)
        solar_names = [p['name'] for p in solar]
        solar_bus = {p['name']: p['bus'] for p in solar}

        reactors = []
        for b in data['buses']:
            try:
                val = float(model.get_value(f"b_shunt_{b['name']}"))
            except Exception:
                continue
            if val < -1e-9:
                reactors.append((b['name'], val))
        reactors.sort(key=lambda x: x[1])
        seed = reactors[:max(0, n_seed_reactors)]

        model.run(t_settle, method='trapezoidal')
        for bus, _ in seed:
            model.set_value(f"b_shunt_{bus}", 0.0)
        model.recalculate_algebraics(tol=1e-8)
        print(f"[V_trip={v_trip}] seed {len(seed)} reactors at t={t_settle}s; "
              f"Q-V droop K_qv={K_qv} on {len(solar)} solar")

        tripped, trip_log, collapsed = set(), [], False
        above = {nm: 0.0 for nm in solar_names}   # definite-time timers (s over threshold)
        t = t_settle
        while t < t_end - 1e-9:
            t_next = min(t + dt_check, t_end)
            try:
                model.run(t_next, method='trapezoidal')
            except RuntimeError as e:
                collapsed = True
                print(f"  VOLTAGE COLLAPSE: solver stopped at t~{t_next:.2f}s -- {e}")
                break
            t = t_next
            # accumulate definite-time timers; reset when back below threshold
            newly = []
            for nm in solar_names:
                if nm in tripped:
                    continue
                V = float(model.get_value(f"V_{solar_bus[nm]}"))
                if V > v_trip:
                    above[nm] += dt_check
                    if above[nm] >= t_delay - 1e-9:
                        newly.append((nm, V))
                else:
                    above[nm] = 0.0
            if newly:
                for nm, V in newly:
                    model.set_value(f"p_s_ppc_{nm}", 0.0)   # disconnect: no P
                    model.set_value(f"K_qv_{nm}", 0.0)      # ... and no Q absorption
                    tripped.add(nm); trip_log.append((t, nm, V))
                model.recalculate_algebraics(tol=1e-8)
                print(f"  t={t:5.2f}s  tripped {len(newly)} solar (V>{v_trip}, "
                      f"{t_delay}s delay) -> {len(tripped)}/{len(solar)} off")

        model.post()
        tt = np.array(model.Time)
        Vmat = np.array([model.get_values(f'V_{b}') for b in region_buses])
        vmax = Vmat.max(axis=0)
        if collapsed and len(tt) > 1:      # drop the last non-converged (garbage) step
            tt, vmax = tt[:-1], vmax[:-1]
        n_off = np.array([sum(1 for (tc, _, _) in trip_log if tc <= ti + 1e-9)
                          for ti in tt])
        print(f"  -> {len(tripped)}/{len(solar)} solar tripped, max V {vmax.max():.3f} pu, "
              f"{'COLLAPSE' if collapsed else 'settled'}")
        return tt, vmax, n_off, collapsed, len(solar)

    tt, vmax, n_off, collapsed, n_solar = scenario(V_trip)

    fig, ax = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    if compare:
        to, vo, no, co, _ = scenario(V_trip_ok)
        ax[0].plot(to, vo, 'g-', lw=1.6, label=f'compliant {V_trip_ok} pu (rides through)')
        ax[1].plot(to, no, 'g-', lw=1.6, drawstyle='steps-post',
                   label=f'compliant {V_trip_ok} pu')
    ax[0].plot(tt, vmax, 'r-', lw=2, label=f'mis-set {V_trip} pu (cascade)')
    ax[0].axhline(V_trip, color='orange', ls='--', lw=1)
    ax[0].axhline(1.10, color='green', ls='--', lw=1)
    ax[0].axhline(1.20, color='red', ls=':', lw=1, label='1.20 pu')
    ax[0].axvline(t_settle, color='gray', ls=':', lw=1, label=f'seed @ {t_settle}s')
    ax[0].set_ylabel('max bus voltage [pu]'); ax[0].grid(True); ax[0].legend(loc='best', fontsize=8)
    status = f'CASCADE / collapse at t={tt[-1]:.2f}s' if collapsed else 'arrested'
    ax[0].set_title(f'Generation-trip overvoltage cascade (Q-V droop K_qv={K_qv}, '
                    f'seed {n_seed_reactors} reactors) -> {status}')
    ax[1].plot(tt, n_off, 'r-', lw=1.8, drawstyle='steps-post', label=f'mis-set {V_trip} pu')
    ax[1].set_ylabel('# solar tripped'); ax[1].set_xlabel('Time [s]')
    ax[1].set_ylim(-0.5, n_solar + 0.5); ax[1].grid(True); ax[1].legend(loc='best', fontsize=8)
    fig.tight_layout(); fig.savefig('esp_preblackout_gentrip.png')
    print('saved esp_preblackout_gentrip.png')
    return fig


if __name__ == '__main__':
    ini()       # steady-state initialization + small-signal analysis
    #run_link()  # step the interconnection power and observe the inter-area response
    model = run()         # time-domain run from a saved operating point

    t = model.Time
    p_tie = -np.array(model.get_values('p_g_SLACK_FR'))*50e9/1e6   
    fig, ax = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    ax[0].plot(t, p_tie, color='navy') 
    ax[0].set_ylabel('Spain -> France [MW]')
    ax[0].grid(True)

    for bus in data['buses']:
       ax[1].plot(t, model.get_values(f'V_{bus["name"]}'), color='navy')  
        
    #for syn in data['syns']:
    #    ax[1].plot(t, model.get_values(f'v_f_{syn["name"]}'), color='navy')  

    # omega_coi = model.get_values('omega_coi')
    # for syn in data['syns']:
    #     ax[1].plot(t, model.get_values(f'omega_{syn["name"]}')-omega_coi)  
                
    fig.tight_layout(); fig.savefig('esp_preblackout_link_dev.png')
    print('saved esp_preblackout_link_dev.png')
 