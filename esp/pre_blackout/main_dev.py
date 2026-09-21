"""
Development / reference harness for the fixed-step implicit-trapezoidal solver,
using CasADi's native linear solver (ca.Linsol) with proper factorization reuse.

Self-contained (does not import main.py). Use it to prototype/benchmark the
linear-solve backend before/against the version integrated in pydae.

    python main_dev.py            # run the interconnection-power test + timing
    python main_dev.py bench      # linear-solver + integrator micro-benchmarks
"""
import sys
import time

import casadi as ca
import hjson
import numpy as np
from matplotlib import pyplot as plt

from pydae.bps import BpsBuilder
from pydae.core.builder.casadi_builder import CasadiBuilder
from pydae.core.model.casadi_model import CasadiModel

DATA = 'esp_preblackout.hjson'


def build():
    grid = BpsBuilder(DATA, use_casadi=True)
    grid.uz_jacs = False
    grid.construct('esp_preblackout')
    return CasadiModel(CasadiBuilder(grid.sys_dict).build())


def warm_ini(model, tol=1e-9):
    """Warm-started load flow (flat start diverges on this 100-bus/130-machine
    model): seed V/theta + rotor angles from the HJSON `results` block and solve
    with a damped Levenberg-Marquardt Newton. Leaves the model at t=0."""
    data = hjson.load(open(DATA))
    model._route_dict({}, phase='ini')
    v = np.concatenate((model.x, model.y_ini)).astype(float)
    p = np.concatenate((model.p_vals, model.u_ini_vals)).astype(float)
    ia = {n: i for i, n in enumerate(list(model.x_names) + list(model.y_ini_names))}

    def sv(name, val):
        if name in ia:
            v[ia[name]] = val

    vth = {b['name']: (b['V_pu'], np.radians(b['theta_deg'])) for b in data['results']['buses']}
    t2r = {t['bus_j']: t['bus_k'] for t in data['transformers']}
    for b, (vv, th) in vth.items():
        sv('V_' + b, vv); sv('theta_' + b, th)
    for term, reg in t2r.items():
        if reg in vth:
            vv, th = vth[reg]; sv('V_' + term, vv); sv('theta_' + term, th)
    for s in data['syns']:
        sv('delta_' + s['name'], vth.get(t2r.get(s.get('bus', ''), ''), (1.0, 0.0))[1])

    resid = lambda x: np.array(model._residual_fn(x, p)).flatten()
    jac = lambda x: np.array(model._jacobian_fn(x, p))
    r = resid(v); lam = 1e-3
    for _ in range(500):
        if np.linalg.norm(r) < tol:
            break
        H = jac(v); A = H.T @ H; g = H.T @ r
        for _ in range(40):
            dv = np.linalg.solve(A + lam * np.diag(np.diag(A) + 1e-10), -g)
            vn = v + dv; rn = resid(vn)
            if np.all(np.isfinite(rn)) and np.linalg.norm(rn) < np.linalg.norm(r):
                v, r = vn, rn; lam = max(lam * 0.7, 1e-10); break
            lam *= 2.5
    model.x, model.y_ini, model.t = v[:model.N_x], v[model.N_x:], 0.0
    model.ini2run()
    return float(np.linalg.norm(r))


class Trapezoidal:
    """Fixed-step implicit-trapezoidal integrator, Jacobian from CasADi AD and the
    linear system solved with **ca.Linsol** using the reuse trio:

        sfact(J)  -- symbolic factorization, ONCE (sparsity is fixed for the sim)
        nfact(J)  -- numeric factorization, only when J is refreshed
        solve(J,b)-- cheap back-substitution, every Newton iteration (reuses nfact)

    Measured on this model (1712x1712, nnz 7308) the 'qr' plugin gives
    nfact ~0.04 ms, solve ~0.014 ms -- ~28x/4x faster than scipy splu, no
    DM<->numpy conversion.
    """

    def __init__(self, model, tol=1e-9, max_iter=10, refactor_every=25,
                 refactor_iters=4, linsol='qr'):
        self.m = model
        self.tol, self.max_iter = tol, max_iter
        self.refactor_every, self.refactor_iters = refactor_every, refactor_iters
        d = model.dae_dict
        x, z, p, f, g = d['x'], d['z'], d['p'], d['ode'], d['alg']
        self.nx = x.numel()
        w = ca.vertcat(x, z)
        Dt, xn, fn = ca.SX.sym('Dt'), ca.SX.sym('xn', self.nx), ca.SX.sym('fn', self.nx)
        R = ca.vertcat(x - xn - (Dt / 2) * (fn + f), g)
        self._f = ca.Function('f', [x, z, p], [f])
        self._R = ca.Function('R', [w, xn, fn, p, Dt], [R])
        self._J = ca.Function('J', [w, xn, fn, p, Dt], [ca.jacobian(R, w)])
        self._lin = ca.Linsol('trap_lin', linsol, self._J.sparsity_out(0))
        self._sfact_done = False
        self._need_nfact = True
        self._Jdm = None
        self._sc = 0

    def refactor(self):
        """Force a numeric refactorization on the next step (after a set_value)."""
        self._need_nfact = True

    def _newton(self, w, xn, fn, p, Dt):
        for it in range(self.max_iter):
            Rdm = self._R(w, xn, fn, p, Dt)
            if float(ca.norm_inf(Rdm)) < self.tol:
                return it
            if self._need_nfact or (self._sc % self.refactor_every == 0 and it == 0):
                self._Jdm = self._J(w, xn, fn, p, Dt)
                if not self._sfact_done:
                    self._lin.sfact(self._Jdm)          # symbolic factorization ONCE
                    self._sfact_done = True
                self._lin.nfact(self._Jdm)              # numeric factorization
                self._need_nfact = False
            dw = self._lin.solve(self._Jdm, -Rdm)       # reuse factorization
            w += np.array(dw).flatten()                 # in-place: caller sees the update
        return None

    def run(self, t_end):
        m = self.m
        m._route_dict({}, phase='run')
        Dt = m.Dt
        p = np.concatenate((m.p_vals, m.u_run_vals))
        x = np.asarray(m.x, float).flatten()
        z = np.asarray(m.y_run, float).flatten()
        n = int(round((t_end - m.t) / Dt))
        has_h = m._h_fn is not None and 'h_dict' in m.sys_dict
        for _ in range(n):
            fn = np.array(self._f(x, z, p)).flatten()
            w = np.concatenate((x, z))
            it = self._newton(w, x, fn, p, Dt)
            if it is None:                              # retry once with fresh factor
                self._need_nfact = True
                w = np.concatenate((x, z))
                it = self._newton(w, x, fn, p, Dt)
                if it is None:
                    raise RuntimeError(f"trapezoidal did not converge at t={m.t + Dt:.4f} s")
            elif it > self.refactor_iters:              # hard step -> refresh J next step
                self._need_nfact = True
            x, z = w[:self.nx], w[self.nx:]             # _newton updated w in place
            m.t += Dt
            self._sc += 1
            m.Time.append(m.t); m.X.append(x.copy()); m.Y.append(z.copy())
            if has_h:
                m.Z.append(np.array(m._h_fn(x, z, p)).flatten())
        m.x, m.y_run = x, z


def run_link(dP_mw=-1000.0, machine=None, t_step=1.0, t_end=25.0, Dt=0.05, linsol='qr'):
    """Step the Spain->France interconnection power and observe the inter-area
    response, integrated with the ca.Linsol trapezoidal solver."""
    data = hjson.load(open(DATA))
    spain = [s for s in data['syns'] if s['name'] != 'SLACK_FR']
    if machine is None:
        machine = max(spain, key=lambda s: s['S_n'])['name']
    S_n = next(s['S_n'] for s in data['syns'] if s['name'] == machine)
    S_FR = next(s['S_n'] for s in data['syns'] if s['name'] == 'SLACK_FR')

    model = build()
    model.Dt = Dt
    print(f"warm-start load flow residual = {warm_ini(model):.1e}")
    model.Time, model.X, model.Y, model.Z = [], [], [], []

    solver = Trapezoidal(model, refactor_every=10, linsol=linsol)
    t0 = time.perf_counter()
    solver.run(t_step)

    ref = f"p_c_lc_{machine}"
    model.set_value(ref, model.get_value(ref) + dP_mw * 1e6 / S_n)
    model.recalculate_algebraics(tol=1e-8)
    solver.refactor()
    print(f"t={t_step}s: step {machine} by {dP_mw:+.0f} MW -> interconnection power change")

    solver.run(t_end)
    model.post()
    print(f'run_link ({linsol}) complete in {time.perf_counter() - t0:.2f} s '
          f'({len(model.Time)} points over {t_end:.0f} s)')

    t = model.Time
    p_tie = -np.array(model.get_values('p_g_SLACK_FR')) * S_FR / 1e6
    f_es = np.mean([model.get_values(f"omega_{s['name']}") for s in spain], axis=0)
    f_fr = np.array(model.get_values('omega_SLACK_FR'))
    fig, ax = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    ax[0].plot(t, p_tie, color='navy'); ax[0].set_ylabel('Spain -> France [MW]'); ax[0].grid(True)
    ax[0].set_title(f'Interconnection power step {dP_mw:+.0f} MW at {machine}  (ca.Linsol[{linsol}])')
    ax[1].plot(t, (f_es - 1) * 50e3, 'crimson', label='Spain (COI)')
    ax[1].plot(t, (f_fr - 1) * 50e3, 'navy', label='Europe (SLACK_FR)')
    ax[1].set_ylabel('Freq. deviation [mHz]'); ax[1].set_xlabel('Time (s)'); ax[1].legend(); ax[1].grid(True)
    fig.tight_layout(); fig.savefig('esp_preblackout_link_dev.png')
    print('saved esp_preblackout_link_dev.png')
    return model


def bench():
    """Validate + time the ca.Linsol trapezoidal vs an IDAS grid reference over a
    3 s transient after a -1000 MW machine step."""
    data = hjson.load(open(DATA))
    mac = max((s for s in data['syns'] if s['name'] != 'SLACK_FR'), key=lambda s: s['S_n'])
    ref = f"p_c_lc_{mac['name']}"; dP = -1000e6 / mac['S_n']

    def prep(Dt=0.02):
        m = build(); m.Dt = Dt; warm_ini(m)
        m.set_value(ref, m.get_value(ref) + dP); m.recalculate_algebraics(tol=1e-10)
        m.Time, m.X, m.Y, m.Z = [], [], [], []
        return m

    # IDAS grid reference
    mr = prep(); mr._route_dict({}, phase='run')
    pr = np.concatenate((mr.p_vals, mr.u_run_vals)); gr = mr.Dt * np.arange(1, int(3.0 / mr.Dt) + 1)
    ig = ca.integrator('i', 'idas', mr.dae_dict, 0.0, gr,
                       {'reltol': 1e-9, 'abstol': 1e-9, 'linear_solver': 'csparse', 'print_stats': False})
    t = time.perf_counter(); rr = ig(x0=mr.x, z0=mr.y_run, p=pr); tid = time.perf_counter() - t
    xref = np.array(rr['xf'])[:, -1].flatten()
    print(f"IDAS grid 3 s: {tid:.2f} s")
    for ls in ('qr', 'csparse', 'lapacklu'):
        m = prep(); solver = Trapezoidal(m, refactor_every=10, linsol=ls)
        t = time.perf_counter(); solver.run(3.0); dt = time.perf_counter() - t
        print(f"trap ca.Linsol[{ls:9}]: {dt:.2f} s ({dt / (3.0 / 0.02) * 1000:.2f} ms/step)  "
              f"max|x-IDAS|={np.max(np.abs(m.x - xref)):.1e}")


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'bench':
        bench()
    else:
        run_link()
