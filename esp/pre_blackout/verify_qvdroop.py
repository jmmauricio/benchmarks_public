"""Verify the pv_pq_ss Q-V droop preserves the pre-blackout equilibrium + SSA."""
import numpy as np
from pydae import ssa
from pydae.bps import BpsBuilder
from pydae.core.builder.casadi_builder import CasadiBuilder
from pydae.core.model.casadi_model import CasadiModel
from pydae.utils import read_data

DATA = 'esp_preblackout_dev.hjson'
data = read_data(DATA)


def build():
    grid = BpsBuilder(DATA, use_casadi=True)
    grid.uz_jacs = False
    grid.construct('esp_preblackout')
    return CasadiBuilder(grid.sys_dict).build()


def enable_qv_droop(model, K_qv=2.0):
    """Seed V_qref to the equilibrium terminal voltage and turn on the droop.

    Because the droop term K_qv*(V_qref - V_m) is zero at V_m = V_qref, the
    algebraic residual at the operating point is unchanged -> no re-init needed.
    """
    solar = [p for p in data['pvs'] if p['type'] == 'pv_pq_ss']
    for p in solar:
        V = float(model.get_value(f"V_{p['bus']}"))
        model.set_value(f"V_qref_{p['name']}", V)
        model.set_value(f"K_qv_{p['name']}", K_qv)
    return solar


model = CasadiModel(build())
model.ini({}, xy_0='xy_0.json')

solar = [p for p in data['pvs'] if p['type'] == 'pv_pq_ss']
q_before = np.array([float(model.get_value(f"q_s_{p['name']}")) for p in solar])
v_before = np.array([float(model.get_value(f"V_{p['bus']}")) for p in solar])
print(f"[base K_qv=0]  n_solar={len(solar)}  sum|q_s|={np.abs(q_before).sum():.4e} pu  "
      f"max|q_s|={np.abs(q_before).max():.4e} pu")

# SSA at base
model.A_eval()
ss0 = ssa.damp(model.A, sort='damp')
z0 = np.asarray(ss0['zetas'])
print(f"[base]  min damping = {100*np.min(z0):.2f} %   n_unstable = {int((z0<0).sum())}")

# Enable droop and re-check equilibrium residual
enable_qv_droop(model, K_qv=2.0)
model.recalculate_algebraics(tol=1e-10)
q_after = np.array([float(model.get_value(f"q_s_{p['name']}")) for p in solar])
v_after = np.array([float(model.get_value(f"V_{p['bus']}")) for p in solar])
print(f"[droop K_qv=2] max|dq_s| = {np.abs(q_after-q_before).max():.4e} pu   "
      f"max|dV| = {np.abs(v_after-v_before).max():.4e} pu")

# SSA with droop on
model.A_eval()
ss1 = ssa.damp(model.A, sort='damp')
z1 = np.asarray(ss1['zetas'])
f1 = np.asarray(ss1['freqs'])
print(f"[droop] min damping = {100*np.min(z1):.2f} %   n_unstable = {int((z1<0).sum())}")

# locate inter-area mode (~0.2 Hz, oscillatory)
osc = f1 > 0.05
near = np.where(osc)[0][np.argsort(np.abs(f1[osc] - 0.2))[:3]]
for i in near:
    print(f"  mode near 0.2 Hz: f={f1[i]:.3f} Hz  zeta={100*z1[i]:.2f} %")
