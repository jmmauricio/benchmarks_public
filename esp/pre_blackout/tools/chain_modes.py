"""Slow modes of the CE-chain case, with the four-group mode shape."""
import os
import sys
import numpy as np
import hjson

os.chdir('/Users/jmmauricio/workspace/benchmarks_public/esp/pre_blackout')

from pydae.bps import BpsBuilder
from pydae.core.builder.casadi_builder import CasadiBuilder
from pydae.core.model.casadi_model import CasadiModel
from pydae.ssa import eig

DATA = 'esp_preblackout_ce_chain.hjson'
TARGET_MW = 1400.0

data = hjson.load(open(DATA))
S = {s['name']: s['S_n']/1e6 for s in data['syns']}
EXT = ['SLACK_FR', 'CE_CENTRE', 'CE_EAST']
sched = {s['name']: s.get('lc', {}).get('p_c_lc_mw', 0.0)
         for s in data['syns'] if s['name'] not in EXT}
total_sched = sum(sched.values())

grid = BpsBuilder(DATA, use_casadi=True)
grid.uz_jacs = False
grid.construct('esp_chain')
model = CasadiModel(CasadiBuilder(grid.sys_dict).build())

trim = float(sys.argv[1]) if len(sys.argv) > 1 else 625.0
u = {}
scale = (total_sched - trim)/total_sched
for name, p in sched.items():
    if p:
        u[f'p_c_lc_{name}'] = p*scale*1e6/(S[name]*1e6)
model.ini(u, xy_0='xy_0.json')
eig(model)

print(f"ES-FR export = {-model.get_value('p_g_SLACK_FR')*S['SLACK_FR']:,.0f} MW")

ev, V = model.eigenvalues, model.right_eigenvectors
f = ev.imag/(2*np.pi)
z = -ev.real/np.abs(ev)
all_syn = [s['name'] for s in data['syns']]
rotor = np.array([sum(model.participation[model.x_list.index(f'omega_{m}'), k]
                      for m in all_syn) for k in range(len(ev))])

# Rotor participation separates electromechanical modes from the AVR / governor /
# load-controller modes, of which this model has many sitting around 0.1 Hz with
# ~45% damping and no rotor content at all.
slow = [k for k in range(len(ev)) if ev[k].imag > 0.05 and 0.05 <= f[k] <= 0.90
        and rotor[k] >= 0.15]
slow.sort(key=lambda k: f[k])

iber = [s['name'] for s in data['syns'] if s['name'] not in EXT]
print(f'\nslow modes 0.10-0.60 Hz.  Each group: |d.delta| (normalised to the '
      f'largest) and its angle in degrees.')
print(f"{'f [Hz]':>7} {'zeta%':>7} | {'Iberia':>15} {'West':>15} {'Centre':>15} "
      f"{'East':>15} | {'part: Ib':>8} {'W':>6} {'C':>6} {'E':>6}")


def group_vec(k, names):
    """Inertia-weighted mean of the group's rotor-angle components."""
    num = sum(V[model.x_list.index(f'delta_{m}'), k]*S[m]*data_H[m] for m in names)
    den = sum(S[m]*data_H[m] for m in names)
    return num/den


def group_part(k, names):
    return sum(model.participation[model.x_list.index(f'omega_{m}'), k] for m in names)


data_H = {s['name']: s['H'] for s in data['syns']}
groups = [('Iberia', iber), ('West', ['SLACK_FR']),
          ('Centre', ['CE_CENTRE']), ('East', ['CE_EAST'])]

for k in slow:
    vecs = [group_vec(k, names) for _lbl, names in groups]
    scale = max(abs(v) for v in vecs)
    ref = vecs[3]                                  # phases relative to East
    cells = ''.join(f'{abs(v)/scale:7.3f}@{np.angle(v/ref, deg=True):7.1f}'
                    for v in vecs)
    parts = ''.join(f'{group_part(k, names):7.3f}' for _lbl, names in groups)
    print(f'{f[k]:7.3f} {100*z[k]:7.2f} | {cells} |{parts}')
