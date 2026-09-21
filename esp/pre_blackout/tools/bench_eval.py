"""Build the bench at a given lambda, tune its stabilisers, report the modes.

    python bench_eval.py <lambda> [--tune]

Without --tune the stabilisers keep the settings already in the file, which is the
fast path while hunting for the right dispatch.  With --tune they are re-tuned by
the residue method on this dispatch, which is what the final bench must ship with.
"""
import os
import subprocess
import sys

import numpy as np
import hjson

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir('/Users/jmmauricio/workspace/benchmarks_public/esp/pre_blackout')

from pydae.bps import BpsBuilder
from pydae.core.builder.casadi_builder import CasadiBuilder
from pydae.core.model.casadi_model import CasadiModel
from pydae.ssa import eig
from pydae.bps.utils.pss_tuner import tune_psss, machines_with_pss, _spec, _param

PY = '/opt/homebrew/Caskroom/miniforge/base/envs/pydae_dev/bin/python'
DATA = 'esp_preblackout_pod_bench.hjson'
LAM = float(sys.argv[1]) if len(sys.argv) > 1 else 0.6
TUNE = '--tune' in sys.argv
ZETA = float(os.environ.get('BENCH_ZETA', '0.10'))

subprocess.run([PY, os.path.join(HERE, 'make_bench.py'), str(LAM)], check=True)

grid = BpsBuilder(DATA, use_casadi=True)
grid.uz_jacs = False
grid.construct('esp_bench')
model = CasadiModel(CasadiBuilder(grid.sys_dict).build())

if TUNE:
    rep = tune_psss(model, DATA, zeta_target=ZETA, xy_0='xy_0.json', verbose=False)
    print(f'\ntuned {len(rep.designs)} stabiliser(s); worst e-mech damping '
          f'{100*rep.floor_before:.2f} -> {100*rep.floor_after:.2f} %')
else:
    model.ini({}, xy_0='xy_0.json')
    eig(model)

data = hjson.load(open(DATA))
S = {s['name']: s['S_n']/1e6 for s in data['syns']}
names = list(S)
print(f"ES-FR export = {-model.get_value('p_g_SLACK_FR')*S['SLACK_FR']:,.0f} MW")


def emech():
    eig(model)
    ev = model.eigenvalues
    f = ev.imag/(2*np.pi)
    z = -ev.real/np.abs(ev)
    rot = np.array([sum(model.participation[model.x_list.index(f'omega_{m}'), k]
                        for m in names) for k in range(len(ev))])
    ks = [k for k in range(len(ev)) if ev[k].imag > 0.05
          and 0.10 <= f[k] <= 2.5 and rot[k] >= 0.15]
    return ev, f, z, ks


ev, f, z, ks = emech()
tgt = [k for k in ks if 0.25 <= f[k] <= 0.55]
ks.sort(key=lambda k: z[k])
print('\nworst electromechanical modes:')
for k in ks[:6]:
    mark = '  <-- target band' if k in tgt else ''
    print(f'   {f[k]:6.3f} Hz   zeta = {100*z[k]:6.2f} %{mark}')

if tgt:
    k = min(tgt, key=lambda k: z[k])
    print(f'\ntarget mode: {f[k]:.3f} Hz  zeta = {100*z[k]:.2f} %   '
          f'(want 1-3 %)   worst overall {100*min(z[k] for k in ks):.2f} %')

# --- sensitivity: is the bench marginal *with* stabilisers tuned? -------------
pss = [(m, t) for m, t in machines_with_pss(DATA)]
gains = {m: model.get_value(_param(_spec(t)['gain'], m)) for m, t in pss}
if any(abs(v) > 0 for v in gains.values()):
    for m, t in pss:
        model.set_value(_param(_spec(t)['gain'], m), 0.0)
    ev2, f2, z2, ks2 = emech()
    t2 = [k for k in ks2 if 0.25 <= f2[k] <= 0.55]
    if t2:
        k2 = min(t2, key=lambda k: z2[k])
        print(f'stabilisers OFF: target {f2[k2]:.3f} Hz  zeta = {100*z2[k2]:.2f} %   '
              f'worst overall {100*min(z2[k] for k in ks2):.2f} %')
    for m, t in pss:
        model.set_value(_param(_spec(t)['gain'], m), gains[m])
