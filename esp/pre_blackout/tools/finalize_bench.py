"""Write the final benchmark: dispatch + N-1 + stabilisers tuned into the file.

    BENCH_N1=1 BENCH_EXPORT_TRIM=1000 BENCH_ZETA=0.05 python finalize_bench.py <lambda>

Unlike bench_eval.py this passes write= to the tuner, so the stabiliser designs end
up in the hjson.  Without that the file does not reproduce the reported numbers --
the tuning would live only in the process that computed them.
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
from pydae.bps.utils.pss_tuner import (tune_psss, machines_with_pss, loop_residue,
                                       _spec, _param)

PY = '/opt/homebrew/Caskroom/miniforge/base/envs/pydae_dev/bin/python'
DATA = 'esp_preblackout_pod_bench.hjson'
LAM = float(sys.argv[1]) if len(sys.argv) > 1 else 0.56
ZETA = float(os.environ.get('BENCH_ZETA', '0.05'))

subprocess.run([PY, os.path.join(HERE, 'make_bench.py'), str(LAM)], check=True)

grid = BpsBuilder(DATA, use_casadi=True)
grid.uz_jacs = False
grid.construct('esp_bench_final')
model = CasadiModel(CasadiBuilder(grid.sys_dict).build())

rep = tune_psss(model, DATA, zeta_target=ZETA, xy_0='xy_0.json',
                write=DATA, verbose=False)
print(f'\ntuned {len(rep.designs)} stabiliser(s), written into {DATA}')
for d in rep.designs:
    print(f'   {d.machine:20} K {d.gain:9.4f}  T_1 {d.T_1:.4f}  T_2 {d.T_2:.4f}'
          f'   {d.f_hz:.3f} Hz  {100*d.zeta_before:.2f} -> {100*d.zeta_after:.2f} %')

data = hjson.load(open(DATA))
S = {s['name']: s['S_n']/1e6 for s in data['syns']}
names = list(S)
iber = [n for n in names if n != 'SLACK_FR']
ke = sum(s['S_n']/1e6*s['H'] for s in data['syns'] if s['name'] != 'SLACK_FR')


def modes():
    eig(model)
    ev = model.eigenvalues
    f, z = ev.imag/(2*np.pi), -ev.real/np.abs(ev)
    rot = np.array([sum(model.participation[model.x_list.index(f'omega_{m}'), k]
                        for m in names) for k in range(len(ev))])
    ks = [k for k in range(len(ev)) if ev[k].imag > 0.05
          and 0.10 <= f[k] <= 2.5 and rot[k] >= 0.15]
    return ev, f, z, ks


ev, f, z, ks = modes()
tgt = min((k for k in ks if 0.15 <= f[k] <= 0.55), key=lambda k: z[k], default=None)
print(f'\nES-FR export {-model.get_value("p_g_SLACK_FR")*S["SLACK_FR"]:,.0f} MW'
      f'   Iberian stored energy {ke:,.0f} MW.s ({100*ke/119474:.0f} % of 12:30)')
print(f'worst electromechanical mode : {f[min(ks, key=lambda k: z[k])]:.3f} Hz  '
      f'{100*min(z[k] for k in ks):.2f} %')
if tgt is not None:
    print(f'target inter-area mode       : {f[tgt]:.3f} Hz  {100*z[tgt]:.2f} %')

# --- sensitivity: stabilisers as shipped vs off ------------------------------
pss = machines_with_pss(DATA)
gains = {m: model.get_value(_param(_spec(t)['gain'], m)) for m, t in pss}
for m, t in pss:
    model.set_value(_param(_spec(t)['gain'], m), 0.0)
ev2, f2, z2, ks2 = modes()
t2 = min((k for k in ks2 if 0.15 <= f2[k] <= 0.55), key=lambda k: z2[k], default=None)
print(f'\nsensitivity, stabilisers OFF : worst {100*min(z2[k] for k in ks2):.2f} %'
      + (f'   target mode {f2[t2]:.3f} Hz {100*z2[t2]:.2f} %' if t2 is not None else ''))
for m, t in pss:
    model.set_value(_param(_spec(t)['gain'], m), gains[m])

# --- controllability of the target mode from synchronous plant ---------------
eig(model)
if tgt is not None:
    R = {m: abs(loop_residue(model, m, tgt)) for m, _t in pss}
    top = sorted(R, key=lambda m: -R[m])[:3]
    print('largest synchronous residues on the target mode: '
          + ', '.join(f'{m} {R[m]:.4f}' for m in top))
