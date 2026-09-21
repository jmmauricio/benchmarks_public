"""Check every case file reproduces the modes its own header claims.

    python tools/verify_cases.py [file ...]

The expected values are parsed out of each file's header comment rather than kept
in a table here, so the check cannot drift out of step with what the files say
about themselves.  A fresh model per case: re-initialising one instance twice
converges somewhere else.
"""
import glob
import os
import re
import sys

import numpy as np
import hjson

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(BASE)

from pydae.bps import BpsBuilder
from pydae.core.builder.casadi_builder import CasadiBuilder
from pydae.core.model.casadi_model import CasadiModel
from pydae.ssa import eig

PAT = re.compile(r'([\d.]+)\s*Hz\s*at\s*zeta\s*=\s*(-?[\d.]+)\s*%')


def claimed(path):
    """(f, zeta) the header claims for the inter-area mode, or None."""
    head = [ln for ln in open(path, encoding='utf-8').read().split('\n')
            if ln.strip().startswith('//')][:40]
    # Prefer a line that names the inter-area mode *and* carries a value, then fall
    # back to the Expected block.  Not `or`: a header can mention "inter-area" in
    # prose with no number, which would otherwise mask the value further down.
    inter = [ln for ln in head if 'inter-area' in ln.lower()]
    for ln in inter + [ln for ln in head if 'Expected' in ln]:
        m = PAT.search(ln)
        if m:
            return float(m.group(1)), float(m.group(2))
    return None


files = sys.argv[1:] or sorted(glob.glob('esp_case*.hjson'))
bad = 0

for path in files:
    want = claimed(path)
    grid = BpsBuilder(path, use_casadi=True)
    grid.uz_jacs = False
    grid.construct('esp_verify')
    model = CasadiModel(CasadiBuilder(grid.sys_dict).build())
    model.ini({}, xy_0='xy_0.json')
    eig(model)

    data = hjson.load(open(path))
    names = [s['name'] for s in data['syns']]
    S_FR = next(s['S_n'] for s in data['syns'] if s['name'] == 'SLACK_FR')/1e6
    ev = model.eigenvalues
    f, z = ev.imag/(2*np.pi), -ev.real/np.abs(ev)
    rot = np.array([sum(model.participation[model.x_list.index(f'omega_{m}'), k]
                        for m in names) for k in range(len(ev))])
    ks = [k for k in range(len(ev)) if ev[k].imag > 0.05
          and 0.10 <= f[k] <= 2.5 and rot[k] >= 0.15]
    worst = min(ks, key=lambda k: z[k])
    n_uns = sum(1 for k in ks if z[k] < 0)

    print(f'\n{path}')
    print(f'   export {-model.get_value("p_g_SLACK_FR")*S_FR:8,.0f} MW   '
          f'{model.N_x} states   worst {f[worst]:.3f} Hz {100*z[worst]:6.2f} %   '
          f'unstable {n_uns}')
    if want is None:
        print('   header states no expected mode')
        continue
    f_w, z_w = want
    near = [k for k in ks if abs(f[k] - f_w) < 0.05]
    if not near:
        print(f'   MISMATCH: no mode near {f_w:.3f} Hz')
        bad += 1
        continue
    k = min(near, key=lambda k: abs(f[k] - f_w))
    ok = abs(100*z[k] - z_w) < 0.5
    bad += 0 if ok else 1
    print(f'   claimed {f_w:.3f} Hz {z_w:6.2f} %  ->  found {f[k]:.3f} Hz '
          f'{100*z[k]:6.2f} %   {"OK" if ok else "MISMATCH"}')

print(f'\n{len(files)} case(s), {bad} mismatch(es)')
sys.exit(1 if bad else 0)
