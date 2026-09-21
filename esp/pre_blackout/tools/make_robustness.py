"""Robustness family for POD design, plus the low-irradiance case.

    python tools/make_robustness.py

Cases 4a/4b/4c hold the dispatch of case 2 fixed and vary only the strength of the
ES-FR interconnection, placing the inter-area mode across the band a POD controller
has to cover.  This is the X_L sweep NTS-SEPE prescribes: a controller phase-
compensated at one frequency has no demonstrated validity at another, and the
lead-lag networks in this system already carry 24-36 degrees of residual phase
error because two stages cannot cover a wide band at once.

Case 5 is the other axis: what is left when the resource is not there.

Stabilisers are re-tuned on every variant, to the same 5% criterion.  Re-using a
tuning from a different network would make any margin an artefact of stale tuning
rather than a property of the case -- the failure mode this whole benchmark set is
built to avoid.
"""
import os
import subprocess
import sys

import numpy as np
import hjson

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
os.chdir(BASE)

from pydae.bps import BpsBuilder
from pydae.core.builder.casadi_builder import CasadiBuilder
from pydae.core.model.casadi_model import CasadiModel
from pydae.ssa import eig
from pydae.bps.utils.pss_tuner import tune_psss

PY = sys.executable
TMP = 'esp_robustness_tmp.hjson'

# (file, label, lambda, N-1, tie X in ohms or None, export trim MW, description)
VARIANTS = [
    ('esp_case4a_tie_strong.hjson', '4a', 0.52, False, None, 1000.0,
     'both ES-FR circuits in service'),
    ('esp_case4b_tie_base.hjson',   '4b', 0.52, True,  None, 1000.0,
     'one circuit (N-1) at nominal reactance -- identical to case 2'),
    ('esp_case4c_tie_weak.hjson',   '4c', 0.52, True,  73.0, 1000.0,
     'one circuit, reactance raised to 73 ohm'),
    ('esp_case5_low_irradiance.hjson', '5', 1.20, True, None, 0.0,
     'low irradiance: PV backed off, synchronous commitment raised to cover it'),
]


def build(path, lam, n1, tie_x, trim):
    env = dict(os.environ, BENCH_OUT=os.path.join(BASE, TMP),
               BENCH_N1='1' if n1 else '0',
               BENCH_EXPORT_TRIM=str(trim))
    if tie_x is not None:
        env['BENCH_TIE_X'] = str(tie_x)
    else:
        env.pop('BENCH_TIE_X', None)
    out = subprocess.run([PY, os.path.join(HERE, 'make_bench.py'), str(lam)],
                         check=True, env=env, capture_output=True, text=True)

    grid = BpsBuilder(TMP, use_casadi=True)
    grid.uz_jacs = False
    grid.construct('esp_rob')
    model = CasadiModel(CasadiBuilder(grid.sys_dict).build())
    rep = tune_psss(model, TMP, zeta_target=0.05, xy_0='xy_0.json',
                    write=TMP, verbose=False)
    return model, rep, out.stdout


def modes(model, data):
    eig(model)
    names = [s['name'] for s in data['syns']]
    ev = model.eigenvalues
    f, z = ev.imag/(2*np.pi), -ev.real/np.abs(ev)
    rot = np.array([sum(model.participation[model.x_list.index(f'omega_{m}'), k]
                        for m in names) for k in range(len(ev))])
    ks = [k for k in range(len(ev)) if ev[k].imag > 0.05
          and 0.10 <= f[k] <= 2.5 and rot[k] >= 0.15]
    # Wide band: weakening the tie pushes the inter-area mode well below 0.2 Hz,
    # and raising synchronous commitment pushes it down too.  A narrow window
    # silently reports "no mode" when the mode has simply moved.
    tgt = min((k for k in ks if 0.08 <= f[k] <= 0.60), key=lambda k: z[k],
              default=None)
    worst = min(ks, key=lambda k: z[k])
    return f, z, tgt, worst, sum(1 for k in ks if z[k] < 0)


only = sys.argv[1:]
rows = []
for path, label, lam, n1, tie_x, trim, desc in VARIANTS:
    if only and label not in only:
        continue
    model, rep, gen = build(path, lam, n1, tie_x, trim)
    data = hjson.load(open(TMP))
    f, z, tgt, worst, n_uns = modes(model, data)

    S = {s['name']: s['S_n']/1e6 for s in data['syns']}
    syn = [s for s in data['syns'] if s['name'] != 'SLACK_FR']
    ke = sum(s['S_n']/1e6*s['H'] for s in syn)
    pv = sum(p['S_n']/1e6*p.get('p_s_ppc', 0.0) for p in data['pvs'])
    exp_mw = -model.get_value('p_g_SLACK_FR')*S['SLACK_FR']
    f_t = f[tgt] if tgt is not None else float('nan')
    z_t = 100*z[tgt] if tgt is not None else float('nan')

    hdr = [
        f'CASE {label} - ' + ('robustness family: ' if label.startswith('4')
                              else '') + desc,
        '',
    ]
    if label.startswith('4'):
        hdr += [
            'One of three variants that hold the dispatch of case 2 fixed and vary',
            'only the interconnection strength, so a POD controller can be designed',
            'and checked against a locus of inter-area frequencies rather than a',
            'single point.  This is the X_L sweep NTS-SEPE prescribes for robustness.',
            '',
        ]
    else:
        hdr += [
            'The other axis: PV output is low, so a PV-based damping controller has',
            'little authority.  Documents where the approach stops being applicable.',
            'Note the system is also less stressed here - fewer converters and more',
            'synchronous plant - so the reduced authority is not as costly as it',
            'looks.',
            '',
        ]
    hdr += [
        f'Expected:  inter-area mode {f_t:.3f} Hz at zeta = {z_t:.2f} %',
        f'           worst mode {f[worst]:.3f} Hz at zeta = {100*z[worst]:.2f} %',
        f'           unstable modes: {n_uns}',
        '',
        f'Synchronous commitment {100*lam:.0f} % of base ({ke:,.0f} MW.s), '
        f'PV {pv:,.0f} MW,',
        f'ES-FR {"one circuit (N-1)" if n1 else "both circuits"}'
        + (f' at {tie_x:.0f} ohm' if tie_x else '') + f', export {exp_mw:,.0f} MW.',
        f'Stabilisers tuned by the residue method to the 5 % criterion '
        f'({len(rep.designs)} committed).',
        '',
        'Generated by tools/make_robustness.py.  See CASES.md.',
    ]
    body = open(TMP, encoding='utf-8').read()
    open(path, 'w', encoding='utf-8').write(
        ''.join(f'// {ln}\n' if ln else '//\n' for ln in hdr) + body)

    rows.append((label, path, f_t, z_t, f[worst], 100*z[worst], n_uns, exp_mw, pv))
    print(f'{label:3} {path:34} inter-area {f_t:.3f} Hz {z_t:6.2f} %   '
          f'worst {f[worst]:.3f} Hz {100*z[worst]:6.2f} %   unstable {n_uns}')

if os.path.exists(TMP):
    os.remove(TMP)

print('\nmodal locus for POD design:')
for label, _p, f_t, z_t, *_r in rows:
    if label.startswith('4'):
        print(f'   {label}   {f_t:.3f} Hz   zeta {z_t:5.2f} %')
