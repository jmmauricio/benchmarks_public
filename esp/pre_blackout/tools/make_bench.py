"""Generate the POD test bench: a deeper converter-penetration dispatch.

    python make_bench.py <lambda>      lambda = committed synchronous capacity,
                                       as a fraction of the base case

The synchronous machines in this case are plant aggregates, so scaling a machine's
S_n *and* its scheduled MW together is exactly "fewer identical units of that plant
committed": the per-unit machine data, and therefore the physics of each unit, are
untouched.  PV takes up the released energy.

Two things are deliberately NOT used to make this system marginal: stabiliser gains
and stabiliser tuning.  Those stay correct -- they are re-tuned by the residue
method afterwards, on this dispatch.  Marginality has to come from the dispatch
alone, or the benchmark cannot be published without it reading as a claim about
somebody's tuning practice.

Text insertion, not an hjson round-trip: dumping reformats the file and discards
the commented-out pss blocks.
"""
import os
import re
import sys

import hjson

BASE = '/Users/jmmauricio/workspace/benchmarks_public/esp/pre_blackout'
SRC = os.environ.get('BENCH_SRC', f'{BASE}/esp_preblackout_pss_tuned.hjson')
OUT = os.environ.get('BENCH_OUT', f'{BASE}/esp_preblackout_pod_bench.hjson')

LAM = float(sys.argv[1]) if len(sys.argv) > 1 else 0.6

# The external equivalent is an area aggregate, not a machine.  Its rotor fluxes
# are held effectively constant over the period of a sub-Hz mode -- the classical
# form `gencls` exists for.  Leaving it with SLACK_FR's turbo data is not neutral:
# T1q0 = 0.8 s puts the q-axis corner at 1.25 rad/s, right on top of these modes,
# and that one parameter was worth 12 points of damping.
T_FLUX, T_SUB, D_EXT = 100.0, 10.0, 0.0

text = open(SRC, encoding='utf-8').read()
data = hjson.loads(text)
syns = [s for s in data['syns'] if s['name'] != 'SLACK_FR']


def set_field(txt, anchor_name, key, value, fmt='{}'):
    """Replace one numeric field inside the entry named *anchor_name*."""
    a = txt.index(f'"name": "{anchor_name}"')
    m = re.compile(rf'("{key}":\s*)(-?[\d.eE+]+)').search(txt, a)
    if m is None:
        raise ValueError(f'{key} not found for {anchor_name}')
    return txt[:m.start(2)] + fmt.format(value) + txt[m.end(2):]


def comment_out_line(txt, bus_j, bus_k):
    """Comment out one line entry -- an N-1 outage of that circuit.

    The ES-FR interconnection is two circuits; opening one roughly doubles the tie
    reactance and acts directly on the Iberia-versus-external mode.  A credible
    contingency every operator plans for, which is why it can carry the stress in
    this benchmark without the case implying anything about anyone's practice.
    """
    key = re.search(rf'"bus_j":\s*"{bus_j}",\s*\n\s*"bus_k":\s*"{bus_k}"', txt)
    if key is None:
        raise ValueError(f'line {bus_j}-{bus_k} not found')
    start = txt.rfind('{', 0, key.start())
    start = txt.rfind('\n', 0, start) + 1
    depth, i = 0, txt.index('{', start)
    while True:
        if txt[i] == '{':
            depth += 1
        elif txt[i] == '}':
            depth -= 1
            if depth == 0:
                break
        i += 1
    end = i + 1
    tail = re.match(r'\s*,', txt[end:])
    if tail:
        end += tail.end()
    body = txt[start:end]
    indent = re.match(r'[ \t]*', body).group(0)
    out = '\n'.join(f'{indent}// {ln[len(indent):]}' if ln.startswith(indent) and ln.strip()
                    else (f'// {ln}' if ln.strip() else ln)
                    for ln in body.split('\n'))
    return txt[:start] + out + txt[end:]


# --- N-1 on the interconnection ----------------------------------------------
if os.environ.get('BENCH_N1', '0') == '1':
    text = comment_out_line(text, 'ES213', 'FR')

# --- interconnection strength -------------------------------------------------
# BENCH_TIE_X overrides the reactance of the ES512-FR circuit, in ohms.  This is
# the X_L sweep NTS-SEPE prescribes for robustness assessment: raising it lowers
# the frequency of the slow inter-area mode, lowering it raises the frequency, so
# a controller can be designed and checked against a locus rather than a point.
TIE_X = os.environ.get('BENCH_TIE_X')
if TIE_X:
    m = re.search(r'"bus_j":\s*"ES512",\s*\n\s*"bus_k":\s*"FR".*?"X":\s*([\d.eE+]+)',
                  text, re.S)
    if m is None:
        raise ValueError('ES512-FR circuit not found')
    text = text[:m.start(1)] + f'{float(TIE_X):.5f}' + text[m.end(1):]

# --- external equivalent -----------------------------------------------------
for key, val in (('T1d0', T_FLUX), ('T1q0', T_FLUX),
                 ('T2d0', T_SUB), ('T2q0', T_SUB), ('D', D_EXT)):
    text = set_field(text, 'SLACK_FR', key, val)

# --- decommit synchronous capacity -------------------------------------------
# Merit order, not a uniform scaling.  Nuclear is must-run and the large CCGTs are
# the thermal backbone; what comes off at midday is mid-merit hydro, biomass and
# run-of-river.  Scaling every machine by the same factor instead shrinks the very
# units that host the effective stabilisers at the same rate as everything else, so
# inertia and stabiliser authority collapse together -- which produced a ~30-point
# swing in damping across a 10-point change in committed capacity, far too sharp an
# edge for a benchmark anyone else has to reproduce.
# Off by default, and the arithmetic says why: nuclear and the large CCGTs carry
# 67,082 of the 106,754 MW.s in this case -- 63% of the inertia.  Holding all of
# them at full commitment therefore floors stored energy at 63% of base however
# much mid-merit plant is decommitted, which is close to the uniform-60% case that
# came out comfortably damped.  Merit order alone cannot reach the marginal regime;
# getting there means taking some baseload off too.
BASELOAD = os.environ.get('BENCH_BASELOAD', '0') == '1'


def committed(name):
    if BASELOAD and (name.startswith('nuclear') or name.startswith('CCGT')):
        return 1.0
    return LAM


released = 0.0
for s in syns:
    lam = committed(s['name'])
    text = set_field(text, s['name'], 'S_n', s['S_n']*lam, '{:.1f}')
    p = s.get('lc', {}).get('p_c_lc_mw')
    if p:
        text = set_field(text, s['name'], 'p_c_lc_mw', p*lam, '{:.4f}')
        released += p*(1.0 - lam)

# --- PV takes up the slack ---------------------------------------------------
# With a circuit out, the remaining one has P_max = 400^2/60.4 = 2,647 MW, so a
# ~1,950 MW exchange sits at 74% of it and roughly 48 degrees across the tie --
# past where the load flow will solve.  Reducing the exchange is what an operator
# does when an interconnector trips, so the trim is part of the contingency
# scenario rather than a numerical convenience.
TRIM = float(os.environ.get('BENCH_EXPORT_TRIM', '0'))

pv_mw = sum(p['S_n']/1e6*p.get('p_s_ppc', 0.0) for p in data['pvs'])
gain = (pv_mw + released - TRIM)/pv_mw
worst = 0.0
for p in data['pvs']:
    new = p.get('p_s_ppc', 0.0)*gain
    worst = max(worst, new)
    text = set_field(text, p['name'], 'p_s_ppc', new, '{:.6f}')

open(OUT, 'w', encoding='utf-8').write(text)

ke_base = sum(s['S_n']/1e6*s['H'] for s in syns)
ke0 = sum(s['S_n']/1e6*s['H']*committed(s['name']) for s in syns)
syn_mw = sum(s.get('lc', {}).get('p_c_lc_mw', 0.0)*committed(s['name']) for s in syns)
w_mw = sum(w['S_n']/1e6*w.get('p_s_ppc', 0.0) for w in data['wecs'])
tot = syn_mw + pv_mw + released + w_mw
kept = [s['name'] for s in syns if committed(s['name']) == 1.0]
print(f'wrote {OUT}')
print(f'  rule: {"baseload held, mid-merit at" if BASELOAD else "uniform scaling at"} '
      f'{100*LAM:.0f} %   ({len(kept)} unit(s) held at full commitment)')
print(f'  committed synchronous capacity : {100*ke0/ke_base:.0f} % of base by energy')
print(f'  synchronous  {syn_mw:8,.0f} MW   ({100*syn_mw/tot:4.1f} % of generation)')
print(f'  PV           {pv_mw + released:8,.0f} MW   (+{released:,.0f} MW released, '
      f'max unit {worst:.3f} pu)')
print(f'  wind         {w_mw:8,.0f} MW')
print(f'  Iberian stored energy {ke0:9,.0f} MW.s   '
      f'({100*ke0/119474:.0f} % of the 119,474 MW.s reported at 12:30)')
