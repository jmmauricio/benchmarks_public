"""Build the CE-chain case from the retuned one by text insertion.

Text insertion rather than hjson round-trip: dumping would reformat the whole file
and throw away the commented-out pss blocks.

    Iberia ==(existing ES-FR ties)== FR/West ==X_wc== Centre ==X_ce== East

SLACK_FR keeps its name (the agc block and the p_g_SLACK_FR diagnostics refer to
it) and its role as angle reference, but shrinks from representing all of
Continental Europe to representing the West group alone.  Its 2,000,000 MW.s of
stored energy is redistributed across the three groups, so total CE inertia is
unchanged -- only its *distribution* is, which is the whole point.
"""
import os
import re
import sys

SRC = os.environ.get(
    'CHAIN_SRC',
    '/Users/jmmauricio/workspace/benchmarks_public/esp/pre_blackout/'
    'esp_preblackout_pss_tuned.hjson')
OUT = os.environ.get(
    'CHAIN_OUT',
    '/Users/jmmauricio/workspace/benchmarks_public/esp/pre_blackout/'
    'esp_preblackout_ce_chain.hjson')

# --- chain parameters (the ones the sizing step moves) -----------------------
X_WC = float(sys.argv[1]) if len(sys.argv) > 1 else 22.0     # West-Centre, ohm
X_CE = float(sys.argv[2]) if len(sys.argv) > 2 else 22.0     # Centre-East, ohm
KE_W = float(sys.argv[3]) if len(sys.argv) > 3 else 600e3    # MW.s
KE_C = float(sys.argv[4]) if len(sys.argv) > 4 else 800e3
KE_E = float(sys.argv[5]) if len(sys.argv) > 5 else 600e3

H_EQ = 40.0                                   # keep H, size via S_n
T_FLUX = 100.0                                # rotor fluxes ~constant at 0.2 Hz
T_SUB = 10.0                                  # dampers likewise
D_EQ = 0.0                                    # damping stated explicitly
S_W, S_C, S_E = KE_W/H_EQ, KE_C/H_EQ, KE_E/H_EQ        # MVA

text = open(SRC, encoding='utf-8').read()

# --- 1. new buses ------------------------------------------------------------
buses = ''.join(f'''  {{
   "name": "{n}",
   "P_W": 0,
   "Q_var": 0,
   "U_kV": 400
  }},
''' for n in ('CE_C', 'CE_E'))
text = text.replace(' "buses": [\n', ' "buses": [\n' + buses, 1)

# --- 2. new corridors --------------------------------------------------------
# bus_j first, deliberately: lines.py converts R/X from ohms using *bus_j*'s
# U_kV, and bus FR is declared at U_kV = 1.  Writing this corridor as FR -> CE_C
# makes Z_base 0.01 ohm, turns 22 ohm into 2200 pu, and silently leaves the
# chain open-circuited -- which is exactly what happened on the first attempt.
lines = ''.join(f'''  {{
   "bus_j": "{j}",
   "bus_k": "{k}",
   "R": {x/25.0:.5f},
   "X": {x:.5f},
   "B_pu": 0,
   "S_mva": 8000
  }},
''' for j, k, x in (('CE_C', 'FR', X_WC), ('CE_C', 'CE_E', X_CE)))
text = text.replace(' "lines": [\n', ' "lines": [\n' + lines, 1)


def equivalent(name, bus, s_mva):
    """A CE area equivalent: an aggregate, not a machine.

    The rotor flux time constants are deliberately long so the fluxes are
    effectively constant over the period of a 0.2 Hz mode -- i.e. a classical
    machine, which is the conventional form for an area equivalent and what
    `gencls` exists for.  Giving a 15-20 GVA aggregate one turbo-generator's
    round-rotor data instead is not harmless: SLACK_FR's T1q0 of 0.8 s puts the
    q-axis corner at 1.25 rad/s, right on top of a 0.2 Hz mode, and that single
    parameter was contributing 12 points of damping to it.

    Damping is then set explicitly through D, where it can be seen and justified,
    rather than arriving as a by-product of flux decay.
    """
    return f'''  {{
   "bus": "{bus}",
   "name": "{name}",
   "type": "genrou",
   "S_n": {s_mva*1e6:.1f},
   "F_n": 50,
   "H": {H_EQ},
   "X_d": 1.9,
   "X_q": 1.8,
   "X1d": 0.35,
   "T1d0": {T_FLUX},
   "X1q": 0.55,
   "T1q0": {T_FLUX},
   "X2d": 0.25,
   "T2d0": {T_SUB},
   "X2q": 0.25,
   "T2q0": {T_SUB},
   "R_a": 0.0025,
   "D": {D_EQ},
   "S_10": 0.05,
   "S_12": 0.2,
   "K_sec": 0,
   "K_delta": 0,
   "avr": {{
    "type": "sexs",
    "K_a": 100,
    "T_a": 1,
    "T_b": 1,
    "T_e": 0.2,
    "E_min": 0,
    "E_max": 5,
    "v_ref": 1
   }},
   "gov": {{
    "type": "tgov1",
    "R": 0.05,
    "T_1": 0.5,
    "V_max": 100,
    "V_min": -100,
    "T_2": 2.1,
    "T_3": 7,
    "D_t": 0,
    "p_c": 0
   }}
  }},
'''


# --- 3. new machines ---------------------------------------------------------
syns = equivalent('CE_CENTRE', 'CE_C', S_C) + equivalent('CE_EAST', 'CE_E', S_E)
text = text.replace(' "syns": [\n', ' "syns": [\n' + syns, 1)

# --- 4. shrink SLACK_FR to the West group, preserving its MW schedule ---------
S_OLD = 50000.0
m = re.search(r'"name":\s*"SLACK_FR".*?"S_n":\s*([\d.eE+]+)', text, re.S)
assert m, 'SLACK_FR not found'
text = text[:m.start(1)] + f'{S_W*1e6:.1f}' + text[m.end(1):]

# gov p_c is pu on the machine base: rescale so the MW setpoint is unchanged
m = re.search(r'"name":\s*"SLACK_FR".*?"p_c":\s*(-?[\d.eE+]+)', text, re.S)
assert m, 'SLACK_FR gov p_c not found'
p_c_old = float(m.group(1))
p_c_new = p_c_old*S_OLD/S_W
text = text[:m.start(1)] + f'{p_c_new:.6f}' + text[m.end(1):]

# SLACK_FR is now the West area equivalent, so it gets the same aggregate
# treatment as the other two -- it participates in the mode as strongly as they do.
anchor = text.index('"name": "SLACK_FR"')
for key, val in (('T1d0', T_FLUX), ('T1q0', T_FLUX),
                 ('T2d0', T_SUB), ('T2q0', T_SUB), ('D', D_EQ)):
    m = re.compile(rf'("{key}":\s*)(-?[\d.eE+]+)').search(text, anchor)
    assert m, f'SLACK_FR {key} not found'
    text = text[:m.start(2)] + f'{val}' + text[m.end(2):]

open(OUT, 'w', encoding='utf-8').write(text)
print(f'wrote {OUT}')
print(f'  West  SLACK_FR  S_n {S_W:8,.0f} MVA  H {H_EQ}  KE {KE_W:10,.0f} MW.s'
      f'   gov p_c {p_c_old:.5f} -> {p_c_new:.5f} pu ({p_c_old*S_OLD:,.0f} MW)')
print(f'  Centre CE_CENTRE S_n {S_C:8,.0f} MVA  H {H_EQ}  KE {KE_C:10,.0f} MW.s')
print(f'  East   CE_EAST   S_n {S_E:8,.0f} MVA  H {H_EQ}  KE {KE_E:10,.0f} MW.s')
print(f'  X_wc {X_WC} ohm  (P_max {160000/X_WC:,.0f} MW)'
      f'   X_ce {X_CE} ohm  (P_max {160000/X_CE:,.0f} MW)')
