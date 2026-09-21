"""Minimize slack power p_A0 by choosing the VSC A4 active-power reference.

The sweep in opt.py shows p_A0(p_vsc_A4_ref) is U-shaped: more VSC injection
reduces grid imports but raises VSC + line losses. We formulate it as

    min_{u, x, y_ini}   p_A0(x, y_ini, u, p)
        s.t.   builder._residual_fn(v_ini, p_ini) = 0
               u_b - u_a = 0,  u_c - u_a = 0       # balanced split
               u_min <= u_ph <= u_max

We reuse builder._residual_fn (the same function CasadiModel.ini() drives
through Newton) instead of rebuilding f, g, and the symbol map by hand.
p_ini is assembled with numeric values everywhere except the three control
slots, where we substitute the decision symbols.
"""
import json
import sys

import casadi as ca
import numpy as np

sys.path.insert(0, "/Users/jmmauricio/workspace/pydae/packages/pydae-uds/src")

from pydae.uds import UdsBuilder
from pydae.core.builder import CasadiBuilder


grid = UdsBuilder("uds8bus.hjson", use_casadi=True)
grid.construct("uds8bus")
builder = CasadiBuilder(grid.sys_dict).build()
sd = builder.sys_dict


# ---------------------------------------------------------------------------
# 1. Decision symbols. v_ini = [x; y_ini] is what _residual_fn expects.
# ---------------------------------------------------------------------------
u_ctrl_names = ["p_vsc_a_ref_A4", "p_vsc_b_ref_A4", "p_vsc_c_ref_A4"]
u_ctrl = ca.SX.sym("u_ctrl", len(u_ctrl_names))

x_sym     = ca.vertcat(*sd["x_list"])
y_ini_sym = ca.vertcat(*sd["y_ini_list"])
v_ini     = ca.vertcat(x_sym, y_ini_sym)


# ---------------------------------------------------------------------------
# 2. Build p_ini = [params; u_ini] with controls left symbolic, rest fixed.
# ---------------------------------------------------------------------------
p_A4 = 200e3
overrides = {
    "p_load_A4_a": p_A4 / 3,
    "p_load_A4_b": p_A4 / 3,
    "p_load_A4_c": p_A4 / 3,
}

p_ini_entries = [sd["params_dict"].get(s.name(), 0.0) for s in sd["p_list"]]
for s in sd["u_ini_list"]:
    name = s.name()
    if name in u_ctrl_names:
        p_ini_entries.append(u_ctrl[u_ctrl_names.index(name)])
    elif name in overrides:
        p_ini_entries.append(overrides[name])
    else:
        p_ini_entries.append(sd["u_ini_dict"].get(name, 0.0))
p_ini = ca.vertcat(*p_ini_entries)


# ---------------------------------------------------------------------------
# 3. Equality constraint and objective.
# ---------------------------------------------------------------------------
F_eq = builder._residual_fn(v_ini, p_ini)

p_A0_expr = ca.substitute(
    sd["h_dict"]["p_A0"],
    ca.vertcat(*sd["p_list"], *sd["u_ini_list"]),
    p_ini,
)

# Balanced split: tie b- and c-phase references to the a-phase.
balance = ca.vertcat(u_ctrl[1] - u_ctrl[0], u_ctrl[2] - u_ctrl[0])
G = ca.vertcat(F_eq, balance)


# ---------------------------------------------------------------------------
# 4. NLP.
# ---------------------------------------------------------------------------
w = ca.vertcat(u_ctrl, v_ini)
nlp = {"x": w, "f": p_A0_expr, "g": G}
solver = ca.nlpsol("solver", "ipopt", nlp,
                   {"ipopt.print_level": 3, "print_time": 0})


# ---------------------------------------------------------------------------
# 5. Bounds + warm start from xy_0.json.
# ---------------------------------------------------------------------------
P_VSC_MAX = 200e3 / 3
with open("xy_0.json") as fobj:
    xy_0 = json.load(fobj)

v_ini_names = [s.name() for s in sd["x_list"]] + [s.name() for s in sd["y_ini_list"]]
v_ini_seed  = [xy_0.get(n, 0.0) for n in v_ini_names]

sol = solver(
    x0  = [0.0] * 3 + v_ini_seed,
    lbx = [-P_VSC_MAX] * 3 + [-ca.inf] * v_ini.size1(),
    ubx = [+P_VSC_MAX] * 3 + [+ca.inf] * v_ini.size1(),
    lbg = [0.0] * G.size1(),
    ubg = [0.0] * G.size1(),
)

w_opt = np.array(sol["x"]).flatten()
p_vsc_A4_ref_opt = 3 * w_opt[0]
print(f"\np_vsc_A4_ref* = {p_vsc_A4_ref_opt/1e3:+.2f} kW")
print(f"p_A0*         = {float(sol['f'])/1e3:+.2f} kW")
