"""Minimize slack power p_A0 with three independent per-phase VSC references.

    min_{u_a, u_b, u_c, x, y_ini}   p_A0(x, y_ini, u, p)
        s.t.   builder._residual_fn(v_ini, p_ini) = 0
               u_min <= u_ph <= u_max,  ph in {a, b, c}

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
p_A4 = 100e3
overrides = {
    "p_load_A4_a": p_A4 * 1.0,
    "p_load_A4_b": p_A4 * 0.2,
    "p_load_A4_c": p_A4 * 0.2,
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
# 3. Equality constraint, objective, and per-phase voltage bounds.
# ---------------------------------------------------------------------------
F_eq = builder._residual_fn(v_ini, p_ini)

p_ini_subs_keys = ca.vertcat(*sd["p_list"], *sd["u_ini_list"])

p_A0_expr = ca.substitute(sd["h_dict"]["p_A0"], p_ini_subs_keys, p_ini)

V_BOUND_NAMES = [f"v_{ph}nm_{bus}"
                 for bus in ("A2", "A3", "A4")
                 for ph in ("a", "b", "c")]
v_bounds_expr = ca.substitute(
    ca.vertcat(*[sd["h_dict"][n] for n in V_BOUND_NAMES]),
    p_ini_subs_keys,
    p_ini,
)

V_NOM = 400.0 / np.sqrt(3)
V_MIN = 0.95 * V_NOM
V_MAX = 1.05 * V_NOM


# ---------------------------------------------------------------------------
# 4. NLP.
# ---------------------------------------------------------------------------
w = ca.vertcat(u_ctrl, v_ini)
G = ca.vertcat(F_eq, v_bounds_expr)
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
    lbg = [0.0] * F_eq.size1() + [V_MIN] * len(V_BOUND_NAMES),
    ubg = [0.0] * F_eq.size1() + [V_MAX] * len(V_BOUND_NAMES),
)

w_opt = np.array(sol["x"]).flatten()
p_a, p_b, p_c = w_opt[0], w_opt[1], w_opt[2]
print(f"\np_vsc_a_ref_A4* = {p_a/1e3:+.2f} kW")
print(f"p_vsc_b_ref_A4* = {p_b/1e3:+.2f} kW")
print(f"p_vsc_c_ref_A4* = {p_c/1e3:+.2f} kW")
print(f"sum             = {(p_a+p_b+p_c)/1e3:+.2f} kW")
print(f"p_A0*           = {float(sol['f'])/1e3:+.2f} kW")

g_opt = np.array(sol["g"]).flatten()
v_opt = g_opt[F_eq.size1():]
print(f"\nvoltage limits: [{V_MIN:.2f}, {V_MAX:.2f}] V (1 pu = {V_NOM:.2f} V)")
for name, v in zip(V_BOUND_NAMES, v_opt):
    print(f"  {name}: {v:.2f} V  ({v/V_NOM:.4f} pu)")
