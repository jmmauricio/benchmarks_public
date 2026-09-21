"""Initialize and run the acdc_7bus example on the CasADi backend.

Prints a per-bus voltage report via pydae.uds.utils.reports.report_v after
ini() and again after a 1-second simulation. Run from this folder:

    python run.py

The pydae_dev conda env has pydae.uds resolving to a different repo, so we
prepend the workspace pydae-uds path before importing. Drop the sys.path
edit if pydae-uds is already first on the import path.
"""
import sys

import numpy as np

sys.path.insert(0, "/Users/jmmauricio/workspace/pydae/packages/pydae-uds/src")

from pydae.uds import UdsBuilder
from pydae.core.builder import CasadiBuilder, CasadiModel
from pydae.uds.utils.reports import report_v
from pydae.uds.utils.model2svg import model2svg

# 1. Construct the symbolic DAE on the CasADi backend and build the runtime.
grid = UdsBuilder("uds8bus.hjson", use_casadi=True)
grid.construct("uds8bus")
model = CasadiModel(CasadiBuilder(grid.sys_dict).build())
model.report_params()

# 2. Initialize from the saved operating-point seed with a 200 kW balanced
#    load split evenly across the three phases at bus A4. ini() walks the
#    fallback chain (newton rootfinder → _newton_solve → calc_ic) automatically.
p_A4 = 200e3
params = {f"p_load_A4_{ph}": p_A4 / 3 for ph in ("a", "b", "c")}
params.update({'K_acdc_A4': 100000, 'p_vsc_a_ref_A4': 50e3})
model.ini(params, xy_0="xy_0.json")
print("ini() converged\n--- bus voltages at t = 0 ---")
report_v(model, "uds8bus.hjson")

# 3. Integrate forward 1 second and re-report.
model.run(1.0, {}); model.post()
print(f"\n--- bus voltages at t = {model.t:.2f} s ---")
report_v(model, "uds8bus.hjson")

print(f"\nrun() ok; X finite={bool(np.isfinite(model.X).all())}, "
      f"Y finite={bool(np.isfinite(model.Y).all())}")

# 4. Annotate the topology SVG with the t=1 s operating point (bus colors by
#    voltage deviation, line colors by current loading, hover tooltips) and
#    write to a new file so the source template stays untouched.
s = model2svg(model, "uds8bus.hjson", "uds8bus.svg")
s.set_tooltips("uds8bus_tooltips.svg")
print("wrote uds8bus_tooltips.svg")
