"""

"""
import sys

import numpy as np

sys.path.insert(0, "/Users/jmmauricio/workspace/pydae/packages/pydae-uds/src")

from pydae.uds import UdsBuilder
from pydae.core.builder import CasadiBuilder, CasadiModel
from pydae.uds.utils.reports import report_v
from pydae.uds.utils.model2svg import model2svg


grid = UdsBuilder("uds8bus.hjson", use_casadi=True)
grid.construct("uds8bus")
model = CasadiModel(CasadiBuilder(grid.sys_dict).build())
model.report_params()


p_A4 = 200e3
params = {f"p_load_A4_{ph}": p_A4 / 3 for ph in ("a", "b", "c")}

for p_vsc_A4_ref in np.arange(-200e3, 200e3 + 1, 5e3):

    params.update({'p_vsc_a_ref_A4': p_vsc_A4_ref/3, 
                   'p_vsc_b_ref_A4': p_vsc_A4_ref/3, 
                   'p_vsc_c_ref_A4': p_vsc_A4_ref/3})

    model.ini(params, xy_0="xy_0.json")
    print(f"p_A4 = {p_A4/1e3:.0f} kW, p_A0 = {model.get_value('p_A0')/1e3:.0f} kW")


model.report_z()
