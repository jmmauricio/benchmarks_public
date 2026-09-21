#!/usr/bin/env python
"""Solve the load flow and export per-generator power-flow results to JSON.

Runs the same steady-state initialisation as ``main.py ini()`` and writes the
*solved* operating point (the numbers ``report_gens`` / ``report_buses`` show):

  * generators: P_MW, Q_Mvar from the model algebraic states — ``p_g/q_g`` for
    synchronous machines, ``p_s/q_s`` for PV/wind inverters (both pu on S_n).
    These are the physical values (PV/wind ~unity PF; machines relaxed), unlike
    the raw PyPSA figures stored in the hjson ``results`` block.
  * buses: solved V_pu.

Output: ``esp_preblackout_pf.json`` — consumed by ``plot_sld_interactive.py``.

Run in the pydae env:  python export_pf.py [--hjson esp_preblackout_dev.hjson]
"""
import argparse
import json
import math
from pathlib import Path

from pydae.bps import BpsBuilder
from pydae.bps.utils.reporter import report_buses, report_gens
from pydae.core.builder.casadi_builder import CasadiBuilder
from pydae.core.model.casadi_model import CasadiModel

HERE = Path(__file__).resolve().parent


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hjson", default=str(HERE / "esp_preblackout_dev.hjson"))
    ap.add_argument("--xy0", default=str(HERE / "xy_0.json"))
    ap.add_argument("--out", default=str(HERE / "esp_preblackout_pf.json"))
    ap.add_argument("--md", default=str(HERE / "esp_preblackout_pf.md"),
                    help="markdown report_buses/report_gens tables")
    args = ap.parse_args()

    d = json.load(open(args.hjson))
    g = BpsBuilder(args.hjson, use_casadi=True)
    g.uz_jacs = False
    g.construct("esp_preblackout")
    m = CasadiModel(CasadiBuilder(g.sys_dict).build())
    m.ini({}, xy_0=args.xy0)

    def gv(name):
        try:
            return float(m.get_value(name))
        except Exception:
            return float("nan")

    gens = {}
    for s in d["syns"]:
        nm, s_mva = s["name"], s["S_n"] / 1e6
        gens[nm] = {"P_MW": round(gv(f"p_g_{nm}") * s_mva, 3),
                    "Q_Mvar": round(gv(f"q_g_{nm}") * s_mva, 3), "kind": "syn"}
    for p in d["pvs"]:
        nm, s_mva = p["name"], p["S_n"] / 1e6
        gens[nm] = {"P_MW": round(gv(f"p_s_{nm}") * s_mva, 3),
                    "Q_Mvar": round(gv(f"q_s_{nm}") * s_mva, 3), "kind": "pv"}
    for w in d.get("wecs", []):     # PMSM wind turbines (grid injection p_s/q_s, on S_n)
        nm, s_mva = w["name"], w["S_n"] / 1e6
        gens[nm] = {"P_MW": round(gv(f"p_s_{nm}") * s_mva, 3),
                    "Q_Mvar": round(gv(f"q_s_{nm}") * s_mva, 3), "kind": "wind"}

    buses = {}
    for b in d["buses"]:
        v = gv(f"V_{b['name']}")
        if not math.isnan(v):
            buses[b["name"]] = round(v, 4)

    json.dump({"gens": gens, "buses": buses}, open(args.out, "w"), indent=1)

    # markdown report_buses / report_gens tables (same numbers, human-readable)
    buses_md = report_buses(m, args.hjson, print_table=False)
    gens_md = report_gens(m, args.hjson, print_table=False)
    with open(args.md, "w") as f:
        f.write("# Pre-blackout power flow (solved)\n\n"
                f"Source model: `{Path(args.hjson).name}` — steady-state `ini()`.\n\n"
                "## Buses\n\n" + buses_md + "\n\n## Generators\n\n" + gens_md + "\n")

    print(f"exported {len(gens)} generators, {len(buses)} bus voltages -> {args.out}")
    print(f"markdown tables -> {args.md}")
    print(f"  PV/wind |Q|>0: {sum(1 for k,v in gens.items() if v['kind']=='pv' and abs(v['Q_Mvar'])>1e-6)}  "
          f"(should be ~0, unity PF)")
    # a couple of spot checks
    for nm in ("CCGT_ES620", "nuclear_ES424", "biomass_ES432", "solar_ES112"):
        if nm in gens:
            print(f"  {nm:16} P={gens[nm]['P_MW']:9.1f} MW  Q={gens[nm]['Q_Mvar']:9.1f} Mvar")


if __name__ == "__main__":
    main()
