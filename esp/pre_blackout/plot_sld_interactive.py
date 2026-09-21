#!/usr/bin/env python
"""Interactive geographic single-line diagram of the esp pre-blackout model (Plotly).

Produces a self-contained ``esp_preblackout_sld.html`` you can open in any
browser. Features:

  * hover on any generator -> name, bus, S_n, P, Q and current;
  * a "size by" menu to drive the circle radius from installed capacity (S_n),
    active power P, reactive power |Q| or current |I|;
  * a "labels" menu to switch labels off / region code / active power;
  * click legend entries to show/hide carriers (nuclear, CCGT, hydro, solar,
    wind, ...); double-click to isolate one;
  * standard Plotly zoom / pan / box-select / PNG export.

Marker area scales with the chosen metric (radius ~ sqrt(value)). Values are
read from the model's ``results`` block (per-unit P_MW, Q_Mvar, bus V_pu).

Run in an env with geopandas + plotly (e.g. conda ``base``):

    python plot_sld_interactive.py
    python plot_sld_interactive.py --hjson esp_preblackout.hjson --out esp_preblackout_sld
"""
import argparse
import json
import math
from collections import defaultdict
from pathlib import Path

import geopandas as gpd
import numpy as np
import plotly.graph_objects as go

HERE = Path(__file__).resolve().parent
DEFAULT_HJSON = HERE / "esp_preblackout.hjson"
DEFAULT_REGIONS = HERE / "nuts3_shapes.geojson"
FALLBACK_REGIONS = Path(
    "/Users/jmmauricio/workspace/pypsa/pypsa-spain/resources/ES_test/nuts3_shapes.geojson"
)
FR_POS = (1.8, 43.6)
U_KV = 400.0

CARRIER_COL = {"nuclear": "#7b3fa0", "CCGT": "#d95f02", "hydro": "#1f78b4",
               "ror": "#66c2ff", "biomass": "#33a02c", "solar": "#f5c000", "wind": "#1b9e77"}
# small per-carrier offsets so different carriers at the same region don't overlap
OFFSET = {"nuclear": (0.0, 0.09), "CCGT": (0.09, 0.0), "hydro": (0.0, -0.09),
          "ror": (-0.09, 0.0), "biomass": (0.07, 0.07), "solar": (-0.12, -0.06),
          "wind": (0.12, -0.06)}
METRICS = [("S_n", "installed capacity [MVA]"), ("P", "active power P [MW]"),
           ("Q", "reactive power |Q| [Mvar]"), ("I", "current |I| [kA]")]
SMIN, SMAX = 6.0, 42.0


def load_network(path):
    try:
        with open(path) as f:
            return json.load(f)
    except json.JSONDecodeError:
        import hjson
        with open(path) as f:
            return hjson.load(f)


def carrier(name):
    c = name.split("_")[0]
    return "wind" if c in ("onwind", "wind") else ("solar" if c == "solar" else c)


def build_records(d, pf=None):
    """One record per generator with position + metric values.

    If ``pf`` (the solved load flow from export_pf.py) is given, P/Q/V are the
    physical power-flow values (PV/wind unity PF, machines relaxed). Otherwise
    fall back to the hjson ``results`` block, zeroing PV/wind Q and clamping
    synchronous Q to the MVA capability circle to strip the raw PyPSA artifacts.
    """
    sn = {g["name"]: g["S_n"] / 1e6 for g in d["syns"] + d["pvs"]}
    vpu = dict(pf["buses"]) if pf else {}
    for b in d["results"]["buses"]:
        vpu.setdefault(b["name"], b["V_pu"])
    gpf = pf["gens"] if pf else {}
    recs = []
    for src in ("syns", "pvs"):
        for g in d["results"][src]:
            name, region = g["name"], g["bus"]
            if name == "SLACK_FR":
                continue
            cat = carrier(name)
            S = sn.get(name, 0.0)
            if name in gpf:                                  # solved load-flow values
                P, Q = float(gpf[name]["P_MW"]), float(gpf[name]["Q_Mvar"])
            elif src == "pvs":                               # fallback: unity PF
                P, Q = float(g["P_MW"]), 0.0
            else:                                            # fallback: clamp to capability
                P = float(g["P_MW"])
                qmax = math.sqrt(max(S ** 2 - P ** 2, 0.0))
                Q = max(-qmax, min(qmax, float(g["Q_Mvar"])))
            v = vpu.get(region, 1.0)
            I = math.hypot(P, Q) / (math.sqrt(3) * U_KV * max(v, 0.1))  # kA
            recs.append(dict(name=name, region=region, cat=cat, S_n=S, P=P, Q=Q, I=I))
    return recs


def _metric_val(r, m):
    """Value used for marker sizing: |Q| for reactive, non-negative otherwise."""
    return abs(r["Q"]) if m == "Q" else max(r[m], 0.0)


def size_arrays(recs_by_cat):
    """For each metric, a size array per category (radius ~ sqrt(value))."""
    vmax = {m: max((_metric_val(r, m) for cat in recs_by_cat.values() for r in cat),
                   default=1.0) for m, _ in METRICS}
    out = {}
    for m, _ in METRICS:
        vm = vmax[m] or 1.0
        out[m] = {cat: [SMIN + (SMAX - SMIN) * math.sqrt(_metric_val(r, m) / vm)
                        for r in cat_recs]
                  for cat, cat_recs in recs_by_cat.items()}
    return out


def _compute_vpu_by_region(d, pf, region_of):
    vpu = dict(pf["buses"]) if pf else {}
    for b in d["results"]["buses"]:
        vpu.setdefault(b["name"], b["V_pu"])
    region_v = defaultdict(list)
    for bus, v in vpu.items():
        region_v[region_of(bus)].append(v)
    return {r: np.mean(vs) for r, vs in region_v.items()}


def _compute_region_power(d, pf, region_of):
    gpf = pf["gens"] if pf else {}
    region_p = defaultdict(float)
    region_q = defaultdict(float)
    for src in ("syns", "pvs"):
        for g in d["results"][src]:
            name, region = g["name"], g["bus"]
            if name == "SLACK_FR":
                continue
            if name in gpf:
                P, Q = float(gpf[name]["P_MW"]), float(gpf[name]["Q_Mvar"])
            elif src == "pvs":
                P, Q = float(g["P_MW"]), 0.0
            else:
                P = float(g["P_MW"])
                sn = float(g.get("S_n", 0)) / 1e6
                qmax = math.sqrt(max(sn**2 - P**2, 0.0))
                Q = max(-qmax, min(qmax, float(g["Q_Mvar"])))
            region_p[region] += P
            region_q[region] += Q
    return dict(region_p), dict(region_q)


def _compute_line_flows(d, pf, region_of):
    vpu = {}
    theta = {}
    if pf:
        for bus, v in pf["buses"].items():
            vpu[bus] = v
    for b in d["results"]["buses"]:
        vpu.setdefault(b["name"], b["V_pu"])
        theta.setdefault(b["name"], b.get("theta_deg", 0.0))
    flows = []
    for l in d["lines"]:
        bj, bk = l["bus_j"], l["bus_k"]
        vi = vpu.get(bj, 1.0)
        vj = vpu.get(bk, 1.0)
        ti = math.radians(theta.get(bj, 0.0))
        tj = math.radians(theta.get(bk, 0.0))
        X = l["X"]
        P = vi * vj * math.sin(ti - tj) / X if X > 0 else 0.0
        Q = (vi**2 - vi * vj * math.cos(ti - tj)) / X if X > 0 else 0.0
        S = math.hypot(P, Q)
        flows.append(dict(bus_j=bj, bus_k=bk, P=P, Q=Q, S=S, S_mva=l["S_mva"]))
    return flows


def _diverging_color(val, val_abs_max, pos_color="red", neg_color="blue"):
    t = val / val_abs_max if val_abs_max > 0 else 0
    t = max(-1.0, min(1.0, t))
    if t >= 0:
        r, g, b = 220, int(80 * (1 - t)), int(60 * (1 - t))
    else:
        at = abs(t)
        r, g, b = int(60 * (1 - at)), int(80 * (1 - at)), 220
    return f"rgba({r},{g},{b},0.5)"


def _voltage_color(v, v_min, v_max):
    t = (v - v_min) / (v_max - v_min) if v_max > v_min else 0.5
    t = max(0.0, min(1.0, t))
    r = int(255 * t)
    b = int(255 * (1 - t))
    return f"rgba({r},{int(80*(1-abs(2*t-1)))},{b},0.45)"


def boundary_trace(gdf, fill_mode=None, data_by_region=None):
    traces = []
    vals = []
    if data_by_region:
        vals = [data_by_region.get(code, 0) for code in gdf.index]
    v_min = min(vals) if vals else 0
    v_max = max(vals) if vals else 1
    v_abs_max = max(abs(v_min), abs(v_max)) if vals else 1
    for code, geom in gdf.geometry.items():
        polys = geom.geoms if geom.geom_type == "MultiPolygon" else [geom]
        for poly in polys:
            x, y = poly.exterior.xy
            xs, ys = list(x), list(y)
            if fill_mode == "voltage" and code in (data_by_region or {}):
                color = _voltage_color(data_by_region[code], v_min, v_max)
                ht = f"{code}<br>V = {data_by_region[code]:.4f} pu<extra></extra>"
            elif fill_mode in ("active_power", "reactive_power") and code in (data_by_region or {}):
                color = _diverging_color(data_by_region[code], v_abs_max)
                unit = "MW" if fill_mode == "active_power" else "Mvar"
                ht = f"{code}<br>{data_by_region[code]:+.1f} {unit}<extra></extra>"
            else:
                color = "rgba(200,200,200,0.25)"
                ht = None
            traces.append(go.Scatter(
                x=xs, y=ys, fill="toself",
                fillcolor=color,
                line=dict(color="#cccccc", width=0.6),
                hoverinfo="skip" if ht is None else "all",
                showlegend=False, name="regions",
                hovertemplate=ht,
            ))
    return traces


def build_figure(hjson_path=None, regions_path=None, pf_path=None,
                 province_fill=None, line_color_by=None,
                 size_by="S_n", labels="off"):
    hjson_path = Path(hjson_path or DEFAULT_HJSON)
    regions_path = Path(regions_path) if regions_path else (DEFAULT_REGIONS if DEFAULT_REGIONS.exists() else FALLBACK_REGIONS)
    pf_path = Path(pf_path) if pf_path else HERE / "esp_preblackout_pf.json"

    d = load_network(hjson_path)
    pf = json.load(open(pf_path)) if pf_path.exists() else None
    gdf = gpd.read_file(regions_path).set_index("index")
    cen = {code: (g.centroid.x, g.centroid.y) for code, g in gdf.geometry.items()}
    cen["FR"] = FR_POS
    region_of = (lambda b: {t["bus_j"]: t["bus_k"] for t in d["transformers"]}.get(b, b))

    recs = build_records(d, pf)
    recs_by_cat = defaultdict(list)
    for r in recs:
        if r["region"] in cen:
            recs_by_cat[r["cat"]].append(r)
    sizes = size_arrays(recs_by_cat)

    vpu_by_region = _compute_vpu_by_region(d, pf, region_of)
    region_p, region_q = _compute_region_power(d, pf, region_of)
    line_flows = _compute_line_flows(d, pf, region_of)

    fill_data = None
    if province_fill == "voltage":
        fill_data = vpu_by_region
    elif province_fill == "active_power":
        fill_data = region_p
    elif province_fill == "reactive_power":
        fill_data = region_q

    fig = go.Figure()
    for trace in boundary_trace(gdf, fill_mode=province_fill, data_by_region=fill_data):
        fig.add_trace(trace)

    flow_key = {"active_power": "P", "reactive_power": "Q", "apparent_power": "S"}.get(line_color_by)
    flow_unit = {"active_power": "MW", "reactive_power": "Mvar", "apparent_power": "MVA"}.get(line_color_by, "")
    flow_abs_max = 1.0
    if flow_key:
        flow_vals = [f[flow_key] for f in line_flows]
        flow_abs_max = max(abs(min(flow_vals)), abs(max(flow_vals))) if flow_vals else 1.0

    gx, gy, rx, ry = [], [], [], []
    mid_x, mid_y, mid_txt, mid_color = [], [], [], []
    for i, l in enumerate(d["lines"]):
        ra, rb = region_of(l["bus_j"]), region_of(l["bus_k"])
        if ra not in cen or rb not in cen:
            continue
        (x0, y0), (x1, y1) = cen[ra], cen[rb]
        is_fr = "FR" in (ra, rb)
        lf = line_flows[i] if i < len(line_flows) else {}
        if flow_key and not is_fr:
            fig.add_trace(go.Scatter(
                x=[x0, x1, None], y=[y0, y1, None], mode="lines",
                line=dict(color=_diverging_color(lf[flow_key], flow_abs_max), width=2),
                hoverinfo="skip", showlegend=False, name="lines"))
        else:
            (rx if is_fr else gx).extend([x0, x1, None])
            (ry if is_fr else gy).extend([y0, y1, None])
        mid_x.append((x0 + x1) / 2); mid_y.append((y0 + y1) / 2)
        val_str = f"<br>{flow_key}: {lf[flow_key]:+,.1f} {flow_unit}" if flow_key else ""
        mid_color.append(_diverging_color(lf[flow_key], flow_abs_max) if flow_key else "rgba(0,0,0,0)")
        mid_txt.append(f"Line {l['bus_j']} — {l['bus_k']}<br>rating {l['S_mva']:,.0f} MVA "
                       f"@ 400 kV<br>X = {l['X']:.2f} ohm{val_str}")

    fig.add_trace(go.Scatter(x=gx, y=gy, mode="lines", line=dict(color="#888", width=1),
                             hoverinfo="skip", name="400 kV lines", legendgroup="grid"))
    fig.add_trace(go.Scatter(x=rx, y=ry, mode="lines", line=dict(color="#e31a1c", width=2.5),
                             hoverinfo="skip", name="France interconnection", legendgroup="grid"))
    fig.add_trace(go.Scatter(x=mid_x, y=mid_y, mode="markers",
                             marker=dict(size=7, color=mid_color, line=dict(width=0)),
                             text=mid_txt, hovertemplate="%{text}<extra></extra>",
                             showlegend=False, name="line info"))

    gen_idx = []
    cat_order = [c for c in CARRIER_COL if c in recs_by_cat]
    for cat in cat_order:
        cat_recs = recs_by_cat[cat]
        ox, oy = OFFSET[cat]
        xs = [cen[r["region"]][0] + ox for r in cat_recs]
        ys = [cen[r["region"]][1] + oy for r in cat_recs]
        cd = [[r["name"], r["region"], r["S_n"], r["P"], r["Q"], r["I"]] for r in cat_recs]
        if labels == "region":
            txt = [r["region"] for r in cat_recs]
            mode = "markers+text"
        elif labels == "power":
            txt = [f"{r['P']:.0f}" for r in cat_recs]
            mode = "markers+text"
        else:
            txt = [""] * len(cat_recs)
            mode = "markers"
        gen_idx.append(len(fig.data))
        fig.add_trace(go.Scatter(
            x=xs, y=ys, mode=mode, name=cat, legendgroup="gen",
            marker=dict(size=sizes[size_by][cat], color=CARRIER_COL[cat],
                        opacity=0.6, line=dict(width=0)),
            customdata=cd,
            hovertemplate=("<b>%{customdata[0]}</b><br>bus %{customdata[1]}<br>"
                           "S_n %{customdata[2]:,.0f} MVA<br>P %{customdata[3]:,.1f} MW<br>"
                           "Q %{customdata[4]:,.1f} Mvar<br>I %{customdata[5]:,.2f} kA"
                           "<extra></extra>"),
            text=txt, textposition="top center", textfont=dict(size=8)))

    s_fr = next(s["S_n"] for s in d["syns"] if s["name"] == "SLACK_FR") / 1e6
    fig.add_trace(go.Scatter(x=[FR_POS[0]], y=[FR_POS[1]], mode="markers+text",
                             marker=dict(size=26, color="#e31a1c", symbol="square",
                                         line=dict(color="black", width=1.5)),
                             text=["FR (slack / Europe)"], textposition="top center",
                             textfont=dict(size=11, color="#e31a1c"), name="France (slack)",
                             legendgroup="grid",
                             hovertemplate=f"<b>SLACK_FR</b><br>France / Europe equivalent<br>"
                                           f"slack machine {s_fr:,.0f} MVA<extra></extra>"))

    n_wind = sum(len(v) for k, v in recs_by_cat.items() if k == "wind")
    n_solar = sum(len(v) for k, v in recs_by_cat.items() if k == "solar")

    annotations = []
    if province_fill == "voltage":
        v_vals = [vpu_by_region[code] for code in gdf.index if code in vpu_by_region]
        v_lo, v_hi = (min(v_vals), max(v_vals)) if v_vals else (0.95, 1.05)
        annotations = [
            dict(text=f"V (pu): {v_lo:.3f}", x=1.0, xref="paper", xanchor="right",
                 y=0.0, yref="paper", yanchor="bottom", showarrow=False,
                 font=dict(size=10, color="blue")),
            dict(text=f"{v_hi:.3f}", x=1.0, xref="paper", xanchor="right",
                 y=1.0, yref="paper", yanchor="top", showarrow=False,
                 font=dict(size=10, color="red")),
        ]
    elif province_fill in ("active_power", "reactive_power"):
        pv = fill_data or {}
        p_lo, p_hi = min(pv.values()) if pv else 0, max(pv.values()) if pv else 1
        unit = "MW" if province_fill == "active_power" else "Mvar"
        annotations = [
            dict(text=f"{unit}: {p_lo:+,.0f}", x=1.0, xref="paper", xanchor="right",
                 y=0.0, yref="paper", yanchor="bottom", showarrow=False,
                 font=dict(size=10, color="blue")),
            dict(text=f"{p_hi:+,.0f}", x=1.0, xref="paper", xanchor="right",
                 y=1.0, yref="paper", yanchor="top", showarrow=False,
                 font=dict(size=10, color="red")),
        ]

    fig.update_layout(
        title=(f"Spanish pre-blackout equivalent — {len(d['buses'])} buses, {len(d['lines'])} lines, "
               f"{len(d['syns'])} machines, {n_solar} solar + {n_wind} wind"),
        annotations=annotations,
        legend=dict(title="carrier / grid  (click to toggle)", itemsizing="constant"),
        plot_bgcolor="white", height=1000, margin=dict(t=90))
    fig.update_xaxes(visible=False)
    fig.update_yaxes(visible=False, scaleanchor="x", scaleratio=1.30)

    return fig


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hjson", default=str(DEFAULT_HJSON))
    ap.add_argument("--regions", default=None)
    ap.add_argument("--pf", default=str(HERE / "esp_preblackout_pf.json"),
                    help="solved load flow from export_pf.py (P/Q/V); "
                         "falls back to the hjson results block if absent")
    ap.add_argument("--out", default=str(HERE / "esp_preblackout_sld"))
    args = ap.parse_args()

    fig = build_figure(args.hjson, args.regions, args.pf)

    out = args.out + ".html"
    fig.write_html(out, include_plotlyjs=True, full_html=True,
                   config={"displaylogo": False, "responsive": True})
    print(f"saved {out}")


if __name__ == "__main__":
    main()
