#!/usr/bin/env python
"""Geographic single-line diagram of the esp pre-blackout pydae-bps model.

Draws the NUTS3 buses over the Spanish region map, the 400 kV lines (width
proportional to thermal rating), the synchronous units coloured by carrier and
sized by installed MVA, and the solar / onshore-wind fleets as translucent
circles with per-region MW labels. The France slack (continental-Europe
equivalent) and the PSS-tuned units are highlighted.

The .svg output carries a native ``<title>`` on every generator and line, so
hovering in a browser / SVG viewer shows the plant (or line) nominal power.
matplotlib cannot emit ``<title>`` directly, so each artist is tagged with a
``gid`` and the saved SVG is post-processed to insert the titles.

Requires geopandas + matplotlib (e.g. the conda ``base`` env, which has
geopandas but not pydae):

    python plot_sld.py
    python plot_sld.py --hjson esp_preblackout.hjson --alpha 0.6 --out esp_preblackout_sld

Note: tooltips are interactive — they appear only in the .svg opened in a
browser, not in the .png or a static preview.
"""
import argparse
import io
import json
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

import geopandas as gpd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

HERE = Path(__file__).resolve().parent
DEFAULT_HJSON = HERE / "esp_preblackout.hjson"
DEFAULT_REGIONS = HERE / "nuts3_shapes.geojson"
FALLBACK_REGIONS = Path(
    "/Users/jmmauricio/workspace/pypsa/pypsa-spain/resources/ES_test/nuts3_shapes.geojson"
)

# France single-equivalent node, placed north of Catalonia / the Pyrenees
FR_POS = (1.8, 43.6)

CARRIER_COL = {
    "nuclear": "#7b3fa0", "CCGT": "#d95f02", "hydro": "#1f78b4",
    "ror": "#66c2ff", "biomass": "#33a02c", "SLACK": "#e31a1c",
}
SOLAR_COL, WIND_COL = "#f5c000", "#1b9e77"


def load_network(path):
    """Read the model file (plain JSON, or hjson if the file uses hjson syntax)."""
    try:
        with open(path) as f:
            return json.load(f)
    except json.JSONDecodeError:
        import hjson  # only needed if the file is hjson-formatted (not in base env)
        with open(path) as f:
            return hjson.load(f)


def carrier(name):
    return name.split("_")[0]


def cap_label(mva):
    return f"{mva / 1000:.1f}G" if mva >= 950 else f"{mva:.0f}"


def aggregate(d, region_of):
    """Per-region installed MVA and per-plant lists for synchronous / solar / wind."""
    sync_mva = defaultdict(lambda: defaultdict(float))   # region -> carrier -> MVA
    sync_plants = defaultdict(list)                       # region -> [(name, MVA)]
    for s in d["syns"]:
        if s["name"] == "SLACK_FR":
            continue
        r, mva = region_of(s.get("bus", "")), s["S_n"] / 1e6
        sync_mva[r][carrier(s["name"])] += mva
        sync_plants[r].append((s["name"], mva))

    solar_mva, wind_mva = defaultdict(float), defaultdict(float)
    solar_plants, wind_plants = defaultdict(list), defaultdict(list)
    # solar in `pvs`; wind may be in `pvs` (pv_pq_ss) or `wecs` (pmsm_pq_ss)
    for p in d["pvs"] + d.get("wecs", []):
        r, mva = region_of(p["bus"]), p["S_n"] / 1e6
        if p["name"].startswith(("onwind", "wind")):
            wind_mva[r] += mva
            wind_plants[r].append((p["name"], mva))
        else:
            solar_mva[r] += mva
            solar_plants[r].append((p["name"], mva))
    return (sync_mva, sync_plants, solar_mva, solar_plants, wind_mva, wind_plants)


def tuned_units(d):
    """Machines whose PSS gain was tuned away from the K_stab=1 default.

    Falls back to the ten largest by S_n if none differ (so the highlight still
    makes sense on an untuned model)."""
    tuned = {s["name"] for s in d["syns"]
             if s["name"] != "SLACK_FR" and abs(s.get("pss", {}).get("K_stab", 1.0) - 1.0) > 1e-9}
    if not tuned:
        spain = [s for s in d["syns"] if s["name"] != "SLACK_FR"]
        tuned = {s["name"] for s in sorted(spain, key=lambda s: -s["S_n"])[:10]}
    return tuned


def plist(items):
    return "\n".join(f"  {n}: {m:,.0f} MVA" for n, m in sorted(items, key=lambda x: -x[1]))


def build_figure(d, gdf, alpha):
    cen = {code: (g.centroid.x, g.centroid.y) for code, g in gdf.geometry.items()}
    cen["FR"] = FR_POS
    t2r = {t["bus_j"]: t["bus_k"] for t in d["transformers"]}
    region_of = lambda bus: t2r.get(bus, bus)

    (sync_mva, sync_plants, solar_mva, solar_plants,
     wind_mva, wind_plants) = aggregate(d, region_of)
    tuned = tuned_units(d)
    tuned_regions = {region_of(next(s["bus"] for s in d["syns"] if s["name"] == nm)) for nm in tuned}

    tips = {}   # svg gid -> hover tooltip text
    fig, ax = plt.subplots(figsize=(13, 13))
    gdf.plot(ax=ax, color="0.96", edgecolor="0.8", lw=0.3, zorder=0)
    gdf.boundary.plot(ax=ax, color="0.75", lw=0.5, zorder=1)

    # transmission lines (grey; France tie in red), width ~ sqrt(rating)
    fr_lines = 0
    for i, l in enumerate(d["lines"]):
        ra, rb = region_of(l["bus_j"]), region_of(l["bus_k"])
        if ra not in cen or rb not in cen:
            continue
        (x0, y0), (x1, y1) = cen[ra], cen[rb]
        lw = 0.6 + 2.4 * (l["S_mva"] / 8000.0) ** 0.5
        is_fr = "FR" in (ra, rb)
        col, z = ("#e31a1c", 3) if is_fr else ("#555555", 2)
        (h,) = ax.plot([x0, x1], [y0, y1], color=col, lw=max(lw, 2.2) if is_fr else lw,
                       zorder=z, alpha=0.9 if is_fr else 0.6)
        fr_lines += is_fr
        gid = f"line_{i}"
        h.set_gid(gid)
        tips[gid] = (f"Line {l['bus_j']} — {l['bus_k']}\n"
                     f"  rating: {l['S_mva']:,.0f} MVA @ 400 kV\n  X = {l['X']:.2f} ohm")

    # synchronous machines: dominant-carrier dot sized by total MVA
    for r, carr in sync_mva.items():
        if r not in cen:
            continue
        x, y = cen[r]
        tot = sum(carr.values())
        dom = max(carr, key=carr.get)
        is_tuned = r in tuned_regions
        h = ax.scatter([x], [y], s=40 + tot * 0.9, color=CARRIER_COL.get(dom, "grey"),
                       edgecolors="k" if is_tuned else "none",
                       linewidths=2.2 if is_tuned else 0, zorder=5, alpha=alpha)
        gid = f"sync_{r}"
        h.set_gid(gid)
        tips[gid] = f"{r} — synchronous machines ({tot:,.0f} MVA total)\n" + plist(sync_plants[r])

    # solar (gold, offset left) and wind (teal, offset right): filled translucent circles
    for mva_map, plants, col, dx, lab_col, kind in [
        (solar_mva, solar_plants, SOLAR_COL, -0.10, "#b8860b", "solar"),
        (wind_mva, wind_plants, WIND_COL, +0.10, "#147a5a", "wind"),
    ]:
        for r, mva in mva_map.items():
            if r not in cen or mva <= 0:
                continue
            x, y = cen[r]
            h = ax.scatter([x + dx], [y], s=50 + mva * 0.22, color=col,
                           edgecolors="none", zorder=4, alpha=alpha)
            gid = f"{kind}_{r}"
            h.set_gid(gid)
            label = "solar PV" if kind == "solar" else "onshore wind"
            tips[gid] = f"{r} — {label} ({mva:,.0f} MVA)\n" + plist(plants[r])
            ax.annotate(cap_label(mva), (x + dx, y), textcoords="offset points",
                        xytext=(0, -10), ha="center", fontsize=5.5, color=lab_col, zorder=7)

    # France slack node
    s_fr = next(s["S_n"] for s in d["syns"] if s["name"] == "SLACK_FR") / 1e6
    h = ax.scatter(*FR_POS, s=900, marker="s", color="#e31a1c",
                   edgecolors="k", linewidths=1.5, zorder=6)
    h.set_gid("slack_FR")
    tips["slack_FR"] = ("SLACK_FR — France / continental Europe equivalent\n"
                        f"  slack machine: {s_fr:,.0f} MVA")
    ax.annotate("FR\n(slack / Europe)", FR_POS, textcoords="offset points", xytext=(0, 22),
                ha="center", fontsize=11, fontweight="bold", color="#e31a1c")

    # region codes on the tuned units
    for r in tuned_regions:
        if r in cen:
            ax.annotate(r, cen[r], textcoords="offset points", xytext=(4, 4),
                        fontsize=7.5, fontweight="bold")

    # legend + titles
    leg = [Line2D([0], [0], marker="o", color="w", markerfacecolor=c, markersize=11, label=k)
           for k, c in CARRIER_COL.items() if k != "SLACK"]
    leg += [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=SOLAR_COL, alpha=0.55,
               markeredgecolor="none", markersize=13,
               label=f"solar PV — {sum(solar_mva.values()) / 1000:.1f} GW (MW labels)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=WIND_COL, alpha=0.55,
               markeredgecolor="none", markersize=13,
               label=f"onshore wind — {sum(wind_mva.values()) / 1000:.1f} GW (MW labels)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="grey", markeredgecolor="k",
               markeredgewidth=2, markersize=12, label="PSS-tuned units"),
        Line2D([0], [0], color="#e31a1c", lw=2.5, label="France interconnection"),
        Line2D([0], [0], color="#555555", lw=1.5, label="400 kV line (width ~ rating)"),
    ]
    ax.legend(handles=leg, loc="lower left", fontsize=10, framealpha=0.95,
              title="Single-line diagram — esp pre-blackout")
    vre = d["pvs"] + d.get("wecs", [])
    n_wind = sum(1 for p in vre if p["name"].startswith(("onwind", "wind")))
    n_solar = len(vre) - n_wind
    ax.set_title(f"Spanish pre-blackout equivalent — {len(d['buses'])} buses, "
                 f"{len(d['lines'])} lines, {len(d['syns'])} synchronous machines, "
                 f"{n_solar} solar + {n_wind} wind plants", fontsize=13)
    ax.set_xlabel("lon"); ax.set_ylabel("lat")
    ax.set_xlim(-10, 5.2); ax.set_ylim(35.5, 44.4); ax.set_aspect(1.3)
    fig.tight_layout()
    stats = dict(sync=len(sync_mva), solar=len(solar_mva), wind=len(wind_mva),
                 fr_lines=int(fr_lines), tuned=sorted(tuned_regions))
    return fig, tips, stats


def save_with_tooltips(fig, tips, out_base):
    """Save PNG + SVG; inject a <title> child into every gid'd group in the SVG."""
    fig.savefig(out_base + ".png", dpi=140, bbox_inches="tight")
    buf = io.BytesIO()
    fig.savefig(buf, format="svg", bbox_inches="tight")
    buf.seek(0)
    ns = "http://www.w3.org/2000/svg"
    ET.register_namespace("", ns)
    ET.register_namespace("xlink", "http://www.w3.org/1999/xlink")
    tree = ET.parse(buf)
    added = 0
    for g in tree.getroot().iter(f"{{{ns}}}g"):
        gid = g.get("id")
        if gid in tips:
            title = ET.Element(f"{{{ns}}}title")
            title.text = tips[gid]
            g.insert(0, title)
            added += 1
    tree.write(out_base + ".svg", xml_declaration=True, encoding="utf-8")
    return added


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hjson", default=str(DEFAULT_HJSON), help="model file (default: %(default)s)")
    ap.add_argument("--regions", default=None,
                    help="NUTS3 shapes geojson (default: local nuts3_shapes.geojson)")
    ap.add_argument("--out", default=str(HERE / "esp_preblackout_sld"),
                    help="output basename, without extension")
    ap.add_argument("--alpha", type=float, default=0.60,
                    help="circle opacity 0-1 (default 0.60 = 40%% transparency)")
    args = ap.parse_args()

    regions = args.regions or (DEFAULT_REGIONS if DEFAULT_REGIONS.exists() else FALLBACK_REGIONS)
    d = load_network(args.hjson)
    gdf = gpd.read_file(regions).set_index("index")

    fig, tips, stats = build_figure(d, gdf, args.alpha)
    added = save_with_tooltips(fig, tips, args.out)

    print(f"regions  sync:{stats['sync']}  solar:{stats['solar']}  wind:{stats['wind']}  "
          f"FR lines:{stats['fr_lines']}")
    print(f"tooltips: {added}/{len(tips)} <title> elements injected")
    print(f"tuned regions: {stats['tuned']}")
    print(f"saved {args.out}.svg (+ .png)")


if __name__ == "__main__":
    main()
