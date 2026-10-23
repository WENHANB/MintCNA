#!/usr/bin/env python3
"""Does every axis title fit beside the panel it belongs to?

A title that overruns its own panel is the defect this catches -- it reads as
belonging to the neighbour, or runs off the page.

Runs are assigned to a panel by GEOMETRY, never by clustering on a fixed gap.
The gap approach merged "UMAP 1" from the left UMAP with "UMAP 1" from the right
into one 67 mm title, and the two stacked y titles of panel G into "log shiftlog
shift", every time the page size changed the spacing. A run belongs to the panel
it sits directly below (x title) or directly left of (y title), or to none.

Text positions come from svgtext.parse_runs, which honours `text-anchor`; svglite
anchors centred text at its middle, so reading x as the left edge shifts every
axis title right by half its own width.
"""
import os, re, sys, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import svgtext

F = sys.argv[1] if len(sys.argv) > 1 else "figure5-multiome/output/fig5_multiome.svg"
PW, PH = svgtext.page_size()
BELOW  = 14.0   # mm a horizontal title may sit under its panel's bottom edge
LEFTOF = 18.0   # mm a rotated title may sit left of its panel's left edge
MIN_PT = 6.5    # tick labels are 6.0 and below; titles are 7.0 and 7.4

svg = open(F).read()
S = PW / float(re.search(r"viewBox='0 0 ([\d.]+) ", svg).group(1))
runs  = svgtext.parse_runs(svg, S, min_fs=MIN_PT * S)
rects = svgtext.panel_rects(svg, S)

xs, ys = collections.defaultdict(list), collections.defaultdict(list)
for r in runs:
    x0, y0, x1, y1 = r["box"]
    cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
    # An IN-PANEL annotation is not an axis title. Panel F's "rho = 0.51" sits just
    # inside F's top edge, which is also within 14 mm of panel E's bottom edge, and
    # it was being glued onto E's x title. Anything whose centre lands inside a
    # plot rectangle belongs to that plot's interior; qa_text_overlap.py checks it.
    if any(p["x0"] <= cx <= p["x1"] and p["ytop"] <= cy <= p["ybot"] for p in rects):
        continue
    for i, p in enumerate(rects):
        if abs(r["rot"]) < 1:
            # an x title sits BELOW its own panel and horizontally inside it
            if 0 < cy - p["ybot"] < BELOW and p["x0"] - 2 <= cx <= p["x1"] + 2:
                xs[i].append(r); break
        elif abs(r["rot"] + 90) < 1:
            # a y title sits LEFT of its own panel and vertically inside it
            if 0 < p["x0"] - cx < LEFTOF and p["ytop"] - 2 <= cy <= p["ybot"] + 2:
                ys[i].append(r); break

def report(kind, groups, extent, avail):
    out = []
    for i, rs in sorted(groups.items()):
        rs.sort(key=lambda r: r["box"][0] if kind == "x" else -r["box"][3])
        lo = min(extent(r)[0] for r in rs); hi = max(extent(r)[1] for r in rs)
        out.append((kind, i, hi - lo, avail(rects[i]), "".join(r["txt"] for r in rs)))
    return out

rows  = report("y", ys, lambda r: (r["box"][1], r["box"][3]),
               lambda p: p["ybot"] - p["ytop"])
rows += report("x", xs, lambda r: (r["box"][0], r["box"][2]),
               lambda p: p["x1"] - p["x0"])

print(f"{len(rects)} panels; a title must fit the panel it belongs to\n")
print(f"{'panel':>5}  {'axis':>4}  {'title':<40} {'title mm':>9} {'panel mm':>9}   verdict")
bad = 0
for kind, i, ln, avail, txt in sorted(rows, key=lambda r: (r[1], r[0])):
    ok = ln <= avail + 0.15
    bad += (not ok)
    print(f"{i:>5}  {kind:>4}  {txt[:40]:<40} {ln:9.2f} {avail:9.2f}   "
          f"{'fits' if ok else 'OVERRUNS'}")
if not rows:
    print("  no axis titles matched a panel -- check the geometry, not the figure")
    bad = 1
print(f"\n  {bad} overrun(s)")
sys.exit(1 if bad else 0)
