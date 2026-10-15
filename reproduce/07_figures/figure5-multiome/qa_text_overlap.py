#!/usr/bin/env python3
"""Text that collides with other text, or that its own panel cuts off.

Eyeballing a 183 mm page at screen size hides both. Glyph positions come from
svgtext.parse_runs, which honours `text-anchor` -- without that every centred run
is placed half a word too far right.

Three checks:
  1. text on text   -- boxes overlapping by more than TOL mm in BOTH directions
  2. text off panel -- ggplot clips at the panel boundary, so an in-panel label
                       centred on a cluster near the edge is truncated silently.
                       "OPC-like" reached the page as "OPC-lik" and nothing said so.
  3. glyph floor    -- Nature's minimum is 5 pt AS RENDERED. Plotmath sub- and
                       superscripts come out at ~0.7x their parent, so a declared
                       7 pt superscript reaches the page at 4.9 pt. That has been
                       missed four times by reading the source instead of the SVG.
"""
import sys, itertools, re, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import svgtext

F = sys.argv[1] if len(sys.argv) > 1 else "figure5-multiome/output/fig5_multiome.svg"
PW, PH = svgtext.page_size()
TOL      = 0.25     # mm of genuine text-on-text overlap before it counts
TOL_CLIP = 0.40     # mm a run may sit past its panel edge before a glyph is lost

svg = open(F).read()
S = PW / float(re.search(r"viewBox='0 0 ([\d.]+) ", svg).group(1))
runs  = svgtext.parse_runs(svg, S)
rects = svgtext.panel_rects(svg, S)

print(f"{len(runs)} text runs, {len(rects)} panels; "
      f"overlap {TOL} mm, clip {TOL_CLIP} mm\n")

# ---- 1. text on text -------------------------------------------------------
hits = []
for a, b in itertools.combinations(runs, 2):
    ax, ay = a["box"], b["box"]
    dx = min(ax[2], ay[2]) - max(ax[0], ay[0])
    dy = min(ax[3], ay[3]) - max(ax[1], ay[1])
    if dx <= TOL or dy <= TOL:
        continue
    # same-word neighbours: identical baseline, touching horizontally only.
    # 1.2 mm, not 0.6: plotmath splits a title into one run per glyph group and a
    # wide glyph such as the Delta of "Dcopy number" over-reports its advance, so
    # consecutive groups appear to overlap by ~0.9 mm. A real collision on this
    # figure measured 5.4 mm, so the threshold has room.
    if abs(ax[1] - ay[1]) < 0.05 and dx < 1.2:
        continue
    # two runs at the same non-zero angle are glyphs of one rotated label, whose
    # axis-aligned hulls always intersect
    if abs(a["rot"]) > 1 and abs(a["rot"] - b["rot"]) < 1:
        continue
    hits.append((dx * dy, dx, dy, a, b))

# ---- 2. text clipped by its own panel --------------------------------------
clipped = []
for r in runs:
    if abs(r["rot"]) > 1:
        continue
    x0, y0, x1, y1 = r["box"]
    cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
    # a run BELONGS to the panel its centre falls in; tick labels and axis titles
    # sit outside every panel and are correctly skipped
    for p in rects:
        if p["x0"] <= cx <= p["x1"] and p["ytop"] <= cy <= p["ybot"]:
            over = max(p["x0"] - x0, x1 - p["x1"])
            if over > TOL_CLIP:
                clipped.append((r["txt"], over))
            break

if clipped:
    print(f"  {len(clipped)} run(s) CLIPPED by the panel edge:")
    for t, over in sorted(clipped, key=lambda c: -c[1]):
        print(f'    "{t}" past its panel by {over:.2f} mm')
    print()
else:
    print("  no text clipped by a panel edge")

# ---- 3. glyph floor --------------------------------------------------------
FLOOR = 5.0
small = sorted(((r["fs"] / S, r["txt"]) for r in runs), key=lambda q: q[0])[:3]
under = [q for q in small if q[0] < FLOOR - 1e-6]
if under:
    print(f"\n  {len(under)} run(s) BELOW the {FLOOR} pt floor as rendered:")
    for pt, t in under:
        print(f'    {pt:.2f} pt  "{t}"')
else:
    print(f"\n  glyph floor OK: smallest rendered text {small[0][0]:.2f} pt "
          f'("{small[0][1][:24]}"), floor {FLOOR} pt')

hits.sort(key=lambda h: -h[0])
if not hits:
    print("  no text-on-text overlaps")
else:
    print(f"\n  {len(hits)} overlapping pairs, worst first:\n")
    for area, dx, dy, a, b in hits[:25]:
        print(f"   {dx:5.2f} x {dy:4.2f} mm   \"{a['txt'][:28]:<28}\" vs "
              f"\"{b['txt'][:28]:<28}\"  (y~{a['box'][1]:.0f} mm from top)")
sys.exit(1 if (hits or clipped or under) else 0)
