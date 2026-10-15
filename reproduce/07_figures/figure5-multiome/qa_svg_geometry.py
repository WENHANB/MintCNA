#!/usr/bin/env python3
"""Measure Figure 4's panels from the EXPORTED SVG, in millimetres.

Everything before this measured what the code intends. This measures what is
actually on the page: svglite writes real coordinates, so the axis lines give the
plot rectangles directly. Axis lines in theme_cowplot are stroked polylines of
two points, one long horizontal (the x axis) and one long vertical (the y axis)
per panel, drawn in the theme's grey30 at 0.3 pt.
"""
import os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import svgtext, collections
F = sys.argv[1] if len(sys.argv) > 1 else "figure5-multiome/output/fig5_multiome.svg"
svg = open(F).read()
m = re.search(r"viewBox='0 0 ([\d.]+) ([\d.]+)'", svg)
Wpt, Hpt = float(m.group(1)), float(m.group(2))
PW, PH = svgtext.page_size()
sx, sy = PW / Wpt, PH / Hpt          # pt -> mm

segs = []
for mm_ in re.finditer(r"<polyline[^>]*points='([^']+)'[^>]*/>", svg):
    pts = [tuple(map(float, p.split(','))) for p in mm_.group(1).split()]
    if len(pts) != 2:
        continue
    (x1, y1), (x2, y2) = pts
    segs.append((x1, y1, x2, y2))

HOR = [s for s in segs if abs(s[1] - s[3]) < 0.4 and abs(s[2] - s[0]) > 60]
VER = [s for s in segs if abs(s[0] - s[2]) < 0.4 and abs(s[3] - s[1]) > 60]
print(f"canvas {Wpt:.1f} x {Hpt:.1f} pt = {PW} x {PH} mm")
print(f"long horizontal segments: {len(HOR)}   long vertical: {len(VER)}\n")

# pair each vertical (y axis) with the horizontal (x axis) that shares its corner
panels = []
for vx1, vy1, vx2, vy2 in VER:
    vtop, vbot = min(vy1, vy2), max(vy1, vy2)
    best = None
    for hx1, hy1, hx2, hy2 in HOR:
        if abs(hy1 - vbot) < 2.0 and abs(min(hx1, hx2) - vx1) < 2.0:
            best = (min(hx1, hx2), max(hx1, hx2)); break
    if best:
        panels.append((vx1, vtop, vbot, best[0], best[1]))

rows = []
for x, ytop, ybot, hl, hr in panels:
    left_mm  = hl * sx
    right_mm = hr * sx
    top_mm   = PH - ytop * sy
    bot_mm   = PH - ybot * sy
    rows.append((left_mm, right_mm, bot_mm, top_mm, right_mm - left_mm, top_mm - bot_mm))
rows.sort(key=lambda r: (-r[3], r[0]))

print(f"{'#':>2}  {'left':>7} {'right':>7} {'bottom':>7} {'top':>7}   {'width':>7} {'height':>7}   (mm)")
for i, r in enumerate(rows, 1):
    print(f"{i:>2}  {r[0]:7.2f} {r[1]:7.2f} {r[2]:7.2f} {r[3]:7.2f}   {r[4]:7.2f} {r[5]:7.2f}")

def spread(v): return max(v) - min(v)
if len(rows) >= 5:
    # The left column now holds A, the heatmap and TWO UMAPs, and one UMAP starts
    # at x = 73 -- a 60 mm split misfiled it as a right-column panel and the checks
    # then compared a UMAP against the volcano. The columns are separated at the
    # real gutter, RIGHT_X = 123.
    SPLIT = 123.0
    row1 = [r for r in rows if r[3] > 160]          # plot tops in the first row
    A = min(row1, key=lambda r: r[0]); B = max(row1, key=lambda r: r[0])
    right = sorted([r for r in rows if r[0] > SPLIT], key=lambda r: -r[3])
    umaps = sorted([r for r in rows if r[0] < SPLIT and r[3] < 60], key=lambda r: r[0])
    print("\n--- checks (measured on the page, tolerance 0.15 mm) ---")
    def chk(lab, v):
        s = spread(v)
        print(f"  {lab:<44} {'OK' if s < 0.15 else 'MISMATCH'}  (spread {s:.3f} mm)")
    # B and E are DELIBERATELY square and therefore narrower than D and F: both
    # carry the same quantity on both axes, and a rectangular panel would stretch
    # one against the other. So width and left/right edges are checked only on the
    # free panels, while height -- which sets the column's rhythm -- must hold for
    # all four, and B and E must actually be square.
    # right column, top to bottom: 0 B, 1 E, 2 F, 3 G. G was briefly two stacked
    # plots; it is one again (a per-gene dumbbell), so E and G are once more the
    # two free-shape panels. Their HEIGHTS are no longer equal -- G is the taller
    # slot because it carries 21 gene rows -- so only the widths and edges are
    # asserted, which is what the column alignment actually guarantees.
    # right column, top to bottom: 0 B, 1 E (cohort, square since 2026-09-05),
    # 2 F, 3 G. Three squares and one free-shape panel.
    # right column, top to bottom: 0 B, 1 E (motif volcano), 2 F -- all square
    free = []
    sq   = [right[i] for i in (0, 1, 2)]       # B, E, F
    # B and E are square and get taller boxes; D and F give that height up. So the
    # column no longer has one height -- what must hold is that the two squares
    # match each other and the two free panels match each other.
    chk("squares B,E,F: plot height",           [r[5] for r in sq])

    if len(umaps) == 2:
        chk("panel D UMAPs: equal width",       [r[4] for r in umaps])
        chk("panel D UMAPs: equal height",      [r[5] for r in umaps])
        chk("panel D UMAPs: share top edge",    [r[3] for r in umaps])
        chk("panel D UMAPs: share bottom edge", [r[2] for r in umaps])
        for i, r in enumerate(umaps):
            d = abs(r[4] - r[5])
            print(f"  {'panel D UMAP ' + str(i+1) + ': square':<44} "
                  f"{'OK' if d < 0.15 else 'MISMATCH'}  (w {r[4]:.2f}, h {r[5]:.2f})")
    for nm, r in zip(("b", "e", "f"), sq):
        d = abs(r[4] - r[5])
        print(f"  {'panel ' + nm + ': square (w == h)':<44} "
              f"{'OK' if d < 0.15 else 'MISMATCH'}  (w {r[4]:.2f}, h {r[5]:.2f}, diff {d:.3f} mm)")
    # With no free-shape panel there is no column reference; the squares must
    # simply share one horizontal centre (they are all centred in the same box).
    chk("squares B,E,F: share a horizontal centre",
        [(r[0]+r[1])/2 for r in sq])
    chk("row 1 (A,B): plot height",        [A[5], B[5]])
    chk("row 1 (A,B): top edge",           [A[3], B[3]])
    chk("row 1 (A,B): bottom edge",        [A[2], B[2]])
    gaps = [round(right[i][2] - right[i+1][3], 2) for i in range(len(right)-1)]
    print(f"  right column vertical gaps between plots: {gaps} mm")
