"""Anchor-aware parsing of svglite text runs, shared by the QA scripts.

svglite writes a run's ANCHOR point in x/y, and centred text carries
`text-anchor='middle'` as an ATTRIBUTE, not as a style property. Reading x as the
left edge therefore shifts every centred run right by half its own width -- half a
word for a direct label, which is enough to invent a panel overrun that is not
there and to hide a text collision that is. Both checkers parse through here so
they cannot disagree about where a glyph sits.

`textLength` is authoritative: svglite emits the advance width it actually used,
so no font metric is estimated.
"""
import re, html, math

ANCHOR = {"start": 0.0, "middle": 0.5, "end": 1.0}

def parse_runs(svg, S, min_fs=0.0):
    """Every text run as dict(txt, fs, rot, box) with box in mm, y from page top."""
    runs = []
    for t in re.finditer(r"<text([^>]*)>([^<]*)</text>", svg):
        at, txt = t.group(1), html.unescape(t.group(2))
        if not txt.strip():
            continue
        fs = float(re.search(r"font-size: *([\d.]+)", at).group(1)) * S
        if fs < min_fs:
            continue
        tl = re.search(r"textLength='([\d.]+)", at)
        tr = re.search(r"transform='translate\(([-\d.]+),([-\d.]+)\)"
                       r"(?:\s*rotate\(([-\d.]+)\))?'", at)
        xm, ym = re.search(r"\bx='([-\d.]+)'", at), re.search(r"\by='([-\d.]+)'", at)
        if tr:
            x, y = float(tr.group(1)) * S, float(tr.group(2)) * S
            rot = float(tr.group(3)) if tr.group(3) else 0.0
        elif xm and ym:
            x, y, rot = float(xm.group(1)) * S, float(ym.group(1)) * S, 0.0
        else:
            continue
        L = float(tl.group(1)) * S if tl else 0.5 * fs * len(txt)
        am = re.search(r"text-anchor='([a-z]+)'", at)
        off = ANCHOR.get(am.group(1) if am else "start", 0.0) * L   # anchor -> run start
        # The start of the run, walking BACK along the text direction from the anchor.
        # An axis-aligned box is only meaningful for axis-aligned text: a 45-degree
        # run's axis-aligned hull is far larger than its glyphs, so those keep a true
        # rotated hull and are compared only across angles.
        if abs(rot + 90) < 1:                     # runs upward: toward decreasing y
            sy = y + off
            box = (x - fs * 0.8, sy - L, x + fs * 0.25, sy)
        elif abs(rot) < 1:                        # runs rightward
            sx = x - off
            box = (sx, y - fs * 0.8, sx + L, y + fs * 0.25)
        else:
            c, sn = math.cos(math.radians(rot)), math.sin(math.radians(rot))
            pts = [(-off, 0), (L - off, 0), (-off, -fs * 0.8), (L - off, -fs * 0.8)]
            xs = [x + px * c - py * sn for px, py in pts]
            ys = [y + px * sn + py * c for px, py in pts]
            box = (min(xs), min(ys), max(xs), max(ys))
        runs.append(dict(txt=txt, box=box, fs=fs, rot=rot, L=L,
                         start=(x - off if abs(rot) < 1 else x,
                                y + off if abs(rot + 90) < 1 else y)))
    return runs

def page_size(code_dir=None):
    """PAGE_W and PAGE_H in mm, parsed out of grid_fig5.R.

    The assembly and the QA have to agree on the page or the QA is measuring a
    figure that does not exist -- both hardcoded 170 while the page became 195.
    """
    import os
    d = code_dir or os.path.dirname(os.path.abspath(__file__))
    src = open(os.path.join(d, "grid_fig5.R")).read()
    g = lambda k: float(re.search(rf"^{k}\s*<-\s*([\d.]+)", src, re.M).group(1))
    return g("PAGE_W"), g("PAGE_H")

def panel_rects(svg, S, minlen=60.0):
    """Plot areas, found as an L of a tall axis line meeting a long one at bottom-left."""
    segs = []
    for m in re.finditer(r"<polyline[^>]*points='([^']+)'[^>]*/>", svg):
        pts = [tuple(map(float, q.split(','))) for q in m.group(1).split()]
        if len(pts) == 2:
            segs.append(pts)
    ver = [q for q in segs if abs(q[0][0]-q[1][0]) < .4 and abs(q[1][1]-q[0][1]) > minlen]
    hor = [q for q in segs if abs(q[0][1]-q[1][1]) < .4 and abs(q[1][0]-q[0][0]) > minlen]
    rects = []
    for (vx, vy1), (_, vy2) in ver:
        vtop, vbot = min(vy1, vy2), max(vy1, vy2)
        for (hx1, hy1), (hx2, hy2) in hor:
            if abs(hy1 - vbot) < 2 and abs(min(hx1, hx2) - vx) < 2:
                rects.append(dict(x0=min(hx1, hx2)*S, x1=max(hx1, hx2)*S,
                                  ytop=vtop*S, ybot=vbot*S)); break
    rects.sort(key=lambda p: (p["ytop"], p["x0"]))
    return rects
