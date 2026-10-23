#!/usr/bin/env python3
"""QA for Nature Communications figure PDFs, stdlib only.
Prints page size (mm), every font object (name / subtype / how it is embedded)
and a histogram of EFFECTIVE text sizes: size = |Tm scale| * Tf * |cm scale|,
so rotated y-axis titles and cairo's size-1 Tf convention are measured correctly.
  python3 pdf_text_sizes.py file.pdf [more.pdf ...]
"""
import re, sys, zlib, math, collections

def streams(raw):
    for m in re.finditer(rb'stream\r?\n(.*?)\r?\nendstream', raw, re.S):
        body = m.group(1)
        try:
            yield zlib.decompress(body)
        except Exception:
            yield body            # already uncompressed, or not Flate

NUM = rb'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?'
TM  = re.compile(rb'(' + NUM + rb')\s+(' + NUM + rb')\s+(' + NUM + rb')\s+(' + NUM + rb')\s+(' + NUM + rb')\s+(' + NUM + rb')\s+(Tm|cm)\b')
TF  = re.compile(rb'/([^\s/]+)\s+(' + NUM + rb')\s+Tf\b')
TOK = re.compile(rb'(?:' + TM.pattern + rb')|(?:' + TF.pattern + rb')|\b(q|Q|BT|ET)\b|(\[[^\]]*\]\s*TJ|\([^)]*\)\s*Tj)')

def sizes(raw):
    out = collections.Counter()
    for s in streams(raw):
        if b'Tf' not in s and b'Tj' not in s and b'TJ' not in s:
            continue
        # an inflated IMAGE stream is random bytes and can contain 'BT'/'Tf' by
        # chance (it did once: a 1 pt phantom run); a content stream is ASCII
        if sum(32 <= c < 127 or c in (9, 10, 13) for c in s[:4000]) < 0.95 * min(len(s), 4000):
            continue
        cm_stack = [1.0]; cm = 1.0; tf = 1.0; tm = 1.0; shown = False
        for m in TOK.finditer(s):
            g = m.groups()
            if g[6] is not None:                                   # Tm or cm
                a, b = float(g[0]), float(g[1]); sc = math.hypot(a, b)
                if g[6] == b'cm': cm *= sc
                else: tm = sc
            elif g[7] is not None:                                 # Tf
                tf = float(g[8])
            elif g[9] is not None:                                 # q Q BT ET
                op = g[9]
                if op == b'q': cm_stack.append(cm)
                elif op == b'Q': cm = cm_stack.pop() if len(cm_stack) > 1 else cm_stack[0]
                elif op == b'BT': tm = 1.0
            elif g[10] is not None:                                # text shown
                out[round(tm * tf * cm, 2)] += 1
    return out

def fonts(raw):
    rows = []
    for m in re.finditer(rb'/Type\s*/Font\b(.*?)(?=endobj)', raw, re.S):
        d = m.group(1)
        sub = re.search(rb'/Subtype\s*/(\w+)', d); bf = re.search(rb'/BaseFont\s*/([^\s/>]+)', d)
        if sub and sub.group(1) == b'Type0':
            continue          # composite wrapper; its descendant CIDFont carries the descriptor
        rows.append((bf.group(1).decode() if bf else '?', sub.group(1).decode() if sub else '?'))
    desc = {}
    for m in re.finditer(rb'/Type\s*/FontDescriptor\b(.*?)(?=endobj)', raw, re.S):
        d = m.group(1); fn = re.search(rb'/FontName\s*/([^\s/>]+)', d)
        kind = ('FontFile2 (TrueType)' if b'/FontFile2' in d else
                'FontFile3 (CFF/OpenType)' if b'/FontFile3' in d else
                'FontFile (Type1)' if b'/FontFile' in d else 'NOT EMBEDDED')
        if fn: desc[fn.group(1).decode()] = kind
    return rows, desc

for f in sys.argv[1:]:
    raw = open(f, 'rb').read()
    print('==', f)
    # cairo >= 1.17 keeps /Page and font dictionaries inside compressed object
    # streams, so page boxes and fonts are searched in the inflated streams too
    allbytes = raw + b'\n'.join(s for s in streams(raw) if b'/MediaBox' in s or b'/Font' in s)
    boxes = set(re.findall(rb'/MediaBox\s*\[\s*(' + NUM + rb')\s+(' + NUM + rb')\s+(' + NUM + rb')\s+(' + NUM + rb')', allbytes))
    for b in boxes:
        x0, y0, x1, y1 = map(float, b)
        print(f'   page {(x1-x0)/72*25.4:.1f} x {(y1-y0)/72*25.4:.1f} mm')
    rows, desc = fonts(allbytes)
    seen = set()
    for name, sub in rows:
        base = name.split('+', 1)[-1]
        if (base, sub) in seen: continue
        seen.add((base, sub))
        print(f'   font {base:28s} {sub:10s} {desc.get(name, desc.get(base, "?"))}')
    sz = sizes(raw)
    if sz:
        tot = sum(sz.values())
        print('   text sizes (pt: runs): ' + ', '.join(f'{k:g}: {v}' for k, v in sorted(sz.items())) + f'   [min {min(sz):g}, max {max(sz):g}, {tot} runs]')
    else:
        print('   text sizes: NO VECTOR TEXT FOUND (text is rasterised?)')
