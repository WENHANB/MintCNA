# Bin a fragments.tsv.gz onto a fixed grid, one row per (cell, chr, bin, count).
#
# Counting rule -- the decisions, stated once:
#   1. ONE per fragment. Column 5 counts PCR-duplicate read pairs, which is
#      sequencing redundancy, not extra molecules.
#   2. Placed by MIDPOINT, so a fragment lands in exactly one bin. Counting both
#      Tn5 ends instead would double every count and split the ~0.2% of fragments
#      that straddle a boundary; at 500 kb neither changes anything.
#   3. Only barcodes in CELLS, with at least MINFRAGS fragments. The depth filter
#      is applied to the fragment file itself so the count in column 3 of the cell
#      list (which is ArchR's, over all chromosomes) is not what gates it.
#   4. Fragments whose MIDPOINT falls in an ENCODE blacklist interval are dropped.
#      Masking at fragment level rather than deleting the whole bin keeps the
#      other 89% of a partially-blacklisted 500 kb window; the reference nuclei
#      lose the same positions, so the shortened effective length cancels in the
#      ratio.
#
# Usage: zcat frags.gz | awk -F'\t' -v CELLS=... -v BLACKLIST=... -v W=500000 \
#                            -v MINFRAGS=3000 -f this.awk

BEGIN {
  if (W == "") W = 500000
  if (MINFRAGS == "") MINFRAGS = 0

  while ((getline line < CELLS) > 0) {
    split(line, f, "\t")
    if (f[1] == "barcode") continue
    keep[f[1]] = 1
  }
  close(CELLS)

  # Blacklist intervals bucketed by bin, so the per-fragment check looks at a
  # handful of candidates instead of scanning all ~600 of them.
  if (BLACKLIST != "") {
    while ((getline line < BLACKLIST) > 0) {
      nf = split(line, b, "\t")            # not `n`: that name holds the count array
      if (nf < 3 || b[2] !~ /^[0-9]+$/) continue
      for (k = int(b[2] / W); k <= int(b[3] / W); k++) {
        key = b[1] SUBSEP k
        nbl[key]++
        bs[key, nbl[key]] = b[2] + 0
        be[key, nbl[key]] = b[3] + 0
      }
    }
    close(BLACKLIST)
  }
}

{
  bc = $4
  if (!(bc in keep)) next
  depth[bc]++

  mid = int(($2 + $3) / 2)
  k   = int(mid / W)
  key = $1 SUBSEP k
  if (key in nbl) {
    for (i = 1; i <= nbl[key]; i++)
      if (mid >= bs[key, i] && mid < be[key, i]) { masked++; next }
  }
  n[bc SUBSEP $1 SUBSEP k]++
}

END {
  print "barcode\tchr\tbin\tcount"
  for (k in n) {
    split(k, a, SUBSEP)
    if (depth[a[1]] < MINFRAGS) continue
    print a[1] "\t" a[2] "\t" a[3] "\t" n[k]
  }
  nd = 0; nk = 0
  for (bcx in depth) { nd++; if (depth[bcx] >= MINFRAGS) nk++ }   # `b` is the blacklist split array
  printf "[bin_fragments] %d masked by blacklist; %d/%d barcodes pass MINFRAGS=%d\n",
         masked, nk, nd, MINFRAGS > "/dev/stderr"
}
