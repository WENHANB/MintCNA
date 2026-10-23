# One pass over a fragments.tsv.gz to characterise scATAC count structure.
#
# Counting rule: ONE per fragment, placed by MIDPOINT. Column 5 is the number of
# PCR-duplicate read pairs that collapsed onto this fragment -- it is sequencing
# redundancy, not extra molecules -- so it is deliberately ignored.
#
# Emits three files:
#   <out>.cell5k / .cell100k / .cell500k / .cell1m  counts per bin for ONE cell
#   <out>.malig500k                                 malignant pseudobulk, 500 kb
#   <out>.ref500k                                   reference  pseudobulk, 500 kb
#
# Usage: zcat frags.gz | awk -F'\t' -v CELLS=... -v FOCUS=... -v OUT=... -f this.awk

BEGIN {
  while ((getline line < CELLS) > 0) {
    n = split(line, f, "\t")
    if (f[1] == "barcode") continue
    type[f[1]] = f[2]
    isref[f[1]] = (f[5] == "TRUE")
  }
  close(CELLS)
}

{
  bc = $4
  if (!(bc in type)) next
  mid = int(($2 + $3) / 2)

  if (bc == FOCUS) {
    c5[$1 "\t" int(mid/5000)]++
    c100[$1 "\t" int(mid/100000)]++
    c500[$1 "\t" int(mid/500000)]++
    c1m[$1 "\t" int(mid/1000000)]++
  }
  if (type[bc] == "Malignant") m500[$1 "\t" int(mid/500000)]++
  else if (isref[bc])          r500[$1 "\t" int(mid/500000)]++
}

END {
  for (k in c5)   print k "\t" c5[k]   > (OUT ".cell5k")
  for (k in c100) print k "\t" c100[k] > (OUT ".cell100k")
  for (k in c500) print k "\t" c500[k] > (OUT ".cell500k")
  for (k in c1m)  print k "\t" c1m[k]  > (OUT ".cell1m")
  for (k in m500) print k "\t" m500[k] > (OUT ".malig500k")
  for (k in r500) print k "\t" r500[k] > (OUT ".ref500k")
}
