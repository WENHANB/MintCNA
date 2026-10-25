# MintCNA demo

100 simulated cells (CNAsim, 3 clones), scDNA + scRNA, chr1–22. All inputs are included in `demo/input/`.

## Run

```bash
bash scripts/mintcna.sh --config demo/demo_config.yaml
```

About 10–12 minutes on 4 CPU cores. Results land in `demo/output/`.

DNA-only: in `demo/demo_config.yaml` set `mode: "dna"`, `rna_counts: null`, `gene_annotation: null`.

## Inputs (`demo/input/`)

| File | Content |
|---|---|
| `DNA_RC.RData` | `Y` (bins × cells raw counts), `gr` (GRanges with `gc`, `mapp`); `DNA_RC.csv` + `DNA_ref.csv` hold the same data as CSV |
| `rna_counts.tsv` | genes × cells raw counts |
| `genes.tsv` | `gene_id, chr, start, end` |
| `cell_types.tsv` | cell ID, `normal` |

`build_demo.R` regenerates these from the `c3_N100_100kb_l1e4_c0.2_n0.1` setting of the simulated data on Zenodo (doi 10.5281/zenodo.23385650); you do not need to run it.

## Outputs (`demo/output/`)

`mintcna_results.RData`: `cna_cn_matrix` (bins × cells copy number, 2 = diploid), `cna_df_all` (called segments), `subclones` (cell → cluster).

Intermediate files: `l2r_dna.csv`, `l2r_rna.csv`, `ref_dna.csv`, `ref_rna.csv` (preprocessing); `ae_dna.csv`, `latent_z.csv`, `loss_curve.pdf` (denoising); `combined.RData` (segmentation and calling).
