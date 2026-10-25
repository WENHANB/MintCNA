# Reproducing the paper

Two parts. `figures/` draws every figure of the paper from the source data shipped here; `analysis/`
holds the scripts that produced that source data from the raw data, one script per method, organised by
dataset. Everything is run from the root of the repository.

```
reproduce/
├── figures/          source data -> figures      (self-contained; data/ is included)
└── analysis/         data -> source data         (needs the datasets and the compared tools)
    ├── simulation/   CNAsim + Splatter simulation, MintCNA, six scDNA-seq methods, evaluation
    ├── crc01/        CRC01 scTrio-seq2: input matrices, MintCNA, bulk WGS truth
    └── multiome/     CARE glioma 10x Multiome: download, input, MintCNA, six scATAC/scRNA methods, evaluation
```

## Part 1 — figures from source data

```bash
Rscript reproduce/figures/fig2_radar.R
Rscript reproduce/figures/fig3_boundary_rmse.R
Rscript reproduce/figures/fig4_crc/assemble_fig4.R
Rscript reproduce/figures/fig5_multiome/assemble_fig5.R
Rscript reproduce/figures/supp/figS1_denoising.R            # one script per Supplementary Figure, S1-S9
```

Output lands in `reproduce/figures/output/`. The R packages are those of `environment.yml` plus
data.table, ggplot2, cowplot, ggh4x, ggrepel, patchwork, ComplexHeatmap, circlize, svglite, ragg, pdftools.
`shared/fig_palette.R` holds the colours, fonts (Arial) and sizes shared by all figures.

| Figure | Script | Reads from `figures/data/` | Produced by (`analysis/`) |
|---|---|---|---|
| Fig. 2 | `fig2_radar.R` | `simulation_tables/summary_{one,del,dup}.csv` | `simulation/evaluate.R`, `merge_tables.R` |
| Fig. 3 | `fig3_boundary_rmse.R` | `simulation_tables/{bb,RMSE}_percell.csv` | same |
| Fig. 4 | `fig4_crc/assemble_fig4.R` (panels a–c) | `fig4_crc/*.RData` (CRC01 calls of every method, bulk WGS truth `DNAcopy_qc.RData`), `latent_z_dna.csv` | `crc01/run_mintcna.sh`, `simulation/methods/*` on the CRC01 matrix, `crc01/wgs_truth*.R` |
| Fig. 5 | `fig5_multiome/assemble_fig5.R` (panels a–f) | `fig5_multiome/pearson_vs_truth_own.tsv`, `f1_states.tsv`, `P111T2_results_h15a005e020.RData`, `clone_embedding_P111T2*.tsv`, `tf_motif_P111T2.tsv`, `atac_rna_clone_P111T2.tsv`, `fig5e_rna_de_P111T2.tsv`, `cna_dosage_link_P111T2_atac.tsv` | `multiome/run_mintcna.sh`, `multiome/methods/*`, `multiome/evaluation/*` |
| Supp. Fig. 1 | `supp/figS1_denoising.R` | `supp/figS1/*` (raw, truth, CAE + MSE and AG-DCAE profiles of the two plotted cells) | `simulation/denoising_comparison/*` |
| Supp. Figs 2–4 | `supp/figS2-4_f1_by_length.R` | `simulation_tables/by_setting/summary_mss_*.csv` | `simulation/evaluate.R`, `merge_tables.R` |
| Supp. Fig. 5 | `supp/figS5_rna_concordance_cn.R` | `supp/figS5a_concordance_f1_published.csv`, `supp/s5_cn_pr_3cov.csv` | `simulation/cn_state_precision_recall.R` (b–d; panel a plots the values of the published figure) |
| Supp. Fig. 6 | `supp/figS6_chr17_arm_loss.R` | `fig4_crc/MultiCNA_AE_197_d250.RData` | `crc01/run_mintcna.sh` |
| Supp. Fig. 7 | `supp/figS7_cohort_flow.R` | `fig5_multiome/sample_manifest.csv` | — |
| Supp. Fig. 8 | `supp/figS8_subclone_biology.R` | as Fig. 5 | as Fig. 5 |
| Supp. Fig. 9 | `supp/figS9_cpu_time.R` | `supp/figS9_dna_cpu_time.csv`, `supp/figS9_atac_runtime_sacct.tsv` | `multiome/evaluation/runtime_table.R` (ATAC); SLURM accounting of the simulation and CRC01 runs (DNA) |

`figures/data/` also carries the per-figure tables of the paper's Source Data (`fig2_radar_metrics.csv`,
`fig3*_percell.csv`, `fig4c_track_metrics.csv`, `fig5a_*.csv`, `fig5b_*.csv`, `fig5c_*.csv`, `fig5e_*`, `fig5f_*`,
`figS*.csv`); the scripts regenerate them into `figures/output/`.

## Part 2 — source data from the data

Set `MINTCNA_DATA` to a directory laid out as `Simulation/cnasim/<setting>/`, `Data/CRC/`, `Result/CRC/` and
`ATAC_RNA/` (the layout every script assumes; each script's header states what it reads and writes), and run
from the repository root with the `mintcna` conda environment active. R packages beyond that environment are
listed with each dataset below; tool versions are in the Methods of the paper and in `environment.lock.yml`.

### Data

| Dataset | Where | Used by |
|---|---|---|
| Simulated scDNA + scRNA, 9 settings | Zenodo, doi 10.5281/zenodo.23385650 — unpack each archive into `Simulation/cnasim/<setting>/`; the archives hold the CNAsim read counts and truth, the Splatter RNA counts, the processed inputs and the paper's MintCNA calls | simulation, Figs 2–3, Supp. Figs 1–5, Supp. Tables 1–3 |
| CRC01 scTrio-seq2 (Bian et al. 2018) | GEO GSE97693 (single cells); bulk WGS under controlled access (EGA) | Fig. 4, Supp. Fig. 6 |
| CARE IDH-mutant glioma 10x Multiome (Johnson et al. 2026) | GEO GSE327580 (snATAC fragments), GSE326221 (snRNA matrices); bulk WGS segments and cell annotations from Synapse syn63600501 | Fig. 5, Supp. Figs 7–9 |

### Simulation (`analysis/simulation/`)

```
cnasim.sh                                   CNAsim 1.3.5, 100 cells, 3 clones, 100 kb bins, 0.05x / 0.1x / 0.2x
simulate_rna_from_cnasim.R <setting>        paired scRNA counts (Splatter), concordance 25-100 %
make_mintcna_input.R <setting> [sigma]      DNA_RC.RData (counts + GC / mappability, CODEX2) and mintcna_input/;
                                            sigma 0.05 / 0.1 = the _n0.05 / _n0.1 noise settings
run_mintcna.sh <setting>                    MintCNA DNA + RNA and DNA only  (config_mintcna.yaml)
methods/<tool>.R | <tool>_run.sh <setting>  one script per compared method (table below)
evaluate.R <setting>                        precision / recall / F1 by direction and event length, boundary bias, RMSE
merge_tables.R                              the nine settings -> the tables of Figs 2-3, Supp. Figs 2-4 and Supp. Tables 1-3
cn_state_precision_recall.R                 per-state precision / recall (Supp. Fig. 5b-d)
denoising_comparison/                       build_truth.py, run_denoisers.py, cae_mse.py (Supp. Fig. 1)
```

CNAsim was run without a seed, so the read counts are not regenerated exactly: start from the deposited
`readcounts.tsv`, `profiles.tsv` and `cell_types.tsv` of each setting.

| Tool (version) | Script |
|---|---|
| HMMcopy 1.46.0 | `methods/HMMcopy.R <setting>` |
| SCOPE 1.16.0 | `methods/SCOPE.R <setting>` |
| SCCNV (R re-implementation of the published algorithm) | `methods/SCCNV.R <setting>` (+ `SCCNV_functions.R`) |
| SCONCE (commit 1f1e1f9) | `methods/SCONCE_run.sh <setting>` = `SCONCE_input.R` → sconce → `SCONCE_output.R`; needs `SCONCE_REPO` |
| rcCAE (zhyu-lab/rccae, 2026-03) | `methods/rcCAE_run.sh <setting>` = `rcCAE_input.R` → rccae → `rcCAE_output.R`; needs `RCCAE_REPO`, MATLAB |
| SCICoNE 1.0.0 | `methods/SCICoNE_run.sh <setting>` = `SCICoNE_input.R` → SCICoNE → `SCICoNE_output.R`; needs `SCICONE_BIN`, `SCICONE_REPO` |

R packages: CODEX2, BSgenome.Hsapiens.UCSC.hg38, splatter, HMMcopy, SCOPE, GenomicRanges, data.table, dplyr.
The denoising comparison runs in the `mintcna` environment (PyTorch).

### CRC01 (`analysis/crc01/`)

```
bin_counts.R                 single-cell methylome reads -> 100 kb read counts with GC / mappability (DNA_RC.RData, cell_types.tsv)
rna_counts.R                 RSEM expected counts -> rna_counts.tsv, genes.tsv (GENCODE v19)
run_mintcna.sh               MintCNA DNA + RNA and DNA only  (config_mintcna.yaml, hg19)
wgs_truth.R                  bulk WGS of the four tumour sites -> log2 ratio to normal blood -> CBS -> states
wgs_truth_postqc.R           events shorter than 6 Mb set to neutral -> DNAcopy_qc.RData, the truth of Fig. 4
```

The six DNA methods were run on CRC01 with the simulation scripts applied to `Data/CRC/mintcna_input/DNA_RC.RData`
(the same `Y` + `gr` layout). The accuracy metrics of Fig. 4 are computed inside `figures/fig4_crc/panel_b_accuracy.R`
and `panel_c_tracks.R`. R packages: CODEX2, BSgenome.Hsapiens.UCSC.hg19, DNAcopy, modSaRa2, GenomicRanges, data.table.

### Multiome (`analysis/multiome/`)

```
download_geo.sh {gex|atac|all}       snATAC fragments (GSE327580) and snRNA matrices (GSE326221)
download_synapse.py                  bulk-DNA copy-number segments and the authors' cell annotations (Synapse)
export_cell_lists.R                  the authors' QC-passing nuclei, cell types and reference (diploid) nuclei
prep_fragments.sh                    fragments restricted to those nuclei
run_mintcna.sh <care_id>             bin the fragments (scripts/atac_bin_fragments.awk), build the ATAC input
                                     (scripts/make_atac_input.R) and the RNA input (make_rna_input.R), run MintCNA
                                     ATAC + RNA (config_multi.yaml) and ATAC only (config_atac.yaml)
methods/run_<tool>.R <care_id>       one script per compared method (table below)
evaluation/prepare_truth.R           GATK segments -> 100 kb hg38 truth
evaluation/mintcna_profile.R         MintCNA pseudobulk profile on the 500 kb evaluation grid
evaluation/pearson_vs_truth.R        pseudobulk Pearson of every method (Fig. 5a); pearson_helpers.R loads the profiles
evaluation/f1_states.R               bin-level precision / recall / F1 (Fig. 5b)
evaluation/clone_in_embedding.R      the subclones in the RNA embedding (Fig. 5d)
evaluation/tf_motif_subclone.R       motif accessibility between the subclones, chromVAR (Supp. Fig. 8a)
evaluation/rna_de_clones.R           differential expression between the subclones (Fig. 5e)
evaluation/atac_to_rna_clone.R       promoter accessibility vs expression between the subclones (Fig. 5f)
evaluation/cna_dosage_link.R         copy-number dosage between the subclones (Supp. Fig. 8b,c)
evaluation/runtime_table.R           CPU time and memory from the SLURM accounting (Supp. Fig. 9)
```

The 20 samples are P59T2 P59T3 P100T1 P100T3 P103T1 P103T2 P104T1 P104T2 P105T1 P105T3 P106T1 P106T2 P108T1
P108T2 P109T1 P109T2 P110T1 P110T2 P111T2 P111T3.

| Tool (version) | Script |
|---|---|
| AtaCNA 0.1.0 (R package AtaCNV) | `methods/run_atacnv.R <care_id>` |
| Copy-scAT 0.40 | `methods/run_copyscat.R <care_id>` |
| RIDDLER (commit a027120) | `methods/run_riddler.R <care_id>` (+ `fast_window_matrix_dt.R`) |
| epiAneufinder 1.1.5 | `methods/run_epiAneufinder.R <care_id> 1e5 20000 w100kb_mf20000` |
| inferCNV 1.26.0 (RNA arm) | `methods/run_infercnv.R <care_id>`, then `export_infercnv.R` |
| CONGAS+ (Rcongas 0.2.0, congas 0.0.77) | `methods/run_congas.R <care_id>` |

R packages: GenomicRanges, rtracklayer, Matrix, Seurat, ensembldb + EnsDb.Hsapiens.v86, chromVAR, motifmatchr,
JASPAR2020, TFBSTools, FNN, cluster, data.table, and the packages of the compared methods.
