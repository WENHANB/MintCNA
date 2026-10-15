# Reproducing the paper

Scripts that generated every result and figure in Bao, Qin & Xiao (MintCNA). They are the scripts
that were run, copied from the analysis tree; the only edits are that absolute paths were replaced by
placeholders and that paths to other scripts point inside this folder. Run
`bash reproduce/set_root.sh <project-root> <R-library>` once: it fills in the placeholders
(`<MINTCNA_ROOT>`, `<R_LIBS_USER>`, `<REPRODUCE_DIR>` = this folder) and seeds the figure-data
directories under `<project-root>/figure4paper/` from `data/`. The R functions the analysis scripts
source are in `lib/` (`lib/MultiCNA/`, the function set used for the simulation and CRC01 runs;
`lib/eval/`, `lib/benchmark/` for the multiome evaluation); the multiome runs use the package itself
(`../MintCNA/`, `../scripts/run_mintcna.R`). Tool versions are listed in `../environment.lock.yml` and
in the Methods.

Paths inside the scripts assume the layout `Simulation/cnasim/<setting>/`, `Data/CRC/`,
`Result/CRC/`, `ATAC_RNA/` under the project root. The simulated datasets are on Zenodo
(see Data availability); unpack each `<setting>.tar.gz` into `Simulation/cnasim/<setting>/`.

## 01_simulation — simulated scDNA + paired scRNA
1. `CNAsim.slurm` — CNAsim 1.3.5, `cnasim -m 1 -n 100 -n1 0.10 -c 3 -U -p1 6 -k 10000 -l1 10000 -l 5000000 -B 100000 -C <coverage>`
   for coverage 0.05, 0.1 and 0.2. CNAsim was run without a seed, so the generated read counts are
   deposited rather than regenerated.
2. `DNA_process.R` — bin and normalise read counts (100 kb); `DNA_process_n0.05.R` / `DNA_process_n0.1.R`
   add multiplicative Gaussian noise (σ = 0.05 / 0.1, `set.seed(1)`).
3. `simulate_rna_from_cnasim.R` — paired scRNA counts with Splatter 1.28.0 (`set.seed(1)`), copy-number
   effect injected at concordance p; `concordance/` holds the 25/50/75 % variants used in Supplementary
   Fig. 5a. `RNA_process.R`, `DNA_RNA_combine.R` build the combined input.

## 02_mintcna_runs — MintCNA on the simulations and on CRC01 (multi-omics mode)
`MintCNA.R` (DNA + RNA) with its SLURM wrappers; `CAE_SCAtt_DNAmeth.py` is the AG-DCAE denoiser as run
(latent_dim 250). `crc01/` — MintCNA on CRC01 (run of 2026-10-14): `01_ae_d250.py` (AG-DCAE on the 197 cells
with both modalities, latent_dim 250), `02_combine_d250.R` (denoised DNA + raw RNA log2 ratios),
`03_call_multi_d250.R` (change points, h = 15; calling at st = 3 and st = 5), `05_cluster_d250.R` (subclones by
the Supplementary Methods procedure), `06_compare_d250.R` (comparison with the earlier run); `s*.sbatch`
wrappers and `timing/` (the re-timed AG-DCAE stage behind Supplementary Fig. 9a). One cell (CRC01_ML3_085) is
excluded after calling by the QC rule in Methods (smoothed log2-ratio s.d. > 3× cohort median); every
downstream script applies it through `utils.R`.

## 03_benchmark_simulation — compared methods on the simulations
HMMcopy, SCOPE, SCCNV, SCONCE, rcCAE, SCICoNE and FLCNA: one `<tool>.R` (or `_input.R` / `_output.R`
pair for tools run outside R) per method, plus SLURM wrappers.

## 04_benchmark_crc — CRC01 (scTrio-seq2)
Benchmark runs (`HMMcopy.R`, `SCCNV.R` with `SCCNV_Fun.R`, `SCICoNE_*.R`, `rcCAE_output.R`, `FLCNA.R`;
MintCNA itself is under `02_mintcna_runs/crc01/`) and the evaluation against bulk WGS (`WGS_evaluation.R`,
`Pearson.R`, `RMSE.R`, `Boundary_Bias.R`, `Event_F1.R`, `Frequency.R`, `Heatmap.R`). In these scripts the
MintCNA result files carry the `_d250` suffix of the 2026-10-14 run (the only edit besides the placeholders).

## 05_benchmark_multiome — CARE glioma 10x Multiome
`prep/` (input matrices), `mintcna/` (configs per sample and run scripts), `benchmark/` (AtaCNA,
Copy-scAT, RIDDLER, epiAneufinder, inferCNV, CONGAS+), `eval/` (accuracy against bulk WGS).

## 06_evaluation — metrics and tables
`Evaulation_all.R` (event-level precision/recall/F1 with ≥ 20 % overlap, boundary bias, RMSE),
`Extract_tables.R` / `Merge_tables.R` (→ `data/simulation_merged_tables/`), `cluster_compare.R`
(subclone ARI, Supplementary Table 1), `Binsize_F1.R` (Supplementary Figs 2–4 inputs), `Concor_F1.R`
and `s5_cn_pr_3cov.R` (Supplementary Fig. 5), `TIme_evaluation.R` (CPU time).

## 07_figures
`shared/fig_palette.R` (colours, fonts, sizes used by every figure), `figure2-radar/`,
`figure3-boundary-rmse/`, `figure4-CRC/`, `figure5-multiome/`, `supplementary/` (Supplementary Figs
1–9; `pdf_text_sizes.py` checks fonts and text sizes in the exported PDFs). The source data each figure
script reads are in `data/figure_source_data/`. `figure4-CRC/utils.R` is the loader that Fig. 4 and
Supplementary Fig. 6 source (as `Result/CRC/figures/scripts/utils.R`): it selects the CRC01 run
(`CRC_SUFFIX`, default `_d250`) and applies the per-cell QC exclusion (`CRC_EXCLUDE`, default CRC01_ML3_085).

## data
`simulation_merged_tables/` — per-setting precision/recall/F1 and per-cell boundary bias and RMSE for
all methods (what Figs 2–3 and Supplementary Figs 2–4 are drawn from).
`figure_source_data/` — the per-figure CSVs.
