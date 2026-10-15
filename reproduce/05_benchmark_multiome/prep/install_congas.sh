#!/bin/bash
# CONGAS+ needs its own conda env: the published requirements.txt pins
# torch==1.7.1 / pyro-ppl==1.5.2 / numpy==1.18.1 / pandas==1.0.3, which will not
# build against a modern HPC stack. We relax the pins to the oldest versions that
# still install cleanly and check that the API we need still works.
# R's rcongas drives the Python package through reticulate, so the R session must
# be pointed at this env (see Code/benchmark/run_congas.R).
set -eu
ENV=<MINTCNA_ROOT>/ATAC_RNA/Tools/congas_env
module load conda
conda create -y -p "$ENV" python=3.10 pip
conda run -p "$ENV" pip install --quiet "torch>=2.0" "pyro-ppl>=1.8" numpy pandas scikit-learn matplotlib seaborn
# The Python backend is caravagnalab/CONGASp, NOT caravagnalab/congas. The latter
# installs as congas_old-0.0.2 -- the original single-modality CONGAS, whose
# congas.models has MixtureGaussian/HmmSegmenter and no LatentCategorical. Rcongas
# 0.2.0 calls cg_mod$LatentCategorical and cg$Interface, so with the wrong package
# every fit dies inside easypar with "No module named congas.models.HmmSimple" and
# the only symptom that reaches the user is "4/4 computations returned errors".
# The correct repo is named in Rcongas' own DESCRIPTION.
conda run -p "$ENV" pip install --quiet "git+https://github.com/caravagnalab/CONGASp.git"
conda run -p "$ENV" python -c "import congas, torch, pyro; print('congas', getattr(congas,'__version__','?'), '| torch', torch.__version__, '| pyro', pyro.__version__)"
