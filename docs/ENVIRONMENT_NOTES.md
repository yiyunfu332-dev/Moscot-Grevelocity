# Reproducible environment

The verified project run used the `vae_gpu` Python 3.11 environment. Its current
state mixes packages from the Conda environment and `~/.local`, which caused the
earlier `anndata`/`mudata` and Numba cache failures. `environment.yml` records a
clean, single-environment specification for future reconstruction.

Before running the notebook in the existing environment, select the `vae_gpu`
kernel. The first notebook import cell now creates writable cache directories
and the moscot import cell applies the required `anndata 0.10.9` compatibility
shim before importing `mudata`.

Validate the current kernel from a terminal with:

```bash
conda run -n vae_gpu python validate_vae_gpu_environment.py
```

Create a clean environment later with:

```bash
conda env create -f environment.yml
conda activate zebrafish-growth-reproducible
python -m ipykernel install --user --name zebrafish-growth-reproducible
```

Do not install packages from inside the analysis notebook. Recreate the
environment and restart the kernel when versions need to change.
