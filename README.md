# Waveguide Effective Index Prediction 
## Authors: Isabella Morgan & Ergun Simsek

Machine learning surrogate models for predicting the effective refractive index (`n_eff`) of Si₃N₄/SiO₂ rectangular waveguides, trained on FEM (finite element method) mode-solver data — as a fast alternative to running a full mode solver for every new waveguide geometry.

## Overview

Given a wavelength and waveguide cross-section (width, height), this project predicts the **TE-mode effective index** (`n_eff_TE`) using several regression approaches, and compares their accuracy:

| Model | Test R² | Notes |
|---|---|---|
| Linear Regression | Low | Baseline; `n_eff` is nonlinear in (λ, w, h), so this underfits |
| Fully Connected Neural Network (PyTorch) | ~0.9999+ | Best overall accuracy; scales well with more data |
| Support Vector Regression (RBF kernel) | High | Tuned via grid search over `C`, `gamma`, `epsilon` |

The trained surrogate can then be used to compute derived quantities — such as **group index** and **chromatic dispersion** `D(λ)` — via polynomial-smoothed differentiation, without needing to re-run the FEM solver at every wavelength.

## Data

- **Training/testing data** (`training_data.mat`, `testing_data.mat`) is generated from a MATLAB-based finite element mode solver, sweeping:
  - Wavelength: 0.75–1.65 μm
  - Waveguide width: 0.75–3.00 μm
  - Waveguide height: 0.60–1.00 μm
- Each row contains: `lambda_um`, `width_um`, `height_um`, `n_eff_TE`, `n_eff_TM`
- Material dispersion (Si₃N₄ core, SiO₂ cladding) is included via a Sellmeier-based refractive index lookup at each wavelength before the FEM solve.

## Workflow

1. **Data loading** — training/testing `.mat` files are pulled directly from this repo and loaded into pandas DataFrames.
2. **Baseline: Linear Regression** — establishes a lower bound on achievable accuracy.
3. **Fully Connected Neural Network (PyTorch)** — a `SiLU`-activated FCNN with input/output standardization (`StandardScaler`), early stopping on a held-out validation split, and `ReduceLROnPlateau` learning-rate scheduling. This is the primary surrogate model.
4. **Support Vector Regression (RBF kernel)** — an alternative surrogate, hyperparameters selected via grid search (`GridSearchCV`).

On Google CoLab, using a GPU, FCNN implementation is faster.
In terms of accuracy, there is not much difference between FCNN and SVM.


## Requirements
numpy
pandas
scipy
scikit-learn
torch
matplotlib
requests

## Usage

Run the notebook (`WG_n_eff_prediction.ipynb`) top to bottom, or adapt the individual model-training cells. Each model section reports **Test RMSE** and **Test R²**, plus a predicted-vs-true scatter plot against the `y = x` ideal line.
