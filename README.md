# Causal Inference in High-Dimensional Observational Data
### A Comparative Simulation Study of TMLE and C-TMLE Estimators

[![R](https://img.shields.io/badge/R-%3E%3D4.0-blue)](https://www.r-project.org/)
[![Methods](https://img.shields.io/badge/Methods-TMLE%20%7C%20C--TMLE%20%7C%20IPW-green)]()
[![Status](https://img.shields.io/badge/Status-Active%20Research-orange)]()

---

## Overview

This repository contains simulation studies comparing causal inference estimators for the **Average Treatment Effect (ATE)** under observational data settings — where treatment is not randomly assigned and confounding must be carefully accounted for.

The work is motivated by real challenges in healthcare and genomic research: how do we reliably estimate the causal effect of a treatment or exposure when we cannot run a randomised controlled trial, and when the data is high-dimensional and messy?

This simulation study benchmarks four estimators across scenarios of increasing complexity, from low-dimensional correlated covariates to high-dimensional sparse settings with up to 100 covariates.

> **Context:** This work forms part of my PhD research in Statistics at the University of Edinburgh, where I apply these methods to large-scale observational healthcare and genomic data (UK Biobank). The simulations here use fully synthetic data — no real patient data is used or shared.

---

## Background: Why This Problem Matters

In observational studies — clinical registries, electronic health records, genomic cohorts — we observe who received a treatment and what happened to them, but we cannot control who got treated. Patients who receive a particular drug may systematically differ from those who don't, in ways that also affect their outcomes. This is **confounding**, and naively comparing treated and untreated groups gives a biased estimate of the true causal effect.

**Causal inference methods** address this by explicitly modelling:
1. **The outcome mechanism** Q(A, W): how the outcome Y depends on treatment A and covariates W
2. **The treatment mechanism** g(W): the probability of receiving treatment given covariates (the propensity score)

Getting both right — especially in high-dimensional settings — is the central challenge this project addresses.

---

## Methods Compared

| Estimator | Description | Key Property |
|-----------|-------------|--------------|
| **TMLE** | Targeted Maximum Likelihood Estimation | Doubly robust; semiparametrically efficient |
| **C-TMLE (Greedy)** | Collaborative TMLE with greedy covariate selection | Jointly optimises Q and g models |
| **C-TMLE1 (Lasso)** | Collaborative TMLE with Lasso regularisation path | Uses penalised regression for g selection |
| **C-TMLE0** | Collaborative TMLE with gradient-based selection | Optimises along regularisation gradient |

All estimators use **influence function-based inference** for standard errors and confidence intervals, which provides valid uncertainty quantification without relying on parametric assumptions.

In high-dimensional settings (`glmnet_update.R`), outcome and propensity score models are estimated using **cross-validated LASSO/elastic net** (`glmnet`), replacing standard GLMs to handle the curse of dimensionality.

**Key metrics evaluated across Monte Carlo replications:**
- Mean bias
- Empirical variance and mean squared error (MSE)
- 95% confidence interval coverage
- Bias-to-standard-error ratio

---

## Simulation Designs

Three data generating processes (DGPs) of increasing complexity:

### Study 1 — Low-dimensional, Correlated Continuous Covariates
Two correlated Gaussian covariates (ρ = 0.5). Three variants:
- **1a:** Nonlinear outcome mechanism `Q₀ = 1 + A − 0.7W₁ + 0.3·exp(−W₁W₂)`
- **1b:** Interaction in treatment mechanism `g₀ = expit(0.5 − 1.5W₁W₂)`
- **1c:** Randomised treatment (benchmark — all estimators should perform well)

### Study 2 — Binary Covariates with Induced Correlation
Eight binary covariates with a structured dependency chain (W4 depends on W1, W5 depends on W1–W4, etc.), mimicking the kind of correlated binary data common in clinical settings.

### Study 3 — High-Dimensional Sparse Setting
- **p = 100 covariates**, Toeplitz covariance structure with ρ = 0.5–0.9
- Only **k = 10–20 covariates** are truly predictive (sparse signal)
- Separate sparse signals for outcome and treatment mechanisms
- True ATE = 2

This is the most challenging setting and the primary focus of the `glmnet`-based estimators.

---

## Repository Structure

```
causal-inference-simulation/
│
├── README.md
│
├── scripts/
│   ├── 01_data_generating_processes.R   # DGP functions for all simulation studies
│   ├── 02_utilities.R                   # Helper functions, MC evaluation, plotting
│   ├── 03_estimators.R                  # Core estimators: DM, OLS, IPW, TMLE, C-TMLE
│   └── 04_estimators_highdim.R          # glmnet-based estimators for high-dim settings
│
├── notebooks/
│   └── simulation_study.Rmd             # Main analysis notebook (rendered below)
│
└── results/
    └── ctmle_glmnet.csv                 # Saved Monte Carlo results (k=100, n=5000)
```

---

## Key Results

Monte Carlo results based on **k = 100 replications**, **n = 5,000** observations, high-dimensional sparse setting (Study 3):

All four estimators recover the true ATE = 2. The C-TMLE variants demonstrate improved **bias-variance trade-off** compared to standard TMLE in the high-dimensional sparse setting, consistent with the theoretical motivation for collaborative estimation. Full results including coverage rates and CI distributions are in `notebooks/simulation_study.Rmd`.

---

## How to Run

### Prerequisites

```r
install.packages(c(
  "tidyverse", "MASS", "survey", "boot",
  "ggplot2", "skimr", "glm2",
  "glmnet", "randomForest", "SuperLearner"
))
```

### Reproducing the simulation

```r
# 1. Clone the repository
# 2. Open simulation_study.Rmd in RStudio
# 3. Set your working directory to the repo root
# 4. Knit or run chunks sequentially

# To run a quick single-dataset test:
source("scripts/01_data_generating_processes.R")
source("scripts/02_utilities.R")
source("scripts/04_estimators_highdim.R")

data_test <- sim3(n = 1000, p = 100)
estimate_TMLE_glmnet(data_test, true_value = 2)
```

> **Note:** The full Monte Carlo study (`k = 100`, `n = 5000`) is computationally intensive (~180-300 min). Pre-computed results are saved in `results/ctmle_glmnet.csv` for immediate exploration.

---

## Technical Stack

- **Language:** R (≥ 4.0)
- **Core packages:** `glmnet`, `SuperLearner`, `survey`, `glm2`, `tidyverse`
- **Estimation:** LASSO / elastic net via cross-validated `glmnet` for high-dimensional nuisance models
- **Inference:** Influence function-based standard errors (semiparametric efficiency theory)
- **Evaluation:** Monte Carlo simulation with empirical bias, variance, MSE, and coverage

---

## About This Project

This repository is part of my PhD research at the **University of Edinburgh, School of Mathematics**. My thesis focuses on estimating causal effects from observational healthcare data in the presence of non-ignorable missingness, applying C-TMLE methodology to identify genetic variants with causal effects on disease outcomes using **UK Biobank** data.

The simulation work here serves to validate and compare the performance of these estimators under controlled conditions before applying them to real data.

**Related areas of application:**
- Pharmacoepidemiology and clinical trials emulation
- Genomics / Mendelian randomisation
- Health technology assessment
- Any domain where RCTs are infeasible and high-quality causal estimates from observational data are needed

---

## Contact

**Juliet Asantewaa Sarpong**  
PhD Candidate, Statistics — University of Edinburgh  
📧 asantewaahsarpong@gmail.com  
🐙 [github.com/Asantewaah](https://github.com/Asantewaah)  
💼 [linkedin.com/in/asantewaah-sarpong](https://www.linkedin.com/in/asantewaah-sarpong/)

---

*Feedback, questions, and collaboration enquiries welcome.*
