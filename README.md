# OrdCRTPwr

R code accompanying the manuscript:

**Power and Sample Size Calculation for Cluster Randomized Trials with Ordinal Outcomes**

Lan Shi, Sarah Osmundson, Dandan Liu, and Bryan Blette

## Overview

This repository contains the R code used for the simulation studies, power and sample size calculations, and generation of figures and tables for the manuscript *Power and Sample Size Calculation for Cluster Randomized Trials with Ordinal Outcomes*.

The manuscript develops closed-form power and sample size formulas for cluster randomized trials (CRTs) with ordinal outcomes under a generalized estimating equations (GEE) framework. The proposed methods consider both independence and exchangeable working correlation structures and explicitly account for the multidimensional within-cluster correlation structure induced by ordinal outcomes.

The repository includes code for:

- three-category ordinal outcome simulations;
- four-category ordinal outcome simulations with design-stage category collapsing;
- GEE-based analytic power calculations;
- Whitehead design-effect power calculations;
- binary-collapse power calculations;
- superpopulation-based estimation of within-cluster correlation parameters;
- empirical power simulations;
- the chronic hypertension in pregnancy (CHTN) application; and
- generation of the manuscript and supplementary figures and tables.

## Repository structure

```text
OrdCRTPwr/
├── R/
│   ├── functions/
│   ├── 3cat/
│   ├── 4cat/
│   └── application/
│
├── results/
│   ├── 3cat/
│   └── 4cat/
│
├── manuscript_tables_and_figures/
│   └── Supplementary/
│
└── OrdCRTPwr.Rproj
```

### `R/functions/`

Contains functions used throughout the analyses, including routines for:

- GEE power calculations for three-category ordinal outcomes;
- binary-outcome power calculations;
- Whitehead power calculations;
- bridge-distribution calculations;
- marginal ordinal probabilities;
- within-cluster correlation calculations;
- category collapsing;
- variance calculations; and
- solving for the required number of clusters.

### `R/3cat/`

Contains the simulation workflow for three-category ordinal outcomes.

The main workflow is organized approximately as follows:

```text
00_create_superpop_scenarios.R
01_use_superpop_obtain_rhomat_array.R
02_use_superpop_obtain_rankICC_exploratory.R
03_create_empirical_scenarios.R
04_obtain_emp_power_mygee.R
05_analytic_power_estimators_combine_empirical.R
```

Additional scripts generate the corresponding manuscript and supplementary figures and tables.

### `R/4cat/`

Contains the simulation workflow for four-category ordinal outcomes. For design-stage power calculation, adjacent categories are collapsed to obtain a three-category outcome before applying the proposed GEE-Ordinal formula.

The folder contains scripts for constructing simulation scenarios, estimating superpopulation correlation parameters, calculating empirical and analytic power, and generating manuscript and supplementary results.

### `R/application/`

Contains code for the chronic hypertension in pregnancy (CHTN) application and generation of **Figure 3**.

The original patient-level data used for this application are not included in this repository because of privacy restrictions. The code and derived quantities needed to reproduce the reported calculations are provided where applicable.

### `manuscript_tables_and_figures/`

Contains the final figures and tables reported in the manuscript.

The `Supplementary/` subdirectory contains supplementary figures and tables.

## Simulation study

The simulation study considers balanced two-arm CRTs with ordinal outcomes generated from a mixed-effects proportional odds model with a bridge-distributed cluster-level random effect.

The primary simulation settings include:

- number of clusters: \(N = 20\) or \(50\);
- cluster size: \(J = 50\);
- latent ICC: \(0.01\), \(0.02\), or \(0.05\);
- conditional treatment effect: \(\theta_c = \log(1.5)\);
- multiple marginal ordinal outcome distributions; and
- 10,000 simulated datasets per scenario for estimation of empirical power.

The methods compared include:

1. **GEE-Ordinal** — the proposed GEE-based power formula using the within-cluster correlation matrix;
2. **GEE-Binary1** — a binary collapse of the ordinal outcome;
3. **GEE-Binary2** — an alternative binary collapse; and
4. **WH-DE-Latent** — Whitehead's ordinal power formula combined with a design effect based on the latent ICC.

For the four-category simulations, empirical power is calculated using the full four-category outcome, while the proposed analytic approach is evaluated using a design-stage collapse to three adjacent categories.

## Reproducing the analyses

Clone the repository:

```bash
git clone https://github.com/Slan1997/OrdCRTPwr.git
cd OrdCRTPwr
```

Open `OrdCRTPwr.Rproj` in RStudio.

The scripts in `R/3cat/` and `R/4cat/` are numbered to indicate the general order of the simulation workflow. Because some simulations are computationally intensive, particularly the empirical power simulations, runtime can be substantial.

Final manuscript figures and tables are generated by the corresponding figure/table scripts and are stored in:

```text
manuscript_tables_and_figures/
```

and

```text
manuscript_tables_and_figures/Supplementary/
```

## Main manuscript outputs

The repository contains code used to generate:

- **Figure 1:** Estimated versus empirical power for three-category ordinal outcomes with 20 clusters.
- **Figure 2:** Estimated versus empirical power for four-category ordinal outcomes with 20 clusters.
- **Figure 3:** Required number of clusters for the CHTN application under assumed and data-driven within-cluster correlation specifications.

Supplementary scripts and output files reproduce the corresponding supplementary correlation summaries and power results.

## Data availability

The simulation code and code used to generate the manuscript figures and tables are publicly available in this repository.

The original datasets used for the CHTN application are **not publicly available due to privacy restrictions**. No identifiable or individual-level clinical data are included in this repository.

## Citation

If you use the methods or code in this repository, please cite:

> Shi L, Osmundson S, Liu D, Blette B. *Power and Sample Size Calculation for Cluster Randomized Trials with Ordinal Outcomes.*

Citation information will be updated upon publication.

## Authors

- **Lan Shi** — Department of Biostatistics, Vanderbilt University Medical Center
- **Sarah Osmundson** — Department of Obstetrics and Gynecology, Vanderbilt University Medical Center
- **Dandan Liu** — Department of Biostatistics, Vanderbilt University Medical Center
- **Bryan Blette** — Department of Biostatistics, Vanderbilt University Medical Center

## Contact

For questions about the manuscript or code, please contact:

**Bryan Blette**  
Department of Biostatistics  
Vanderbilt University Medical Center  
bryan.blette@vumc.org
