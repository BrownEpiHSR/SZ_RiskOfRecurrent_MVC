# SZ_RiskOfRecurrent_MVC
Sirui Zhang et al. - Risk of a Recurrent Motor Vehicle Crash Following Antihypertensive Medication Changes Among Older Drivers

# Risk of a Subsequent Motor Vehicle Crash Following Antihypertensive 
# Medication Changes Among Older Drivers

**Sirui Zhang, MPH; Nina Joyce, PhD; Adam D'Amico, MPH; Arman Oganisian, PhD;
Andrew R. Zullo, PharmD, PhD; Daniel A. Harris, PhD; Kaleen N. Hayes, PharmD, PhD**

Department of Epidemiology, Brown University School of Public Health, Providence, RI, USA
Center for Gerontology and Healthcare Research, Brown University School of Public Health,
Providence, RI, USA
Department of Biostatistics, Brown University School of Public Health, Providence, RI, USA 
Department of Health Services, Policy, and Practice, Brown University School of Public
Health, Providence, RI, USA
Department of Epidemiology, College of Health Sciences, University of Delaware, Newark, DE, USA


---

## Description

This repository contains data documentation and analytic code supporting
the manuscript titled *Risk of a Subsequent Motor Vehicle Crash Following
Antihypertensive Medication Changes Among Older Drivers*, submitted to the
*Journal of the American Geriatrics Society*.

The analysis applies a Clone-Censor-Weight (CCW) target trial emulation
framework to linked New Jersey police crash records and Medicare claims
(Parts A, B, and D) to estimate the causal effect of antihypertensive
treatment strategy changes (intensification, deprescribing, no change)
on 12-month recurrent motor vehicle crash risk among 54,658 older
Medicare-enrolled drivers in New Jersey, 2008–2017.

**No individual-level data, Protected Health Information (PHI), or
Personally Identifiable Information (PII) are included in this
repository.** Data were obtained under CMS Data Use Agreement
RSCH-2021-56434 and cannot be shared publicly. Researchers interested
in accessing CMS data should visit
[ResDAC](https://resdac.org) to get started.

---

## Repository Contents

| File | Type | Description |
|---|---|---|
| `CCW_template.R` | R script | Generalized, fully annotated CCW analytic template |
| `RiskOfRecurrent_MVC.xlsx` | Excel workbook | Data documentation: project overview, cohort development, variable descriptions, and antihypertensive drug class lists |
| `README.md` | Markdown | This file |

---

## Funding

This research was supported in full by the National Institute on Aging
under award RF1AG087210, totaling $2,612,229, with no amount financed
from non-federal sources. The content is solely the responsibility of
the authors and does not necessarily represent the official views of
the National Institutes of Health.

Additional support was provided by the National Institute on Aging
under award GR5271876 (PI: Nina Joyce), *The Risks and Consequences of
a Motor Vehicle Crash in Older Adults with Alzheimer's Disease and
Related Dementias*.

---

## License

This repository is shared under the MIT License. See `LICENSE` for
full terms. You are free to use, modify, and distribute this code
with attribution.

---
> **Repository DOI:** [To be issued via Zenodo upon manuscript submission]
> **Manuscript Status:** Under preparation for submission to the *Journal of the American Geriatrics Society* (JAGS)

---

## Table of Contents
- [Repository Overview](#repository-overview)
- [Study Background](#study-background)
- [Data Sources](#data-sources)
- [Repository Structure](#repository-structure)
- [Code Execution](#code-execution)
- [Methods Summary](#methods-summary)
- [Key Variables](#key-variables)
- [Software and Package Requirements](#software-and-package-requirements)
- [Data Documentation](#data-documentation)
- [Reproducibility Notes](#reproducibility-notes)
- [License](#license)
- [Citation](#citation)
- [Contact](#contact)

---

## Repository Overview

This repository contains the analytic code supporting the manuscript examining the **causal effect of antihypertensive treatment strategy changes** (Intensification, Deprescribing, No Change) on 12-month risk of recurrent motor vehicle crash (MVC) among older Medicare-enrolled drivers.

The analysis applies a **Clone-Censor-Weight (CCW) target trial emulation** framework to a linked New Jersey crash records and Medicare claims dataset (Parts A, B, D). The analytic cohort consists of **54,658 unique older drivers** (67,845 person-crashes) involved in a police-reported index MVC in New Jersey between 2008 and 2017.

**No individual-level data, Protected Health Information (PHI), or Personally Identifiable Information (PII) are included in this repository.**

---

## Study Background

Nearly 20% of older adult drivers who experience a motor vehicle crash (MVC) will have a second MVC. Antihypertensives are prevalent potentially driver-impairing medications and may affect subsequent MVC risk in this population. We examined whether changes in antihypertensive medications affect the 12-month risk of recurrent crashes in older drivers.

Using a CCW framework, we emulated a hypothetical randomized trial in which eligible older drivers were assigned at the time of their index crash to one of three sustained antihypertensive management strategies and followed for 12 months. This approach addresses key methodological limitations of prior work, including immortal time bias and confounding by indication.

The full target trial protocol is described in Supplementary Table S1 of the manuscript.

---

## Data Sources

| Dataset | Full Name | Coverage | Purpose |
|---|---|---|---|
| NJ-SHO | New Jersey Safety and Health Outcomes records | 2008-2017 | Index and recurrent MVC ascertainment; at-fault status |
| MBSF | Medicare Beneficiary Summary File | 2007-2017 | Demographics; enrollment; mortality dates |
| MedPAR | Medicare Provider Analysis and Review | 2007-2017 | Inpatient hospitalizations |
| Carrier (Part B) | Medicare Carrier claims | 2007-2017 | Outpatient visits |
| Part D | Medicare Prescription Drug Event records | 2007-2017 | Antihypertensive dispensing; PDI drug class ascertainment |
| MDS/SNF | Minimum Data Set / Skilled Nursing Facility claims | 2007-2017 | Long-term care entry dates |

**Data Access:** These data were obtained under a Data Use Agreement (DUA) with the Centers for Medicare and Medicaid Services (CMS). The linked NJ-SHO-Medicare dataset is not publicly available. Only de-identified summary statistics and analytic code are included in this repository.

---

## Repository Structure

```
Zhang_Antihypertensives_MVC/
|
|-- README.md                          <- This file
|-- LICENSE                            <- MIT License
|-- CITATION.cff                       <- Machine-readable citation (added after publication)
|
|-- data_documentation/
|   |-- RiskOfRecurrent_MVC.xlsx       <- Project overview, cohort development,
|                                         variable descriptions, drug class lists
|
|-- code/
|   |-- CCW_template.R                 <- Generalized, fully annotated CCW template
|                                         with placeholder variable names for replication
```

> **Note on code availability:** This repository provides a generalized analytic template (`CCW_template.R`) that documents the full CCW implementation pipeline with annotated placeholder variable names. The template is designed to enable methodological replication of this study design with analogous datasets. Project-specific analysis scripts are available from the corresponding author upon reasonable request.

---

## Code Execution

The `CCW_template.R` file is organized into sequential sections covering the full analytic pipeline:

| Section | Description |
|---|---|
| 0 | Package loading and global configuration |
| 1 | Load input dataset |
| 2 | Data cleaning Step 1: standardize competing/censoring columns |
| 3 | Data cleaning Step 2: construct observed data (obs_data) |
| 4 | Construct analytic dataset: variable selection, factor coding, reference levels |
| 5 | Clone: create arm-specific datasets |
| 6 | Fit IPACW models (months 1-3, per arm, numerator and denominator) |
| 7 | Compute cumulative IPACW; carry forward weights from month 3 |
| 8 | Fit IPCW models for informative right-censoring |
| 9 | Combine weights; truncate; freeze IPACW at month 3 |
| 10 | Apply artificial censoring; pool cloned arm datasets |
| 11 | Fit weighted outcome model; G-computation (point estimates) |
| 12 | Bootstrap resampling with full model re-estimation |
| 13 | Consolidate bootstrap results; compute 95% confidence intervals |
| 14 | Visualization: CIF curves; subgroup forest plots |

---

## Methods Summary

### Clone-Censor-Weight (CCW) Framework

The CCW framework emulates random treatment assignment from observational data through three steps:

**Step 1 - Clone**

Each eligible participant is cloned three times, once per treatment strategy:
- **S1 - Intensification:** adding at least one new antihypertensive drug class or increasing the dose of an existing antihypertensive agent relative to the pre-crash regimen
- **S2 - Deprescribing:** discontinuing at least one antihypertensive drug class or reducing the dose of an existing antihypertensive agent relative to the pre-crash regimen
- **S3 - No Change (reference):** maintaining the same antihypertensive drug classes and doses as the pre-crash regimen

All clones share identical baseline covariate values, outcome history, and competing event history at time zero.

**Step 2 - Censor**

During the 3-month strategy assessment window (months 1-3), a clone is artificially censored at the month when its observed medication dispensing deviates from its assigned strategy. The specific operationalization of adherence and deviation for each arm is described in Supplementary Tables S3 and S4 of the manuscript. Classification follows a clinical priority hierarchy: Intensification > Deprescribing > No Change.

**Step 3 - Weight**

Artificial censoring introduces selection bias that is corrected by two sets of stabilized inverse probability weights:

**(a) Inverse Probability of Artificial Censoring Weights (IPACW)**

Separate pooled logistic regression models are fitted for each arm (S1, S2, S3) and each month of the assessment window (months 1, 2, and 3) to model the probability that a clone remains adherent to its assigned strategy. This corresponds to modeling the probability that a person was not artificially censored given their observed covariate history.

- **Denominator model:** P(adherence = 1 at month t | treatment history [months 2-3 only], baseline and time-varying covariates)
- **Numerator model:** P(adherence = 1 at month t | treatment history [months 2-3 only]) - stabilizes the weights

Combined stabilized weights at the end of month 3 are carried forward (frozen) and applied unchanged throughout months 4-12, reflecting the assumption that treatment strategies were classified based on dispensing episodes during the 3-month grace period and that medication changes after this window were not accounted for in the analysis. We did not incorporate lags between exposure changes and outcomes, anticipating that the primary mechanisms through which antihypertensive medication changes affect MVC risk (e.g., blood pressure changes, dizziness, orthostatic hypotension) would be operative within 24 hours of medication changes.

**(b) Inverse Probability of Censoring Weights (IPCW) for Informative Right-Censoring**

A single pooled logistic regression model across all arms and months estimates the probability of remaining uncensored due to Medicare disenrollment. Unlike competing events (death and long-term care entry), these censoring events remove participants from the risk set and may be informative.

- **Denominator model:** P(not censored at t | month, baseline and time-varying covariates)
- **Numerator model:** P(not censored at t | month) - stabilizes the weights

Note: Death and long-term care (LTC) entry are not included in this model; they are handled as competing events under the total-effects estimand (see below).

The final combined stabilized weight at each person-month is: **SW_AC = IPACW x IPCW**, truncated at the 99th percentile per arm to limit influence of extreme weights.

### Competing Events

All-cause mortality and long-term care entry are treated as competing events under Hernan's total-effects estimand. When a competing event occurs, outcome flags are set to 0 (the outcome is unobserved but not missing) and IPACW weights are frozen and carried forward. This approach reflects the real-world interest in estimating the total causal effect of treatment strategy changes, including pathways operating through competing events.

### Outcome Model

A weighted pooled logistic marginal structural model (MSM) using `speedglm` with `quasibinomial(link = "logit")` family:

- **Main analysis:** `outcome ~ arm + as.factor(month) + arm:as.factor(month)`
- **EM analyses:** `outcome ~ arm * as.factor(month) * factor(em_variable)`

Weights applied: IPACW x IPCW combined (sw_ac_99_te).

### G-Computation

12-month cumulative incidence functions (CIFs) are estimated via g-computation: sequential individual-level survival chaining (product of (1 - monthly hazard)) averaged over the baseline population covariate distribution.

### Confidence Intervals

To construct 95% confidence intervals (CIs), we used non-parametric percentile-based bootstrapping with replacement at the individual level across 500 replicates, following the method described by Hanley and MacGibbon. All IPACW and IPCW models are re-fitted within each bootstrap replicate.

### Effect Modification Analyses

All effect modification analyses are pre-specified. Stratum-specific risk differences (RDs) and risk ratios (RRs) are estimated for:
- **Baseline ADRD status** (yes vs. no)
- **Number of antihypertensive drug classes at baseline** (1 class vs. 2 classes vs. 3 or more classes)

Results are displayed as forest plots of stratum-specific RRs and RDs at 12 months, with 95% bootstrap CIs.

---

## Key Variables

> **Important:** The analytic dataset is structured as a **person-month panel** (one row per unique person per month of follow-up). All variable values in this dataset reflect measurements or status for that specific monthly interval. Baseline (time-fixed) variables retain the same value across all months for a given person; time-varying variables may change from month to month.

Full variable definitions are provided in `data_documentation/RiskOfRecurrent_MVC.xlsx`. Key variables are summarized below.

### Strategy Adherence Indicators

| Variable | Description |
|---|---|
| `S1` | Per-person-month binary adherence indicator for the Intensification strategy (1 = clone is adherent in this month-interval; 0 = clone has deviated and is artificially censored) |
| `S2` | Per-person-month binary adherence indicator for the Deprescribing strategy |
| `S3` | Per-person-month binary adherence indicator for the No Change strategy (reference) |

### Primary Outcome Variable

| Variable | Description | Ascertainment |
|---|---|---|
| `final_crash_flag` | Per-person-month binary indicator for recurrent MVC (1 = crash occurred in this interval; 0 = no crash; NA = administratively censored) | NJ-SHO crash records probabilistically linked to Medicare |

### Censoring and Competing Event Variables

| Variable | Description |
|---|---|
| `censor_disenrol` | Per-person-month indicator for Medicare disenrollment censoring (1 = disenrolled in this interval) |
| `censor_study_end` | Per-person-month indicator for end-of-study censoring (1 = reached December 31, 2017 in this interval) |
| `competing_death` | Per-person-month indicator for all-cause mortality competing event (1 = died in this interval) |
| `competing_ltc` | Per-person-month indicator for long-term care entry competing event (1 = entered LTC in this interval) |

### Baseline (Time-Fixed) Covariates

| Variable | Description |
|---|---|
| `age_cat` | Age category at index crash: 66-69 (reference), 70-74, 75-79, 80-84, >=85 |
| `sex_base` | Sex (1=Male, 2=Female) |
| `race_base` | Race/ethnicity (categorical; reference: Non-Hispanic White) |
| `fault_base` | At-fault for index crash (binary) |
| `gagne_base` | Gagne comorbidity index (continuous) |
| `dual_base` | Dual Medicare-Medicaid eligibility (binary) |
| `meds_antihyp_cat` | Number of antihypertensive drug classes at baseline (1 [reference], 2, >=3) |
| `meds_pdi_cat` | Number of PDI drug classes at baseline (0 [reference], 1, 2, >=3) |
| `hosp_intensity_base` | Hospitalization frequency in 12-month lookback (None [reference], One, Multiple) |
| `ed_intensity_base` | ED visit frequency in 12-month lookback (None [reference], One, Multiple) |
| `op_intensity_base` | Outpatient visit frequency in 12-month lookback (None [reference], 1-9, 10-19, >=20) |
| `adrd_base` | Baseline Alzheimer disease and related dementias (binary) |
| `chf_base` | Baseline heart failure (binary) |
| `stroke_base` | Baseline stroke/TIA (binary) |
| `ami_base`, `anxi_base`, `chrkid_base`, `depr_base`, `diab_base`, `ischhd_base` | Additional baseline comorbidities (binary) |

### Time-Varying Covariates (Updated Monthly)

| Variable | Description |
|---|---|
| `hosp_updated`, `ed_updated`, `op_updated` | Updated healthcare utilization indicators |
| `ami_updated`, `anxi_updated`, `chrkid_updated`, `depr_updated`, `diab_updated`, `chf_updated`, `ischhd_updated`, `stroke_updated`, `adrd_updated` | Updated comorbidity indicators (binary) |

---

## Software and Package Requirements

**Primary language:** R (version 4.5.3 or higher recommended)
**Data management:** SAS (used for upstream data management and outcome ascertainment; SAS code available from corresponding author upon request)

### Required R Packages

```r
install.packages(c(
  "dplyr",       # Data manipulation
  "data.table",  # Fast in-memory data operations; required for weight computation
  "tidyr",       # Pivot operations for bootstrap output consolidation
  "ggplot2",     # Visualization
  "patchwork",   # Multi-panel figure composition
  "speedglm",    # Fast GLM fitting for large datasets
  "haven",       # Reading SAS datasets
  "stringr",     # String manipulation
  "survey",      # Survey-weighted analyses
  "tableone"     # Descriptive statistics and Table 1 generation
))
```

---

## Data Documentation

The file `data_documentation/RiskOfRecurrent_MVC.xlsx` contains the following sheets:

| Sheet | Contents |
|---|---|
| **Project Overview** | Study title, authors, funding, DUA number, data sources and years |
| **Cohort Development** | Step-by-step inclusion and exclusion criteria with participant counts at each step |
| **Variable Descriptions** | Variable name, source dataset, type, person-month level definition, use in IPACW/IPCW models |
| **Drug Class Lists** | Antihypertensive drug entries across 10 classes used to classify strategy changes; PDI drug class entries |

---

## Reproducibility Notes

1. **Bootstrap seed:** Set at the beginning of the bootstrap loop for reproducibility. See `CCW_template.R` Section 12.
2. **Path dependencies:** All absolute file paths must be updated to the user's local or server directory before execution.
3. **Data availability:** Input data are not publicly available due to CMS DUA restrictions. The code is provided to enable methodological replication with analogous datasets.
4. **R version:** Code was developed and tested under R 4.5.3. Minor syntax differences may exist in earlier versions.
5. **SAS upstream:** Data management, exposure construction (S1/S2/S3 strategy classification), and outcome ascertainment were performed in SAS on the Brown University Gerontology server. SAS code is available from the corresponding author upon reasonable request.
6. **Person-month structure:** All analytic steps assume a person-month panel as input. Each row represents one calendar month of follow-up for one participant. The temporal structure (ordering by person then month) is critical for correct computation of treatment history lags, cumulative weight products, and the competing event carry-forward logic.

---

## License

This repository is licensed under the **MIT License**. See `LICENSE` for full terms. You are free to use, modify, and distribute this code with attribution.

---

## Citation

> [To be completed after manuscript acceptance]
> Zhang S, Joyce N, D'Amico A, Oganisian A, Zullo AR, Harris DA, Hayes KN. Antihypertensive Treatment Strategy Changes Following a Motor Vehicle Crash and Risk of Recurrent Crash Among Older Medicare Drivers: A Target Trial Emulation. *Journal of the American Geriatrics Society*. [Year]; doi: [DOI]

**Repository citation:**
> Zhang S, et al. Code and data documentation for: Antihypertensive Treatment Strategy Changes Following a Motor Vehicle Crash and Risk of Recurrent Crash Among Older Medicare Drivers. GitHub repository. doi: [Zenodo DOI - to be issued]
