################################################################################
#                                                                              #
#   CLONE-CENSOR-WEIGHT (CCW) TARGET TRIAL EMULATION                          #
#   Generalized Analytic Template                                              #
#                                                                              #
#   Purpose:                                                                   #
#     Fully annotated, generalized template for implementing the CCW           #
#     framework for target trial emulation using observational person-time     #
#     panel data. Covers all steps from data cleaning through bootstrap        #
#     confidence interval estimation and visualization.                        #
#                                                                              #
#   Estimand:                                                                  #
#     Per-protocol total effect of sustained treatment strategies on           #
#     time-to-event outcomes. Competing events handled per Hernan's total-     #
#     effects framework: outcome is considered unobserved (not missing) when   #
#     a competing event occurs, and weights are carried forward.               #
#                                                                              #
#   Three-strategy application (adapt K as needed):                            #
#     S1 = Intensification (adding/increasing antihypertensive agents)         #
#     S2 = Deprescribing   (removing/reducing antihypertensive agents)         #
#     S3 = No Change       (reference strategy)                                #
#                                                                              #
#   Assumptions for use of this template:                                      #
#     - Input: person-month panel dataset (one row per person per month)       #
#     - A single primary outcome (generalize J as needed)                      #
#     - Competing events are present (e.g., death, long-term care entry)       #
#       and are handled per Hernan's total-effects/mediator theory             #
#     - Assessment window of T_ASSESS months; maximum follow-up T_MAX months   #
#     - All IPACW and IPCW models are re-fitted within each bootstrap replicate#
#                                                                              #
#   NOTE ON TEMPORAL ALIGNMENT:                                                #
#     Two valid approaches exist for person-month panel construction:          #
#     (A) Lagged covariates: covariates measured in the prior interval,        #
#         outcome and censoring in the current interval.                       #
#     (B) Current-month covariates: covariates and outcome in the same         #
#         interval.                                                             #
#     Both approaches should yield the same estimates; the choice depends on   #
#     the team's data structure and clinical assumptions. Document which        #
#     approach is used. This template assumes approach (B) by default.         #
#                                                                              #
#   HOW TO ADAPT THIS TEMPLATE:                                                #
#     1. Set all USER CONFIGURATION variables in Section 0 to your column     #
#        names (character strings).                                            #
#     2. Update BL_COVS and TV_COVS with your own covariate names.            #
#     3. List your covariate column names in BL_COV_NAMES and TV_COV_NAMES.   #
#     4. Set T_ASSESS, T_MAX, B_STEPS, SEED, PATH_DATA, PATH_OUT.             #
#                                                                              #
#   Reference for bootstrap CI method:                                         #
#     Hanley JA, MacGibbon B. Creating non-parametric bootstrap samples using  #
#     Poisson frequencies. Comput Methods Programs Biomed. 2006;83:57-62.      #
#                                                                              #
#   Authors:  [Sirui (Rita) Zhang]                                             #
#   Affiliation: [Brown University]                                            #
#   Date:     [2026-09-17]                                                     #
#   R Version:  4.6.0                                                          #
#                                                                              #
################################################################################


# =============================================================================
# SECTION 0: PACKAGE LOADING AND GLOBAL CONFIGURATION
# =============================================================================

required_pkgs <- c(
  "dplyr",       # Data manipulation
  "data.table",  # Fast in-memory operations; required throughout pipeline
  "tidyr",       # Pivot operations for bootstrap output
  "ggplot2",     # Visualization
  "patchwork",   # Multi-panel figure composition
  "speedglm",    # Fast GLM fitting for large cloned datasets
  "haven",       # Read SAS datasets (if applicable)
  "stringr",     # String utilities
  "tableone"     # Table 1 / descriptive statistics
)
for (pkg in required_pkgs) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg); library(pkg, character.only = TRUE)
  }
}

options(scipen = 999)  # Suppress scientific notation

# ---- File paths --------------------------------------------------------------
PATH_DATA <- "path/to/your/input/data/"
PATH_OUT  <- "path/to/your/output/folder/"

# ---- Analysis parameters -----------------------------------------------------
T_ASSESS <- 3      # Length of strategy assessment / grace period (months)
T_MAX    <- 12     # Maximum follow-up duration (months)
B_STEPS  <- 500    # Number of bootstrap replicates
SEED     <- 12345  # Random seed for reproducibility

# ---- Strategy labels ---------------------------------------------------------
STRAT_1 <- "Intensification"  # S1
STRAT_2 <- "Deprescribing"    # S2
STRAT_3 <- "No Change"        # S3 (reference)

# ---- USER CONFIGURATION: set column names as character strings ---------------
# Replace the values on the right-hand side with your actual column names.
# Do NOT use brackets [ ] -- use plain quoted strings.
# These are used via get() and .data[[]] throughout the pipeline.

PERSON_ID_VAR        <- "person_id"         # e.g., "bene_id_18900"
MONTH_VAR            <- "month"             # integer month (1 to T_MAX)
INDEX_DATE_VAR       <- "index_date"        # date of time zero (index event)
OUTCOME_FLAG         <- "outcome_flag"      # e.g., "final_crash_flag"
OUTCOME_LABEL        <- "outcome"           # e.g., "crash"
STRATEGY_1_COL       <- "S1"               # S1 adherence: 1=adherent, 0=deviated
STRATEGY_2_COL       <- "S2"               # S2 adherence: 1=adherent, 0=deviated
STRATEGY_3_COL       <- "S3"               # S3 adherence: 1=adherent, 0=deviated
CENSOR_DISENROL_VAR  <- "censor_disenrol"  # 1 = disenrolled this month
CENSOR_STUDY_END_VAR <- "censor_study_end" # 1 = reached study end date
COMPETING_DEATH_VAR  <- "competing_death"  # 1 = died this month
COMPETING_LTC_VAR    <- "competing_ltc"    # 1 = entered long-term care this month

# ---- Covariate strings -------------------------------------------------------
# BL_COVS: time-fixed baseline covariates (same value across all months per person)
# TV_COVS: time-varying covariates (may change month to month)
# Both are included in IPACW denominator models (months 1, 2, 3) AND
# in the IPCW denominator model (all months 1to T_MAX).
#
# Example from the antihypertensive MVC study:
#
# BL_COVS -- time-fixed baseline:
#   age_cat              age group: 66-69 (ref), 70-74, 75-79, 80-84, >=85
#   sex_base             sex (Male ref, Female)
#   race_base            race/ethnicity (Non-Hispanic White ref, ...)
#   fault_base           at-fault for index crash (binary)
#   gagne_base           Gagne comorbidity index (continuous)
#   dual_base            dual Medicare-Medicaid eligibility (binary)
#   meds_antihyp_cat     antihypertensive drug class count (1 ref, 2, >=3)
#   meds_pdi_cat         PDI drug class count (0 ref, 1, 2, >=3)
#   hosp_intensity_base  hospitalization frequency in 12-month lookback (None ref, One, Multiple)
#   ed_intensity_base    ED visit frequency in 12-month lookback (None ref, One, Multiple)
#   op_intensity_base    outpatient visit frequency (None ref, 1-9, 10-19, >=20)
#   ami_base, anxi_base, chrkid_base, depr_base, diab_base,
#   chf_base, ischhd_base, stroke_base, adrd_base  -- binary comorbidities
#
# TV_COVS -- time-varying (updated monthly):
#   hosp_updated, ed_updated, op_updated   healthcare utilization indicators
#   ami_updated, anxi_updated, chrkid_updated, depr_updated, diab_updated,
#   chf_updated, ischhd_updated, stroke_updated, adrd_updated ( those are comorbidity updates)
#
# Replace with your own variable names:

BL_COVS <- paste(
  "age_cat + sex_base + race_base + fault_base + gagne_base + dual_base +",
  "meds_antihyp_cat + meds_pdi_cat + hosp_intensity_base + ed_intensity_base + op_intensity_base +",
  "ami_base + anxi_base + chrkid_base + depr_base + diab_base + chf_base + ischhd_base +",
  "stroke_base + adrd_base"
)

TV_COVS <- paste(
  "hosp_updated + ed_updated + op_updated +",
  "ami_updated + anxi_updated + chrkid_updated + depr_updated + diab_updated +",
  "chf_updated + ischhd_updated + stroke_updated + adrd_updated"
)

# Column name vectors for dplyr::select -- list all covariate names as strings
BL_COV_NAMES <- c(
  "age_cat", "sex_base", "race_base", "fault_base", "gagne_base", "dual_base",
  "meds_antihyp_cat", "meds_pdi_cat", "hosp_intensity_base", "ed_intensity_base",
  "op_intensity_base", "ami_base", "anxi_base", "chrkid_base", "depr_base",
  "diab_base", "chf_base", "ischhd_base", "stroke_base", "adrd_base"
)
TV_COV_NAMES <- c(
  "hosp_updated", "ed_updated", "op_updated",
  "ami_updated", "anxi_updated", "chrkid_updated", "depr_updated", "diab_updated",
  "chf_updated", "ischhd_updated", "stroke_updated", "adrd_updated"
)


# =============================================================================
# SECTION 1: LOAD INPUT DATASET
# =============================================================================
# Input dataset structure:
#   - One row per PERSON_ID_VAR per follow-up MONTH_VAR
#   - Sorted by PERSON_ID_VAR, then MONTH_VAR
#   - Contains: outcome flags, censoring/competing event indicators,
#               time-fixed baseline covariates (BL_COV_NAMES),
#               time-varying covariates (TV_COV_NAMES),
#               strategy adherence indicators (S1, S2, S3)
#
# In the antihypertensive MVC study:
#   PERSON_ID_VAR = "bene_id_18900"   (Medicare beneficiary identifier)
#   MONTH_VAR     = "month"           (integer 1 to 12)

raw_data <- haven::read_sas(paste0(PATH_DATA, "YOUR_DATASET.sas7bdat"))
# Alternative: raw_data <- readRDS(paste0(PATH_DATA, "YOUR_DATASET.rds"))


# =============================================================================
# SECTION 2: DATA CLEANING STEP 1
# Standardize competing event and censoring column names for downstream use.
# =============================================================================

raw_data1 <- raw_data %>%
  arrange(.data[[PERSON_ID_VAR]], .data[[MONTH_VAR]]) %>%
  group_by(.data[[PERSON_ID_VAR]]) %>%
  mutate(
    # Map source columns to standard pipeline names used below
    # Replace right-hand side with your actual source column names if different
    competing_death = .data[[COMPETING_DEATH_VAR]],
    competing_ltc   = .data[[COMPETING_LTC_VAR]]
  ) %>%
  ungroup()


# =============================================================================
# SECTION 3: DATA CLEANING STEP 2
# Construct the observed data (obs_data) by applying the event priority hierarchy.
#
# Priority order (highest to lowest):
#
#   HIGHEST -- Competing events (death, LTC entry)
#     When competing_death = 1 or competing_ltc = 1 is first met:
#       -- Outcome flags are set to 0 (unobserved, NOT missing)
#       -- Per Hernan's total-effects framework, follow-up continues
#       -- Rows for subsequent months are NOT discarded
#       -- Weights are frozen and carried forward from the competing event month
#
#   MIDDLE -- Primary outcome events
#     When outcome_flag = 1:
#       -- Follow-up terminates after this month
#       -- Rows for subsequent months are discarded
#
#   LOWEST -- Study exit censoring (disenrollment, end of study period)
#     When censor_disenrol = 1 or censor_study_end = 1:
#       -- Outcome flags set to NA (observation not available)
#       -- Follow-up terminates after this month
#       -- Rows for subsequent months are discarded
#
# IMPORTANT: Step 1 (censoring -> NA) is applied BEFORE Step 2
# (competing events -> 0) so that competing events correctly overwrite NA
# when both occur in the same month.
# =============================================================================

obs_data <- raw_data1 %>%
  arrange(.data[[PERSON_ID_VAR]], .data[[MONTH_VAR]]) %>%
  group_by(.data[[PERSON_ID_VAR]]) %>%
  mutate(

    # Step 1: Study exit censoring -- set outcome to NA (applied first)
    !!sym(OUTCOME_FLAG) := ifelse(
      .data[[CENSOR_DISENROL_VAR]] == 1 | .data[[CENSOR_STUDY_END_VAR]] == 1,
      NA_real_, .data[[OUTCOME_FLAG]]
    ),

    # Step 2: Competing events -- set outcome to 0 (overwrites NA above)
    # Outcome = 0 signals event is unobserved (not missing) per Hernan's framework.
    # Example: when competing_death = 1 or competing_ltc = 1 is first met,
    # set outcome_flag to 0 for this month and all remaining months.
    has_had_competing   = cumsum(competing_death == 1 | competing_ltc == 1),
    !!sym(OUTCOME_FLAG) := ifelse(has_had_competing >= 1, 0, .data[[OUTCOME_FLAG]]),

    # Step 3: Identify months after which follow-up should stop
    # (outcome event OR study exit; NOT after competing events)
    stop_after_this_month = ifelse(
      .data[[OUTCOME_FLAG]] == 1 |
        is.na(.data[[OUTCOME_FLAG]]) |
        .data[[CENSOR_DISENROL_VAR]] == 1 |
        .data[[CENSOR_STUDY_END_VAR]] == 1,
      1, 0
    ),
    stop_after_this_month = coalesce(stop_after_this_month, 0),
    past_stops            = lag(cumsum(stop_after_this_month), default = 0)
  ) %>%

  # Step 4: Remove rows following an outcome or study exit event
  # Rows following a competing event are NOT removed
  filter(past_stops == 0) %>%
  dplyr::select(-stop_after_this_month, -past_stops, -has_had_competing) %>%
  ungroup()


# =============================================================================
# SECTION 4: CONSTRUCT ANALYTIC DATASET (anal_data)
# Select variables needed for CCW analysis.
# Convert categoricals to factors; set reference levels.
# =============================================================================

anal_data <- obs_data %>%
  dplyr::select(all_of(c(
    PERSON_ID_VAR,       # Unique participant identifier
    INDEX_DATE_VAR,      # Time-zero date (index event)
    MONTH_VAR,           # Month of follow-up (integer 1 to T_MAX)
    OUTCOME_FLAG,        # Primary outcome: 1=event, 0=no event, NA=censored
    "competing_death",   # Competing event: death
    "competing_ltc",     # Competing event: long-term care entry
    CENSOR_DISENROL_VAR, # Study exit: disenrollment
    CENSOR_STUDY_END_VAR,# Study exit: end of study period
    BL_COV_NAMES,        # Time-fixed baseline covariates (defined in Section 0)
    TV_COV_NAMES,        # Time-varying covariates (defined in Section 0)
    STRATEGY_1_COL,      # S1 adherence indicator (Intensification): 1=adherent, 0=deviated
    STRATEGY_2_COL,      # S2 adherence indicator (Deprescribing):    1=adherent, 0=deviated
    STRATEGY_3_COL       # S3 adherence indicator (No Change):        1=adherent, 0=deviated
  )))

# Convert to data.table for efficient operations throughout pipeline
anal_data <- as.data.table(anal_data)
setkeyv(anal_data, c(PERSON_ID_VAR, MONTH_VAR))

# ---- Factor coding examples from the antihypertensive MVC study --------------
# Replace variable names, levels, and labels for your covariates.
# All example lines below are commented out -- uncomment and adapt as needed.

# Age category (derived from continuous age_base, or pre-computed in SAS)
# anal_data[, age_cat := cut(
#   age_base,
#   breaks = c(65, 69, 74, 79, 84, Inf),
#   labels = c("66-69", "70-74", "75-79", "80-84", ">=85"),
#   right = TRUE, include.lowest = TRUE
# )]

# Sex
# anal_data[, sex_base := factor(sex_base,
#   levels = c(1, 2), labels = c("Male", "Female"))]

# Race/ethnicity
# anal_data[, race_base := factor(race_base,
#   levels = c(0, 1, 2, 3, 4, 5, 6),
#   labels = c("Unknown", "Non-Hispanic White", "Black/African-American",
#              "Other", "Asian/Pacific Islander", "Hispanic",
#              "American Indian/Alaska Native"))]

# Healthcare utilization intensity
# anal_data[, hosp_intensity_base := factor(hosp_intensity_base,
#   levels = c(0, 1, 2), labels = c("None", "One", "Multiple"))]
# anal_data[, ed_intensity_base := factor(ed_intensity_base,
#   levels = c(0, 1, 2), labels = c("None", "One", "Multiple"))]
# anal_data[, op_intensity_base := factor(op_intensity_base,
#   levels = c(0, 1, 2, 3), labels = c("None", "1-9", "10-19", ">=20"))]

# Antihypertensive drug class count (grouped)
# anal_data[, meds_antihyp_cat := factor(
#   cut(meds_antihyp_base, breaks = c(0, 1, 2, Inf),
#       labels = c("1", "2", ">=3"), right = TRUE),
#   levels = c("1", "2", ">=3")
# )]

# ---- Set reference levels ----------------------------------------------------
# anal_data$age_cat             <- relevel(anal_data$age_cat,             ref = "66-69")
# anal_data$sex_base            <- relevel(anal_data$sex_base,            ref = "Male")
# anal_data$race_base           <- relevel(anal_data$race_base,           ref = "Non-Hispanic White")
# anal_data$hosp_intensity_base <- relevel(anal_data$hosp_intensity_base, ref = "None")
# anal_data$ed_intensity_base   <- relevel(anal_data$ed_intensity_base,   ref = "None")
# anal_data$op_intensity_base   <- relevel(anal_data$op_intensity_base,   ref = "None")
# anal_data$meds_antihyp_cat    <- relevel(anal_data$meds_antihyp_cat,    ref = "1")

cat("Person-months in analytic dataset:", nrow(anal_data), "\n")
cat("Unique persons:", uniqueN(anal_data[[PERSON_ID_VAR]]), "\n")


# =============================================================================
# SECTION 5: CLONE -- CREATE ARM-SPECIFIC DATASETS
# Each eligible person is cloned once per strategy.
# All clones share identical covariate and outcome values at this stage.
# The 'arm' column labels which strategy each clone is assigned to.
# =============================================================================

arm_s1 <- copy(anal_data)[, arm := STRAT_1]  # Intensification clones
arm_s2 <- copy(anal_data)[, arm := STRAT_2]  # Deprescribing clones
arm_s3 <- copy(anal_data)[, arm := STRAT_3]  # No Change clones (reference)

# Verify clone counts at month 1 (all three should equal total unique N)
cat("Clone counts at month 1:\n")
cat("  S1 (Intensification):", sum(arm_s1[[MONTH_VAR]] == 1), "\n")
cat("  S2 (Deprescribing):",   sum(arm_s2[[MONTH_VAR]] == 1), "\n")
cat("  S3 (No Change):",       sum(arm_s3[[MONTH_VAR]] == 1), "\n")


# =============================================================================
# SECTION 6: INVERSE PROBABILITY OF ARTIFICIAL CENSORING WEIGHTS (IPACW)
# Also referred to as IPTW in some CCW literature. These weights correct for
# confounding by indication that arises when a clone deviates from its
# assigned strategy during the 3-month assessment window (artificial censoring).
#
# The adherence indicator (STRATEGY_k_COL = 1) means the person's observed
# medication behavior in that month-interval is consistent with strategy k.
# When STRATEGY_k_COL = 0, the clone is artificially censored.
#
# Rationale:
#   Artificial censoring is informative: sicker patients may be more or less
#   likely to intensify, deprescribe, or maintain their regimen. IPACW
#   reweights remaining adherent clones to represent the full eligible
#   population at each time point, correcting for this confounding.
#
# Model structure (fitted separately per arm and per assessment month):
#   Denominator: P(STRATEGY_k = 1 | tx_history [months 2-3], BL_COVS, TV_COVS)
#     Month 1: no prior treatment history; BL_COVS + TV_COVS only
#     Months 2-3: prior-month adherence lag + BL_COVS + TV_COVS
#   Numerator:   P(STRATEGY_k = 1 | tx_history [months 2-3])
#     Stabilizes the weights.
#
# SW_A(t) = product_{k=1}^{min(t, T_ASSESS)} [num(k) / den(k)]
# The month-3 cumulative SW_A is frozen and carried forward through months 4-12.
# (See Section 9 for the freeze implementation.)
# =============================================================================

# Helper: create 1-month prior treatment history lag
create_tx_history <- function(df, adherence_col, person_id_var, month_var) {
  adherence_sym <- sym(adherence_col)
  df %>%
    group_by(.data[[person_id_var]]) %>%
    arrange(.data[[month_var]], .by_group = TRUE) %>%
    mutate(tx_history = lag(factor(!!adherence_sym), n = 1)) %>%
    ungroup()
}

arm_s1 <- create_tx_history(arm_s1, STRATEGY_1_COL, PERSON_ID_VAR, MONTH_VAR)
arm_s2 <- create_tx_history(arm_s2, STRATEGY_2_COL, PERSON_ID_VAR, MONTH_VAR)
arm_s3 <- create_tx_history(arm_s3, STRATEGY_3_COL, PERSON_ID_VAR, MONTH_VAR)

# ---- IPACW models for S1 (Intensification) -----------------------------------
# STRATEGY_1_COL = 1: person is adherent to Intensification this month
# STRATEGY_1_COL = 0: person deviated; clone will be artificially censored
cat("Fitting IPACW models for S1 (Intensification)...\n")

# Month 1: no prior-month adherence history; BL_COVS + TV_COVS
den_s1_m1 <- glm(
  as.formula(paste0("(", STRATEGY_1_COL, " == 1) ~ ", BL_COVS, " + ", TV_COVS)),
  family = binomial(),
  data   = filter(arm_s1, .data[[MONTH_VAR]] == 1)
)
num_s1_m1 <- glm(
  as.formula(paste0("(", STRATEGY_1_COL, " == 1) ~ 1")),
  family = binomial(),
  data   = filter(arm_s1, .data[[MONTH_VAR]] == 1)
)

# Month 2: prior-month adherence history + BL_COVS + TV_COVS
den_s1_m2 <- glm(
  as.formula(paste0("(", STRATEGY_1_COL, " == 1) ~ tx_history + ", BL_COVS, " + ", TV_COVS)),
  family = binomial(),
  data   = filter(arm_s1, .data[[MONTH_VAR]] == 2)
)
num_s1_m2 <- glm(
  as.formula(paste0("(", STRATEGY_1_COL, " == 1) ~ tx_history")),
  family = binomial(),
  data   = filter(arm_s1, .data[[MONTH_VAR]] == 2)
)

# Month 3: same structure as month 2
den_s1_m3 <- glm(
  as.formula(paste0("(", STRATEGY_1_COL, " == 1) ~ tx_history + ", BL_COVS, " + ", TV_COVS)),
  family = binomial(),
  data   = filter(arm_s1, .data[[MONTH_VAR]] == 3)
)
num_s1_m3 <- glm(
  as.formula(paste0("(", STRATEGY_1_COL, " == 1) ~ tx_history")),
  family = binomial(),
  data   = filter(arm_s1, .data[[MONTH_VAR]] == 3)
)

# ---- IPACW models for S2 (Deprescribing) -------------------------------------
cat("Fitting IPACW models for S2 (Deprescribing)...\n")

den_s2_m1 <- glm(
  as.formula(paste0("(", STRATEGY_2_COL, " == 1) ~ ", BL_COVS, " + ", TV_COVS)),
  family = binomial(), data = filter(arm_s2, .data[[MONTH_VAR]] == 1)
)
num_s2_m1 <- glm(
  as.formula(paste0("(", STRATEGY_2_COL, " == 1) ~ 1")),
  family = binomial(), data = filter(arm_s2, .data[[MONTH_VAR]] == 1)
)
den_s2_m2 <- glm(
  as.formula(paste0("(", STRATEGY_2_COL, " == 1) ~ tx_history + ", BL_COVS, " + ", TV_COVS)),
  family = binomial(), data = filter(arm_s2, .data[[MONTH_VAR]] == 2)
)
num_s2_m2 <- glm(
  as.formula(paste0("(", STRATEGY_2_COL, " == 1) ~ tx_history")),
  family = binomial(), data = filter(arm_s2, .data[[MONTH_VAR]] == 2)
)
den_s2_m3 <- glm(
  as.formula(paste0("(", STRATEGY_2_COL, " == 1) ~ tx_history + ", BL_COVS, " + ", TV_COVS)),
  family = binomial(), data = filter(arm_s2, .data[[MONTH_VAR]] == 3)
)
num_s2_m3 <- glm(
  as.formula(paste0("(", STRATEGY_2_COL, " == 1) ~ tx_history")),
  family = binomial(), data = filter(arm_s2, .data[[MONTH_VAR]] == 3)
)

# ---- IPACW models for S3 (No Change; reference) ------------------------------
cat("Fitting IPACW models for S3 (No Change)...\n")

den_s3_m1 <- glm(
  as.formula(paste0("(", STRATEGY_3_COL, " == 1) ~ ", BL_COVS, " + ", TV_COVS)),
  family = binomial(), data = filter(arm_s3, .data[[MONTH_VAR]] == 1)
)
num_s3_m1 <- glm(
  as.formula(paste0("(", STRATEGY_3_COL, " == 1) ~ 1")),
  family = binomial(), data = filter(arm_s3, .data[[MONTH_VAR]] == 1)
)
den_s3_m2 <- glm(
  as.formula(paste0("(", STRATEGY_3_COL, " == 1) ~ tx_history + ", BL_COVS, " + ", TV_COVS)),
  family = binomial(), data = filter(arm_s3, .data[[MONTH_VAR]] == 2)
)
num_s3_m2 <- glm(
  as.formula(paste0("(", STRATEGY_3_COL, " == 1) ~ tx_history")),
  family = binomial(), data = filter(arm_s3, .data[[MONTH_VAR]] == 2)
)
den_s3_m3 <- glm(
  as.formula(paste0("(", STRATEGY_3_COL, " == 1) ~ tx_history + ", BL_COVS, " + ", TV_COVS)),
  family = binomial(), data = filter(arm_s3, .data[[MONTH_VAR]] == 3)
)
num_s3_m3 <- glm(
  as.formula(paste0("(", STRATEGY_3_COL, " == 1) ~ tx_history")),
  family = binomial(), data = filter(arm_s3, .data[[MONTH_VAR]] == 3)
)


# =============================================================================
# SECTION 7: COMPUTE CUMULATIVE IPACW (SW_A)
# SW_A(t) = product_{k=1}^{min(t, T_ASSESS)} [num(k) / den(k)]
# For months > T_ASSESS: ratio set to 1 (freeze handled in Section 9).
# =============================================================================

compute_sw_a <- function(arm_df,
                          num_m1, den_m1,
                          num_m2, den_m2,
                          num_m3, den_m3,
                          month_var    = MONTH_VAR,
                          person_id_var = PERSON_ID_VAR) {
  setDT(arm_df)

  arm_df[get(month_var) == 1, p_num_a := predict(num_m1, newdata = .SD, type = "response")]
  arm_df[get(month_var) == 1, p_den_a := predict(den_m1, newdata = .SD, type = "response")]
  arm_df[get(month_var) == 2, p_num_a := predict(num_m2, newdata = .SD, type = "response")]
  arm_df[get(month_var) == 2, p_den_a := predict(den_m2, newdata = .SD, type = "response")]
  arm_df[get(month_var) == 3, p_num_a := predict(num_m3, newdata = .SD, type = "response")]
  arm_df[get(month_var) == 3, p_den_a := predict(den_m3, newdata = .SD, type = "response")]
  # Months > T_ASSESS: ratio = 1; freeze carries forward in Section 9
  arm_df[get(month_var) > T_ASSESS, `:=`(p_num_a = 1, p_den_a = 1)]
  # Cumulative product within person
  arm_df[, sw_a := cumprod(p_num_a / p_den_a), by = person_id_var]
  return(arm_df)
}

arm_s1 <- compute_sw_a(arm_s1, num_s1_m1, den_s1_m1, num_s1_m2, den_s1_m2, num_s1_m3, den_s1_m3)
arm_s2 <- compute_sw_a(arm_s2, num_s2_m1, den_s2_m1, num_s2_m2, den_s2_m2, num_s2_m3, den_s2_m3)
arm_s3 <- compute_sw_a(arm_s3, num_s3_m1, den_s3_m1, num_s3_m2, den_s3_m2, num_s3_m3, den_s3_m3)


# =============================================================================
# SECTION 8: INVERSE PROBABILITY OF CENSORING WEIGHTS (IPCW)
#            for informative right-censoring due to disenrollment
#
# Rationale:
#   Participants who disenroll from Medicare 
#   period may differ systematically from those who remain, making their
#   exit potentially informative. IPCW reweights participants at each month
#   to represent the eligible population that has not yet exited.
#
#   NOTE: Death and long-term care entry are NOT in this model.
#   They are handled as competing events (weight carry-forward; Section 9).
#
# Model structure (pooled across all arms, all months 1-T_MAX):
#   Denominator: P(not censored at t | as.factor(month), BL_COVS, TV_COVS)
#   Numerator:   P(not censored at t | as.factor(month))
#
#   SW_C(t) = product_{k=1}^{t} [num_c(k) / den_c(k)]
#   Unlike IPACW, IPCW is updated at every month throughout follow-up (not frozen).
# =============================================================================

# Create study exit censoring indicator (1 = exited due to disenrollment)
anal_data[, censor_study_exit := fifelse(
  get(CENSOR_DISENROL_VAR) == 1, 1L, 0L
)]

# Pool all three arm datasets for IPCW model fitting
all_arms_pooled <- rbindlist(list(arm_s1, arm_s2, arm_s3), use.names = TRUE, fill = TRUE)
all_arms_pooled[, censor_study_exit := fifelse(
  get(CENSOR_DISENROL_VAR) == 1 | get(CENSOR_STUDY_END_VAR) == 1, 1L, 0L
)]

# Denominator: P(not censored | month, BL_COVS, TV_COVS)
mod_den_c <- glm(
  as.formula(paste0(
    "(censor_study_exit == 0) ~ as.factor(", MONTH_VAR, ") + ", BL_COVS, " + ", TV_COVS
  )),
  family = binomial(),
  data   = as.data.frame(all_arms_pooled)
)

# Numerator: P(not censored | month) -- marginal; stabilizes weights
mod_num_c <- glm(
  as.formula(paste0("(censor_study_exit == 0) ~ as.factor(", MONTH_VAR, ")")),
  family = binomial(),
  data   = as.data.frame(all_arms_pooled)
)

# Predict and attach IPCW probabilities; compute cumulative product
for (arm_name in c("arm_s1", "arm_s2", "arm_s3")) {
  arm_df <- get(arm_name)
  arm_df[, p_num_c := predict(mod_num_c, newdata = as.data.frame(arm_df), type = "response")]
  arm_df[, p_den_c := predict(mod_den_c, newdata = as.data.frame(arm_df), type = "response")]
  arm_df[, sw_c    := cumprod(p_num_c / p_den_c), by = PERSON_ID_VAR]
  assign(arm_name, arm_df)
}


# =============================================================================
# SECTION 9: COMBINE WEIGHTS, TRUNCATE, AND FREEZE
#
# Combined stabilized weight at each person-month:
#   SW_AC = IPACW (SW_A) x IPCW (SW_C)
#
# Truncation:
#   SW_AC is truncated at the 99th percentile per arm to limit influence
#   of extreme weights on the outcome model.
#
# IPACW freeze (carry-forward from month 3):
#   The SW_A component is frozen at its month-3 value and carried forward
#   unchanged through months 4-12. This reflects the assumption that
#   treatment strategy classification is finalized at the end of month 3
#   (the grace period boundary) and medication changes after month 3 are
#   not accounted for. Key assumption: mechanisms linking antihypertensive
#   medication changes to outcomes (e.g., blood pressure changes, dizziness)
#   operate within 24 hours of medication changes.
#
# IPCW is NOT frozen: it is updated every month throughout follow-up.
#
# Competing event carry-forward (total-effects estimand):
#   When a competing event (death or LTC entry) occurs, the COMBINED weight
#   SW_AC at that month is frozen and carried forward for all subsequent
#   months for that person, implementing Hernan's total-effects estimand.
# =============================================================================

# Helper: 99th percentile truncation
truncate_weights <- function(df, col, pct = 0.99) {
  t99 <- quantile(df[[col]], pct, na.rm = TRUE)
  df[[col]] <- pmin(df[[col]], t99)
  return(df)
}

# Helper: freeze combined weight at competing event months and carry forward
freeze_at_competing <- function(df_arm,
                                 person_id_var     = PERSON_ID_VAR,
                                 month_var         = MONTH_VAR,
                                 comp_death_var    = COMPETING_DEATH_VAR,
                                 comp_ltc_var      = COMPETING_LTC_VAR) {
  setDT(df_arm)
  setkeyv(df_arm, c(person_id_var, month_var))

  df_arm[, has_prior_competing := shift(
    cumsum(get(comp_death_var) == 1 | get(comp_ltc_var) == 1), n = 1, fill = 0
  ), by = person_id_var]

  df_arm[, freeze_month := min(
    ifelse(get(comp_death_var) == 1 | get(comp_ltc_var) == 1, get(month_var), Inf)
  ), by = person_id_var]

  df_arm[, weight_at_freeze := ifelse(get(month_var) == freeze_month, sw_ac_99, NA_real_)]
  df_arm[, weight_at_freeze := {
    mx <- suppressWarnings(max(weight_at_freeze, na.rm = TRUE))
    ifelse(is.infinite(mx), NA_real_, mx)
  }, by = person_id_var]

  df_arm[, sw_ac_99_te := fifelse(
    has_prior_competing >= 1 & !is.na(weight_at_freeze),
    weight_at_freeze,
    sw_ac_99
  )]
  df_arm[, c("has_prior_competing", "freeze_month", "weight_at_freeze") := NULL]
  return(df_arm)
}

# Apply to each arm
for (arm_name in c("arm_s1", "arm_s2", "arm_s3")) {
  arm_df <- get(arm_name)
  arm_df[, sw_ac := sw_a * sw_c]           # Combine IPACW and IPCW
  arm_df <- truncate_weights(arm_df, "sw_ac") # Truncate at 99th pct
  setnames(arm_df, "sw_ac", "sw_ac_99")
  arm_df <- freeze_at_competing(arm_df)     # Freeze for total-effects estimand
  assign(arm_name, arm_df)
}

# Weight diagnostics
cat("=== Weight diagnostics (sw_ac_99_te) ===\n")
for (arm_name in c("arm_s1", "arm_s2", "arm_s3")) {
  cat(arm_name, ":\n"); print(summary(get(arm_name)$sw_ac_99_te))
}


# =============================================================================
# SECTION 10: ARTIFICIAL CENSORING AND DATA POOLING
# Set outcome to NA on months where a clone deviates from its strategy.
# Discard all months after the deviation month.
# Stack three arm datasets into one pooled dataset for the outcome model.
# =============================================================================

apply_artificial_censoring <- function(arm_df, strategy_col,
                                        person_id_var = PERSON_ID_VAR,
                                        outcome_flag  = OUTCOME_FLAG) {
  setDT(arm_df)
  arm_df[, my_censor := fifelse(get(strategy_col) == 0, 1L, 0L)]
  arm_df[, post_dev  := cumsum(my_censor), by = person_id_var]
  # Set outcome to NA on the deviation row (artificial censoring)
  arm_df[my_censor == 1, (outcome_flag) := NA_real_]
  # Keep deviation row and all preceding rows; discard subsequent rows
  arm_df <- arm_df[post_dev <= 1]
  arm_df[, post_dev := NULL]
  return(arm_df)
}

arm_s1 <- apply_artificial_censoring(arm_s1, STRATEGY_1_COL)
arm_s2 <- apply_artificial_censoring(arm_s2, STRATEGY_2_COL)
arm_s3 <- apply_artificial_censoring(arm_s3, STRATEGY_3_COL)

# Stack all three cloned arm datasets
final_data <- rbindlist(list(arm_s1, arm_s2, arm_s3), use.names = TRUE, fill = TRUE)
final_data[, arm := relevel(factor(arm), ref = STRAT_3)]
arm_levels <- levels(final_data$arm)


# =============================================================================
# SECTION 11: OUTCOME MODEL AND G-COMPUTATION (POINT ESTIMATES)
#
# Filter criteria for the outcome model pool:
#   Include only rows that satisfy ALL of the following:
#     (1) NOT artificially censored (my_censor == 0)
#     (2) NOT censored due to disenrollment (CENSOR_DISENROL_VAR == 0)
#     (3) NOT censored due to end of study  (CENSOR_STUDY_END_VAR == 0)
#     (4) Outcome flag is NOT missing       (!is.na(outcome_flag))
#   Rows following competing events are RETAINED (outcome = 0, weights frozen).
#
# Outcome model:
#   Weighted pooled logistic MSM using speedglm with quasibinomial(logit).
#   quasibinomial suppresses non-integer weight warnings from sw_ac_99_te.
#   Weights: sw_ac_99_te (IPACW x IPCW, truncated, competing-event frozen).
#
# G-computation rationale:
#   The weighted MSM estimates conditional monthly hazards h(t | arm, month).
#   G-computation converts these to marginal (population-average) cumulative
#   incidence by:
#     (1) Predicting h(t) for each person under each counterfactual strategy
#     (2) Computing cumulative survival: S_i(t) = product_{k<=t}(1 - h_i(k))
#     (3) Averaging over the baseline population: CIF(t) = 1 - mean(S_i(t))
#   When the outcome model includes no person-level covariates (only arm x month),
#   all persons get identical hazard predictions and g-computation equals direct
#   prediction. However, g-computation provides a principled marginal estimator
#   that generalizes correctly if covariates are added to the outcome model,
#   and integrates naturally with bootstrap weights via weighted.mean().
# =============================================================================

pool <- final_data[
  my_censor == 0 &
    get(CENSOR_DISENROL_VAR)  == 0 &
    get(CENSOR_STUDY_END_VAR) == 0 &
    !is.na(get(OUTCOME_FLAG))
]

fit_outcome <- speedglm(
  as.formula(paste0(
    OUTCOME_FLAG, " ~ arm + as.factor(", MONTH_VAR, ") + arm:as.factor(", MONTH_VAR, ")"
  )),
  family  = quasibinomial(link = "logit"),
  data    = pool,
  weights = pool$sw_ac_99_te
)

# G-computation: sequential survival chaining
base_pop    <- final_data[get(MONTH_VAR) == 1]
results_cif <- list()

for (target_arm in arm_levels) {
  surv_prob <- rep(1, nrow(base_pop))
  cif_path  <- numeric(T_MAX)

  for (m in 1:T_MAX) {
    grid <- copy(base_pop)
    grid[, arm         := factor(target_arm, levels = arm_levels)]
    grid[, (MONTH_VAR) := m]
    haz <- tryCatch(
      predict(fit_outcome, newdata = as.data.frame(grid), type = "response"),
      error = function(e) rep(0, nrow(grid))
    )
    haz       <- pmin(pmax(haz, 0), 1)
    surv_prob <- surv_prob * (1 - haz)
    cif_path[m] <- 1 - mean(surv_prob, na.rm = TRUE)
  }

  results_cif[[target_arm]] <- data.frame(
    arm   = target_arm,
    month = 1:T_MAX,
    risk  = cif_path
  )
}

cif_pe <- bind_rows(results_cif)
write.csv(cif_pe, file.path(PATH_OUT, "cif_point_estimates.csv"), row.names = FALSE)
cat("Point estimates saved.\n")


# =============================================================================
# SECTION 12: BOOTSTRAP RESAMPLING WITH FULL MODEL RE-ESTIMATION
#
# To construct 95% confidence intervals (CIs), we use non-parametric
# percentile-based bootstrapping with replacement at the individual level
# across 500 replicates, following the method described by Hanley and MacGibbon.
#
# Within each replicate:
#   (a) A Poisson(lambda=1) frequency weight is drawn per unique person.
#       Persons with weight 0 are excluded (mathematically equivalent to
#       sampling with replacement).
#   (b) All IPACW and IPCW models are re-fitted on the resampled data.
#   (c) Artificial censoring is re-applied; combined weights are recomputed.
#   (d) Outcome model is re-fitted; g-computation is repeated.
#   (e) CIF, risk differences (RD), and risk ratios (RR) are stored.
#
# Reference: Hanley JA, MacGibbon B. Creating non-parametric bootstrap samples
# using Poisson frequencies. Comput Methods Programs Biomed. 2006;83:57-62.
# =============================================================================

# Save and reload analytic dataset (used as the clean base for each replicate)
saveRDS(anal_data, file.path(PATH_OUT, "anal_data.rds"))
anal_data_master <- readRDS(file.path(PATH_OUT, "anal_data.rds"))
setDT(anal_data_master)

unique_ids <- unique(anal_data_master[[PERSON_ID_VAR]])
N_persons  <- length(unique_ids)

# Safe GLM: falls back to intercept-only on convergence failure
safe_glm <- function(formula_str, data_in) {
  tryCatch(
    glm(as.formula(formula_str), family = binomial(), data = as.data.frame(data_in)),
    warning = function(w) suppressWarnings(
      glm(as.formula(formula_str), family = binomial(), data = as.data.frame(data_in))
    ),
    error = function(e) {
      message("GLM failed (", conditionMessage(e), ") -- intercept fallback")
      glm(as.formula(sub("~.*", "~ 1", formula_str)), family = binomial(),
          data = as.data.frame(data_in))
    }
  )
}

boot_results <- list()
set.seed(SEED)

for (b in seq_len(B_STEPS)) {
  cat(sprintf("Bootstrap replicate: %d / %d\n", b, B_STEPS))

  # Draw Poisson(1) weights per person; drop zero-weight persons
  poi_ids         <- data.table(tmp = unique_ids, boot_weight = rpois(N_persons, lambda = 1))
  setnames(poi_ids, "tmp", PERSON_ID_VAR)
  poi_ids         <- poi_ids[boot_weight > 0]
  b_data          <- merge(anal_data_master, poi_ids, by = PERSON_ID_VAR)
  setkeyv(b_data, c(PERSON_ID_VAR, MONTH_VAR))

  # Treatment history lags
  b_data[, tx_hist_s1 := shift(factor(get(STRATEGY_1_COL)), n = 1), by = PERSON_ID_VAR]
  b_data[, tx_hist_s2 := shift(factor(get(STRATEGY_2_COL)), n = 1), by = PERSON_ID_VAR]
  b_data[, tx_hist_s3 := shift(factor(get(STRATEGY_3_COL)), n = 1), by = PERSON_ID_VAR]
  b_data[, censor_study_exit := fifelse(
    get(CENSOR_DISENROL_VAR) == 1 | get(CENSOR_STUDY_END_VAR) == 1, 1L, 0L
  )]

  # Re-fit IPACW models -- S1 (Intensification)
  den_s1_m1_b <- safe_glm(paste0("(", STRATEGY_1_COL, "==1) ~ ", BL_COVS, " + ", TV_COVS),
                            b_data[get(MONTH_VAR) == 1])
  num_s1_m1_b <- safe_glm(paste0("(", STRATEGY_1_COL, "==1) ~ 1"),
                            b_data[get(MONTH_VAR) == 1])
  den_s1_m2_b <- safe_glm(paste0("(", STRATEGY_1_COL, "==1) ~ tx_hist_s1 + ", BL_COVS, " + ", TV_COVS),
                            b_data[get(MONTH_VAR) == 2])
  num_s1_m2_b <- safe_glm(paste0("(", STRATEGY_1_COL, "==1) ~ tx_hist_s1"),
                            b_data[get(MONTH_VAR) == 2])
  den_s1_m3_b <- safe_glm(paste0("(", STRATEGY_1_COL, "==1) ~ tx_hist_s1 + ", BL_COVS, " + ", TV_COVS),
                            b_data[get(MONTH_VAR) == 3])
  num_s1_m3_b <- safe_glm(paste0("(", STRATEGY_1_COL, "==1) ~ tx_hist_s1"),
                            b_data[get(MONTH_VAR) == 3])

  # Re-fit IPACW models -- S2 (Deprescribing)
  den_s2_m1_b <- safe_glm(paste0("(", STRATEGY_2_COL, "==1) ~ ", BL_COVS, " + ", TV_COVS),
                            b_data[get(MONTH_VAR) == 1])
  num_s2_m1_b <- safe_glm(paste0("(", STRATEGY_2_COL, "==1) ~ 1"),
                            b_data[get(MONTH_VAR) == 1])
  den_s2_m2_b <- safe_glm(paste0("(", STRATEGY_2_COL, "==1) ~ tx_hist_s2 + ", BL_COVS, " + ", TV_COVS),
                            b_data[get(MONTH_VAR) == 2])
  num_s2_m2_b <- safe_glm(paste0("(", STRATEGY_2_COL, "==1) ~ tx_hist_s2"),
                            b_data[get(MONTH_VAR) == 2])
  den_s2_m3_b <- safe_glm(paste0("(", STRATEGY_2_COL, "==1) ~ tx_hist_s2 + ", BL_COVS, " + ", TV_COVS),
                            b_data[get(MONTH_VAR) == 3])
  num_s2_m3_b <- safe_glm(paste0("(", STRATEGY_2_COL, "==1) ~ tx_hist_s2"),
                            b_data[get(MONTH_VAR) == 3])

  # Re-fit IPACW models -- S3 (No Change)
  den_s3_m1_b <- safe_glm(paste0("(", STRATEGY_3_COL, "==1) ~ ", BL_COVS, " + ", TV_COVS),
                            b_data[get(MONTH_VAR) == 1])
  num_s3_m1_b <- safe_glm(paste0("(", STRATEGY_3_COL, "==1) ~ 1"),
                            b_data[get(MONTH_VAR) == 1])
  den_s3_m2_b <- safe_glm(paste0("(", STRATEGY_3_COL, "==1) ~ tx_hist_s3 + ", BL_COVS, " + ", TV_COVS),
                            b_data[get(MONTH_VAR) == 2])
  num_s3_m2_b <- safe_glm(paste0("(", STRATEGY_3_COL, "==1) ~ tx_hist_s3"),
                            b_data[get(MONTH_VAR) == 2])
  den_s3_m3_b <- safe_glm(paste0("(", STRATEGY_3_COL, "==1) ~ tx_hist_s3 + ", BL_COVS, " + ", TV_COVS),
                            b_data[get(MONTH_VAR) == 3])
  num_s3_m3_b <- safe_glm(paste0("(", STRATEGY_3_COL, "==1) ~ tx_hist_s3"),
                            b_data[get(MONTH_VAR) == 3])

  # Re-fit IPCW models (pooled, all arms, all months 1-T_MAX)
  mod_den_c_b <- safe_glm(
    paste0("(censor_study_exit==0) ~ as.factor(", MONTH_VAR, ") + ", BL_COVS, " + ", TV_COVS),
    b_data
  )
  mod_num_c_b <- safe_glm(
    paste0("(censor_study_exit==0) ~ as.factor(", MONTH_VAR, ")"),
    b_data
  )

  # Create arm datasets; compute weights; apply artificial censoring
  b_s1 <- copy(b_data)[, arm := STRAT_1]
  b_s2 <- copy(b_data)[, arm := STRAT_2]
  b_s3 <- copy(b_data)[, arm := STRAT_3]

  b_s1 <- create_tx_history(b_s1, STRATEGY_1_COL, PERSON_ID_VAR, MONTH_VAR)
  b_s2 <- create_tx_history(b_s2, STRATEGY_2_COL, PERSON_ID_VAR, MONTH_VAR)
  b_s3 <- create_tx_history(b_s3, STRATEGY_3_COL, PERSON_ID_VAR, MONTH_VAR)

  b_s1 <- compute_sw_a(b_s1, num_s1_m1_b, den_s1_m1_b, num_s1_m2_b, den_s1_m2_b, num_s1_m3_b, den_s1_m3_b)
  b_s2 <- compute_sw_a(b_s2, num_s2_m1_b, den_s2_m1_b, num_s2_m2_b, den_s2_m2_b, num_s2_m3_b, den_s2_m3_b)
  b_s3 <- compute_sw_a(b_s3, num_s3_m1_b, den_s3_m1_b, num_s3_m2_b, den_s3_m2_b, num_s3_m3_b, den_s3_m3_b)

  for (arm_name in c("b_s1", "b_s2", "b_s3")) {
    arm_df <- get(arm_name)
    arm_df[, p_num_c := predict(mod_num_c_b, newdata = as.data.frame(arm_df), type = "response")]
    arm_df[, p_den_c := predict(mod_den_c_b, newdata = as.data.frame(arm_df), type = "response")]
    arm_df[, sw_c    := cumprod(p_num_c / p_den_c), by = PERSON_ID_VAR]
    arm_df[, sw_ac   := sw_a * sw_c]
    arm_df <- truncate_weights(arm_df, "sw_ac")
    setnames(arm_df, "sw_ac", "sw_ac_99")
    arm_df <- freeze_at_competing(arm_df)
    assign(arm_name, arm_df)
  }

  b_s1 <- apply_artificial_censoring(b_s1, STRATEGY_1_COL)
  b_s2 <- apply_artificial_censoring(b_s2, STRATEGY_2_COL)
  b_s3 <- apply_artificial_censoring(b_s3, STRATEGY_3_COL)

  pool_b <- rbindlist(list(b_s1, b_s2, b_s3), use.names = TRUE, fill = TRUE)
  pool_b[, arm := relevel(factor(arm), ref = STRAT_3)]
  pool_b <- pool_b[
    my_censor == 0 &
      get(CENSOR_DISENROL_VAR)  == 0 &
      get(CENSOR_STUDY_END_VAR) == 0 &
      !is.na(get(OUTCOME_FLAG))
  ]

  fit_b <- tryCatch(
    speedglm(
      as.formula(paste0(
        OUTCOME_FLAG, " ~ arm + as.factor(", MONTH_VAR, ") + arm:as.factor(", MONTH_VAR, ")"
      )),
      family  = quasibinomial(link = "logit"),
      data    = pool_b,
      weights = pool_b$sw_ac_99_te
    ),
    error = function(e) glm(
      as.formula(paste0(
        OUTCOME_FLAG, " ~ arm + as.factor(", MONTH_VAR, ") + arm:as.factor(", MONTH_VAR, ")"
      )),
      family  = quasibinomial(link = "logit"),
      data    = as.data.frame(pool_b),
      weights = pool_b$sw_ac_99_te
    )
  )

  base_pop_b <- b_data[get(MONTH_VAR) == 1]
  iter_res   <- list()

  for (target_arm in arm_levels) {
    surv_prob <- rep(1, nrow(base_pop_b))
    cif_path  <- numeric(T_MAX)
    for (m in 1:T_MAX) {
      grid <- copy(base_pop_b)
      grid[, arm         := factor(target_arm, levels = arm_levels)]
      grid[, (MONTH_VAR) := m]
      haz <- tryCatch(
        predict(fit_b, newdata = as.data.frame(grid), type = "response"),
        error = function(e) rep(0, nrow(grid))
      )
      haz       <- pmin(pmax(haz, 0), 1)
      surv_prob <- surv_prob * (1 - haz)
      cif_path[m] <- 1 - weighted.mean(surv_prob, w = base_pop_b$boot_weight, na.rm = TRUE)
    }
    iter_res[[target_arm]] <- data.frame(
      arm = target_arm, month = 1:T_MAX, risk = cif_path, boot_idx = b
    )
  }
  boot_results[[b]] <- bind_rows(iter_res)
}  # end bootstrap loop


# =============================================================================
# SECTION 13: CONSOLIDATE BOOTSTRAP RESULTS AND COMPUTE 95% CIs
# =============================================================================

df_boot <- bind_rows(boot_results)

# CIF summary: percentile CI per arm x month
cif_summary <- df_boot %>%
  group_by(arm, month) %>%
  summarise(
    cif_est   = mean(risk, na.rm = TRUE),
    cif_lower = quantile(risk, 0.025, na.rm = TRUE),
    cif_upper = quantile(risk, 0.975, na.rm = TRUE),
    .groups   = "drop"
  )

# Contrast summary: RD and RR per month with 95% CI
col_s1 <- paste0("r_", STRAT_1)
col_s2 <- paste0("r_", STRAT_2)
col_s3 <- paste0("r_", STRAT_3)

contrast_summary <- df_boot %>%
  pivot_wider(id_cols = c(month, boot_idx), names_from = arm,
              values_from = risk, names_prefix = "r_") %>%
  mutate(
    RD_s1 = .data[[col_s1]] - .data[[col_s3]],
    RR_s1 = .data[[col_s1]] / .data[[col_s3]],
    RD_s2 = .data[[col_s2]] - .data[[col_s3]],
    RR_s2 = .data[[col_s2]] / .data[[col_s3]]
  ) %>%
  pivot_longer(cols = matches("^(r_|RD_|RR_)"),
               names_to = "metric", values_to = "value") %>%
  group_by(month, metric) %>%
  summarise(
    pe    = mean(value, na.rm = TRUE),
    lower = quantile(value, 0.025, na.rm = TRUE),
    upper = quantile(value, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(cif_summary,      file.path(PATH_OUT, "bootstrap_cif_summary.csv"),      row.names = FALSE)
write.csv(contrast_summary, file.path(PATH_OUT, "bootstrap_contrast_summary.csv"), row.names = FALSE)
cat("Bootstrap outputs saved.\n")


# =============================================================================
# SECTION 14: VISUALIZATION
# (a) CIF curves with 95% bootstrap confidence intervals
# (b) Subgroup forest plots for EM analyses (stratum-specific RRs and RDs)
# =============================================================================

library(ggplot2)
library(dplyr)

# ---- Aesthetic mapping -------------------------------------------------------
STRAT_INTENS <- "Intensification"
STRAT_DE     <- "Deprescribing"
STRAT_REF    <- "No Change (ref)"

arm_cols <- c(
  "Intensification"  = "#D55E00",   # Okabe-Ito vermillion (colorblind-safe)
  "Deprescribing"    = "#0072B2",   # Okabe-Ito blue
  "No Change (ref)"  = "#000000"    # Black (clearly distinct from both)
)
arm_ltys <- c(
  "Intensification"  = "solid",
  "Deprescribing"    = "solid",
  "No Change (ref)"  = "12"         # Hex linetype: 1-unit on, 2-unit off
)
arm_code_to_label <- c(
  "Intensification"  = STRAT_INTENS,
  "Deprescribing"    = STRAT_DE,
  "No Change"        = STRAT_REF
)

# ---- (a) CIF curve: recurrent MVC -------------------------------------------
main_plot_data <- cif_summary %>%
  mutate(
    arm_clean = gsub("[_ ]", "-", arm),
    strategy  = factor(
      arm_code_to_label[arm],
      levels = c(STRAT_REF, STRAT_DE, STRAT_INTENS)
    )
  )

p_cif <- ggplot(
  main_plot_data,
  aes(x = month, y = cif_est * 100,
      colour = strategy, fill = strategy, linetype = strategy)
) +
  geom_ribbon(aes(ymin = cif_lower * 100, ymax = cif_upper * 100),
              alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.95) +
  scale_colour_manual(values = arm_cols, name = NULL) +
  scale_fill_manual(  values = arm_cols, name = NULL) +
  scale_linetype_manual(values = arm_ltys, name = NULL) +
  scale_x_continuous(
    breaks = c(2, 4, 6, 8, 10, 12),
    labels = c("M2", "M4", "M6", "M8", "M10", "M12")
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    limits = c(0, NA)
  ) +
  labs(
    x       = "Month of Follow-up",
    y       = "Cumulative Incidence of Second MVC (%)",
    title   = NULL,
    caption = "MVC = motor vehicle crash."
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position   = "bottom",
    legend.key.width  = unit(2.5, "cm"),
    legend.key.height = unit(0.5, "cm"),
    legend.text       = element_text(size = 10),
    panel.grid.minor  = element_blank(),
    plot.caption      = element_text(size = 8, colour = "grey40",
                                     hjust = 0, margin = margin(t = 6))
  ) +
  guides(
    colour   = guide_legend(
      override.aes = list(
        linewidth = 1.5,
        linetype  = c("12", "solid", "solid")  # matches factor level order
      )
    ),
    linetype = "none",
    fill     = "none"
  )

ggsave(file.path(PATH_OUT, "CIF_main_crash.png"), p_cif, width = 7, height = 5, dpi = 300)
ggsave(file.path(PATH_OUT, "CIF_main_crash.pdf"), p_cif, width = 7, height = 5, dpi = 300)
cat("CIF plot saved.\n")


# ---- (b) Subgroup forest plot: EM analyses -----------------------------------
# EM variables in this study:
#   (1) Baseline ADRD status (adrd_base: 0 = No, 1 = Yes)
#   (2) Number of antihypertensive drug classes at baseline
#       (meds_antihyp_cat: 1, 2, >=3)
#
# em_summary should contain stratum-specific RRs and RDs with 95% bootstrap CIs
# from separate EM bootstrap runs. Expected structure:
#   modifier | strata | comparison | metric | pe | lower | upper
#   Where metric is one of: RR_s1, RR_s2, RD_s1, RD_s2

arm_cols_forest <- c(
  "Intensification vs No Change" = "#D55E00",
  "Deprescribing vs No Change"   = "#0072B2"
)

# Uncomment and adapt once em_summary is available from EM bootstrap run:
# p_forest <- ggplot(
#   em_summary %>% filter(metric %in% c("RR_s1", "RR_s2")),
#   aes(x = pe, xmin = lower, xmax = upper, y = strata, colour = comparison)
# ) +
#   geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40", linewidth = 0.8) +
#   geom_errorbarh(height = 0.25, linewidth = 0.9) +
#   geom_point(size = 3.5, shape = 18) +
#   scale_x_continuous(trans = "log10") +
#   scale_colour_manual(values = arm_cols_forest, name = "Comparison") +
#   facet_grid(modifier ~ ., scales = "free_y", space = "free_y") +
#   labs(
#     x       = "Risk Ratio vs No Change (log scale)",
#     y       = NULL,
#     caption = "Diamond = point estimate. Bars = 95% bootstrap CI. Dashed line = null (RR = 1.0)."
#   ) +
#   theme_bw(base_size = 11) +
#   theme(legend.position = "bottom", panel.grid.minor = element_blank())
#
# ggsave(file.path(PATH_OUT, "forest_subgroup_RR.png"), p_forest, width = 9, height = 7, dpi = 300)

cat("\n=== CCW ANALYSIS COMPLETE ===\n")
cat("All outputs written to:", PATH_OUT, "\n")
