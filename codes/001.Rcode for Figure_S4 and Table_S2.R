# ==============================================================================
# Manuscript:   Cultivation facilitates global naturalization
# Description:  Bayesian Structural Equation Models (BSEM) examining direct and
#               indirect effects of native range size on naturalization success
#               across life-form groups (annual herbs, perennial herbs, woody).
# Author:       YW and DBC
# Date:         2026-01-16
# ==============================================================================

# 1. Environment Setup ---------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  tidyr,
  here,          # Reproducible file paths
  readr,         # Fast CSV reading
  dplyr,         # Data manipulation
  stringr,       # String operations (if needed later)
  purrr,         # Functional programming
  tidyr,         # Data reshaping
  ggplot2,       # Visualization
  ggridges,      # Ridge density plots
  ggpubr,        # Publication-ready plots
  patchwork,     # Plot composition
  scales,        # Axis formatting
  brms,          # Bayesian SEM
  broom.mixed,   # Tidy model outputs
  posterior      # Posterior processing
)

# Rationale: Set project root for all relative paths
setwd(here())
getwd()

# 2. Data Loading & Preprocessing ----------------------------------------------
# Rationale: Load matched native-cultivated plant dataset with standardized
#            TDWG Level 3 geographic units and trait data

native_raw <- read_csv(
  here("data", "20250530.native_plant.matched.cultivated_plant_DO11.csv"),
  show_col_types = FALSE
)

# Check raw variable ranges (quality control)
native_raw %>%
  select(planting.China.tdwg3, WorldCuP.n.tdwg3, 
         native.global.tdwg3, nat.extent) %>%
  purrr::map(range)

# --- Transformation Function: Log + Standardize ---
# Rationale: Log-transformation reduces right-skew in count/area data;
#            z-scaling enables comparison of effect sizes across predictors
scale_predictors <- function(df) {
  df %>%
    mutate(
      cul.China.tdwg3.scaled     = scale(log(planting.China.tdwg3 + 1)),
      cul.foreign.tdwg3.scaled   = scale(log(WorldCuP.n.tdwg3 + 1)),
      native.global.tdwg3.scaled = scale(log(native.global.tdwg3 + 1)),
      nat.extent.scaled          = scale(log(nat.extent + 1))
    )
}

# Apply transformations to full dataset
native_all <- scale_predictors(native_raw)

# Subset by life-form groups for hypothesis testing
# Rationale: Life-history strategies (annual/perennial/woody) may exhibit
#            distinct cultivation-naturalization pathway dependencies
native_annual    <- native_raw %>% 
  filter(life.form.integrated == "annual herb") %>% 
  scale_predictors()

native_perennial <- native_raw %>% 
  filter(life.form.integrated == "perennial herb") %>% 
  scale_predictors()

native_woody     <- native_raw %>% 
  filter(life.form.integrated == "woody") %>% 
  scale_predictors()


# 3. Bayesian Structural Equation Modeling (BSEM) -----------------------------
# Hypothesis: Native range size (NRS) affects naturalization success (NatSuc)
#             via three pathways:
#             1. Direct effect: NRS → NatSuc
#             2. Indirect via cultivation within China (CWC): NRS → CWC → NatSuc
#             3. Indirect via cultivation outside China (COC): NRS → CWC → COC → NatSuc

# --- Path Model Specification ---
# Rationale: Gaussian family appropriate for z-scaled continuous responses;
#            set_rescor(FALSE) assumes independence among response residuals
formula_cul_china    <- brms::bf(
  cul.China.tdwg3.scaled ~ native.global.tdwg3.scaled, 
  family = gaussian()
)

formula_cul_foreign  <- brms::bf(
  cul.foreign.tdwg3.scaled ~ cul.China.tdwg3.scaled + native.global.tdwg3.scaled,
  family = gaussian()
)

formula_nat_success  <- brms::bf(
  nat.extent.scaled ~ cul.foreign.tdwg3.scaled + cul.China.tdwg3.scaled + 
    native.global.tdwg3.scaled,
  family = gaussian()
)

# --- Computational Setup ---
# Rationale: Parallel chains improve sampling efficiency; seed ensures reproducibility
total_cores         <- parallel::detectCores(logical = FALSE)
chains_n            <- 4
cores_per_chain     <- max(1, total_cores / chains_n)

# --- Model Fitting (4 Groups) ---
# Rationale: Separate models allow life-form-specific path coefficient estimation
#            File argument auto-saves fitted objects to avoid re-running (2000 iter × 4 chains)

# Model: All species combined
bsem_natsucc_all <- brms::brm(
  formula_cul_china + formula_cul_foreign + formula_nat_success + set_rescor(FALSE),
  data   = native_all,
  chains = chains_n,
  cores  = cores_per_chain,
  iter   = 2000,
  seed   = 2024,
  file   = here("models", "bsem_natsucc_all.rds")
)

# Model: Annual herbs
bsem_natsucc_annual <- brms::brm(
  formula_cul_china + formula_cul_foreign + formula_nat_success + set_rescor(FALSE),
  data   = native_annual,
  chains = chains_n,
  cores  = cores_per_chain,
  iter   = 2000,
  seed   = 2024,
  file   = here("models", "bsem_natsucc_annual.rds")
)

# Model: Perennial herbs
bsem_natsucc_perennial <- brms::brm(
  formula_cul_china + formula_cul_foreign + formula_nat_success + set_rescor(FALSE),
  data   = native_perennial,
  chains = chains_n,
  cores  = cores_per_chain,
  iter   = 2000,
  seed   = 2024,
  file   = here("models", "bsem_natsucc_perennial.rds")
)

# Model: Woody plants
bsem_natsucc_woody <- brms::brm(
  formula_cul_china + formula_cul_foreign + formula_nat_success + set_rescor(FALSE),
  data   = native_woody,
  chains = chains_n,
  cores  = cores_per_chain,
  iter   = 2000,
  seed   = 2024,
  file   = here("models", "bsem_natsucc_woody.rds")
)

# --- Model Summaries ---
summary(bsem_natsucc_all)
summary(bsem_natsucc_annual)
summary(bsem_natsucc_perennial)
summary(bsem_natsucc_woody)

# --- Reproducibility: Save Consolidated Model Archive ---
# Rationale: Single .Rdata file enables batch loading for future analyses
save(
  bsem_natsucc_all,
  bsem_natsucc_annual,
  bsem_natsucc_perennial,
  bsem_natsucc_woody,
  file = here("models", "bsem_natsucc_all_groups_v202601.Rdata")
)

# --- Model Comparison: Bayesian R2 ---
# Rationale: Assess variance explained in each response variable across groups
set.seed(2024)
r2_annual    <- bayes_R2(bsem_natsucc_annual)
r2_perennial <- bayes_R2(bsem_natsucc_perennial)
r2_woody     <- bayes_R2(bsem_natsucc_woody)

# Consolidate R2 estimates
r2_comparison <- bind_rows(
  r2_annual    %>% as.data.frame() %>% mutate(Variable = rownames(.), Group = "Annual"),
  r2_perennial %>% as.data.frame() %>% mutate(Variable = rownames(.), Group = "Perennial"),
  r2_woody     %>% as.data.frame() %>% mutate(Variable = rownames(.), Group = "Woody"),
  .id = "Response"
) %>%
  mutate(
    Variable = str_remove(Variable, "\\.\\.\\.\\d+$"),
    Variable = case_when(
      Variable == "R2culChinatdwg3scaled" ~ "CWC",
      Variable == "R2culforeigntdwg3scaled" ~ "COC",
      Variable == "R2natextentscaled" ~ "NatSuc",
      TRUE ~ Variable
    )
  )


print(r2_comparison)

write_csv(r2_comparison, here("results", "r2_comparison.csv"))

# 4. Posterior Path Analysis ---------------------------------------------------
# Rationale: Quantify direct vs. indirect effects by computing path products
#            from posterior draws (mediational analysis)

# --- Parameter Name Definitions ---
# Rationale: Match coefficient names from brms output (use parnames(fit) to verify)
param_china_native   <- "b_culChinatdwg3scaled_native.global.tdwg3.scaled"
param_foreign_native <- "b_culforeigntdwg3scaled_native.global.tdwg3.scaled"
param_foreign_china  <- "b_culforeigntdwg3scaled_cul.China.tdwg3.scaled"
param_nat_native     <- "b_natextentscaled_native.global.tdwg3.scaled"
param_nat_china      <- "b_natextentscaled_cul.China.tdwg3.scaled"
param_nat_foreign    <- "b_natextentscaled_cul.foreign.tdwg3.scaled"

# --- Function: Extract Path-Specific Posterior Draws ---
# Rationale: Calculate four pathway effects per iteration to preserve uncertainty
extract_path_posteriors <- function(fit, group_name) {
  draws <- as_draws_df(fit)
  
  draws %>%
    transmute(
      Direct               = .data[[param_nat_native]],
      Indirect_via_China   = .data[[param_china_native]] * .data[[param_nat_china]],
      Indirect_via_Foreign = .data[[param_foreign_native]] * .data[[param_nat_foreign]],
      Indirect_via_Both    = .data[[param_china_native]] * .data[[param_foreign_china]] * 
        .data[[param_nat_foreign]]
    ) %>%
    pivot_longer(everything(), names_to = "Path", values_to = "Effect") %>%
    mutate(Group = group_name)
}

# --- Consolidate All Posterior Draws ---
model_list <- list(
  Annual    = bsem_natsucc_annual,
  Perennial = bsem_natsucc_perennial,
  Woody     = bsem_natsucc_woody
)

posteriors_long <- bind_rows(
  lapply(names(model_list), \(g) extract_path_posteriors(model_list[[g]], g))
)

# --- Data Preparation for Visualization ---
# Rationale: Standardize factor levels for consistent ordering in figures
posteriors_long <- posteriors_long %>%
  mutate(
    Group = factor(
      Group,
      levels = c("Annual", "Perennial", "Woody"),
      labels = c("Annual herb.", "Perennial herb.", "Woody")
    ),
    Path = factor(
      Path,
      levels = c("Direct", "Indirect_via_China", "Indirect_via_Foreign", "Indirect_via_Both"),
      labels = c(
        "NRS → NatSuc",
        "NRS → CWC → NatSuc",
        "NRS → COC → NatSuc",
        "NRS → CWC → COC → NatSuc"
      )
    )
  )


# 5. Visualization: Ridge Plots ------------------------------------------------
# Rationale: Ridge plots display full posterior distributions, revealing
#            uncertainty and multi-modality in path effect estimates

# --- Facet Label Formatting ---
path_levels <- levels(posteriors_long$Path)
facet_labels <- paste0("(", LETTERS[seq_along(path_levels)], ") ", path_levels)
names(facet_labels) <- path_levels

# --- Publication Theme ---
# Rationale: Nature/Science journals require high-contrast, serif fonts,
#            and minimal non-data ink
theme_publication <- theme_classic(base_size = 16, base_family = "serif") +
  theme(
    panel.background   = element_rect(fill = "white", colour = NA),
    plot.background    = element_rect(fill = "white", colour = NA),
    panel.grid         = element_blank(),
    axis.line          = element_blank(),
    panel.border       = element_rect(colour = "black", fill = NA, linewidth = 1.5),
    axis.ticks         = element_line(colour = "black", linewidth = 1.2),
    axis.text          = element_text(colour = "black", size = 16),
    axis.title         = element_text(colour = "black", size = 18),
    legend.title       = element_text(size = 16),
    legend.text        = element_text(size = 16),
    legend.position    = "top",
    strip.background   = element_blank(),
    strip.text         = element_text(face = "bold", size = 18, hjust = 0)
  )

# --- Generate Ridge Plot ---
fig_ridge <- ggplot(
  posteriors_long,
  aes(x = Effect, y = Group, fill = Group)
) +
  geom_density_ridges(
    alpha          = 0.7,
    scale          = 1.1,
    rel_min_height = 0.01,
    color          = "grey30"
  ) +
  geom_vline(
    xintercept = 0, 
    linetype   = "dashed", 
    linewidth  = 0.5, 
    color      = "grey40"
  ) +
  facet_wrap(
    ~ Path, 
    scales   = "free_x", 
    ncol     = 1, 
    labeller = labeller(Path = facet_labels)
  ) +
  scale_fill_manual(
    name   = NULL,
    values = c(
      "Annual herb."    = "#E31A1C",
      "Perennial herb." = "#1F78B4",
      "Woody"           = "#F8AF66"
    ),
    breaks = c("Annual herb.", "Perennial herb.", "Woody")
  ) +
  labs(x = "Posterior effect size", y = NULL) +
  theme_publication

print(fig_ridge)

# --- Export Figure ---
# Rationale: High-resolution PNG (300 DPI) suitable for journal submission
ggexport(
  fig_ridge,
  filename  = here("results", "Figure_S4_pathway_effects.png"),
  width     = 1800,
  height    = 4200,
  pointsize = 12,
  res       = 300
)


# 6. Hypothesis Testing: Directional Effects -----------------------------------
# Rationale: Test one-tailed hypotheses (β > 0 or β < 0) for each path
#            coefficient using Bayesian evidence ratios (Evid.Ratio)

# --- Function: Extract and Test Hypotheses ---
# Rationale: Automate hypothesis generation based on coefficient sign
test_model_hypotheses <- function(model_fit, model_name) {
  # Tidy model output
  tidy_output <- tidy(model_fit) %>%
    mutate(
      direction  = ifelse(estimate > 0, " > 0", " < 0"),
      hyp_code   = paste0(response, "_", term, direction),
      Hypothesis = paste0("(", response, "_", term, ")", direction)
    )
  
  # Filter fixed effects (exclude intercepts)
  fixed_effects <- tidy_output %>%
    filter(effect == "fixed", term != "(Intercept)")
  
  # Test each hypothesis
  hyp_tests <- fixed_effects %>%
    mutate(
      test_result = purrr::map(hyp_code, ~ hypothesis(model_fit, .x, seed = 2024)),
      hypothesis_df = purrr::map(test_result, "hypothesis")
    )
  
  # Extract hypothesis test results
  hyp_results <- hyp_tests %>%
    pull(hypothesis_df) %>%
    purrr::map_dfr(~.x)
  
  # Join results with tidy output
  final_table <- tidy_output %>%
    left_join(hyp_results, by = "Hypothesis") %>%
    mutate(Model = model_name)
  
  return(final_table)
}

# --- Test All Models ---
hyp_all       <- test_model_hypotheses(bsem_natsucc_all,       "All species")
hyp_annual    <- test_model_hypotheses(bsem_natsucc_annual,    "Annual herbs")
hyp_perennial <- test_model_hypotheses(bsem_natsucc_perennial, "Perennial herbs")
hyp_woody     <- test_model_hypotheses(bsem_natsucc_woody,     "Woody plants")

# --- Consolidate Results ---
hypothesis_table <- bind_rows(hyp_all, hyp_annual, hyp_perennial, hyp_woody)

# --- Export Summary Table ---
# Rationale: Provides Supplement Table with Bayesian credible intervals and p-directional
write_csv(
  hypothesis_table,
  here("results", "Table_S2_hypothesis_tests.csv")
)

# ==============================================================================
# END OF SCRIPT
# ==============================================================================