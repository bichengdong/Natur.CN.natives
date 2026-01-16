# ==============================================================================
# Manuscript:   Cultivation facilitates global naturalization
# Description:  Two-part negative binomial hurdle models examining how native
#               range size, cultivation outside China, and cultivation within
#               China predict naturalization extent across life forms.
#               Zero component models presence/absence of naturalization;
#               count component models naturalization breadth when present.
# Author:       YW and DBC
# Date:         2026-01-16
# ==============================================================================

# 1. Environment Setup ---------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  here,         # Reproducible file paths
  tidyverse,    # Data manipulation and ggplot2
  readr,        # CSV import
  patchwork,    # Multi-panel figure composition
  ggpubr,       # Publication-ready export (ggexport)
  pscl,         # Hurdle models (hurdle function)
  stringr       # String operations (if needed for future extensions)
)

# Clear workspace to ensure reproducibility
rm(list = ls())
cat("\014")

# 2. Data Loading --------------------------------------------------------------
# Rationale: Import matched native-cultivated plant dataset for hurdle modeling
flora_raw <- read_csv(
  here("data", "20250530.native_plant.matched.cultivated_plant_DO11.csv"),
  show_col_types = FALSE
)

# 3. Data Processing -----------------------------------------------------------
# Hypothesis: Naturalization extent (nat.extent) is jointly determined by:
#             1) Cultivation pressure (inside/outside China)
#             2) Native range size (propagule pressure proxy)
#             3) Life form (ecological strategy)

# --- Extract core variables ---
# Rationale: Focus on response (nat.extent) and predictors (ranges, life form)
flora_core <- flora_raw %>%
  select(
    taxon_name,
    nat.extent,                  # Response: naturalization breadth (count)
    planting.China.tdwg3,        # Cultivation range within China
    WorldCuP.n.tdwg3,            # Cultivation range outside China
    native.global.tdwg3,         # Native range size (TDWG Level-3 count)
    life.form.integrated         # Categorical: annual/perennial/woody
  )

# --- Filter and factorize ---
# Rationale: Remove records with missing life form; set factor levels with
#            "woody" as baseline (reference category) for interpretability
flora_analysis <- flora_core %>%
  filter(!is.na(life.form.integrated)) %>%
  mutate(life.form.integrated = factor(
    life.form.integrated,
    levels = c("woody", "annual herb", "perennial herb")  # Woody = reference
  ))

# --- Standardize predictors ---
# Rationale: Log(x+1) transformation addresses:
#            1) Right-skewed distributions of range sizes
#            2) Multiplicative relationships (e.g., doubling range → additive effect)
#            3) +1 handles zeros in cultivation ranges
#            Subsequent scaling (z-scores) ensures comparable coefficient magnitudes
flora_analysis <- flora_analysis %>%
  mutate(
    scale_WorldCuP       = as.numeric(scale(log(WorldCuP.n.tdwg3 + 1))),
    scale_native_range   = as.numeric(scale(log(native.global.tdwg3 + 1))),
    scale_China_cultiv   = as.numeric(scale(log(planting.China.tdwg3 + 1)))
  )

# 4. Statistical Modeling ------------------------------------------------------
# Rationale: Hurdle model decomposes naturalization into:
#            - Zero component (binomial): Probability of ANY naturalization
#            - Count component (negbin): Extent of naturalization GIVEN presence
#            Negative binomial addresses overdispersion in count data

# --- Model specification ---
# Formula structure: response ~ count_predictors | zero_predictors
# Both components include same predictors + life form interactions
hurdle_model <- hurdle(
  nat.extent ~ 
    scale_native_range + scale_WorldCuP + scale_China_cultiv +
    life.form.integrated +
    scale_native_range:life.form.integrated +
    scale_WorldCuP:life.form.integrated +
    scale_China_cultiv:life.form.integrated |
    scale_native_range + scale_WorldCuP + scale_China_cultiv +
    life.form.integrated +
    scale_native_range:life.form.integrated +
    scale_WorldCuP:life.form.integrated +
    scale_China_cultiv:life.form.integrated,
  dist = "negbin",
  data = flora_analysis
)

# --- Reproducibility: Save model object ---
# Rationale: Hurdle models with interactions are computationally expensive;
#            saving fitted object enables rapid re-plotting without re-fitting
if (!dir.exists(here("results"))) dir.create(here("results"))
saveRDS(hurdle_model, here("results", "hurdle_model_naturalization.rds"))

model_summary <- summary(hurdle_model)

# 5. Visualization Theme -------------------------------------------------------
# Rationale: Consistent publication-grade aesthetics across all panels
theme_publication <- theme_classic(base_size = 14, base_family = "serif") +
  theme(
    panel.background   = element_rect(fill = "white", colour = NA),
    plot.background    = element_rect(fill = "white", colour = NA),
    panel.grid         = element_blank(),
    axis.line          = element_line(colour = "black", linewidth = 1.5),
    axis.ticks         = element_line(colour = "black", linewidth = 1.5),
    axis.text          = element_text(colour = "black", size = 14),
    axis.title         = element_text(colour = "black", size = 14),
    legend.title       = element_text(size = 14),
    legend.text        = element_text(size = 14),
    legend.position    = "top"
  )

# Define life form colors (consistent with Figure S1)
lifeform_colors <- c(
  "woody"          = "#F8AF66",
  "perennial herb" = "#1F78B4",
  "annual herb"    = "#E31A1C"
)

# 6. Prediction Function -------------------------------------------------------
# Rationale: Reduce code duplication by creating reusable function for
#            generating linear predictor plots across different focal variables

create_linear_predictor_plots <- function(
    focal_var,           # Name of focal predictor (string)
    focal_label,         # X-axis label for plots
    ylim_zero,           # Y-axis limits for zero component
    ylim_count,          # Y-axis limits for count component
    ybreaks_zero,        # Y-axis breaks for zero component
    ybreaks_count        # Y-axis breaks for count component
) {
  
  # --- Generate prediction grid ---
  focal_seq <- seq(
    min(flora_analysis[[focal_var]], na.rm = TRUE),
    max(flora_analysis[[focal_var]], na.rm = TRUE),
    length.out = 100
  )
  
  pred_grid <- expand.grid(
    focal_value          = focal_seq,
    life.form.integrated = levels(flora_analysis$life.form.integrated)
  )
  names(pred_grid)[1] <- focal_var
  
  # --- Hold other predictors at mean ---
  other_vars <- c("scale_native_range", "scale_WorldCuP", "scale_China_cultiv")
  other_vars <- setdiff(other_vars, focal_var)
  
  for (var in other_vars) {
    pred_grid[[var]] <- mean(flora_analysis[[var]], na.rm = TRUE)
  }
  
  # --- Compute linear predictors ---
  X_zero  <- model.matrix(
    delete.response(terms(hurdle_model, component = "zero")),
    data = pred_grid
  )
  X_count <- model.matrix(
    delete.response(terms(hurdle_model, component = "count")),
    data = pred_grid
  )
  
  pred_grid$linear_zero  <- as.vector(X_zero %*% coef(hurdle_model, "zero"))
  pred_grid$linear_count <- as.vector(X_count %*% coef(hurdle_model, "count"))
  
  # --- Sort for smooth lines ---
  pred_grid <- pred_grid %>%
    arrange(life.form.integrated, !!sym(focal_var))
  
  # --- Plot: Zero component ---
  plot_zero <- ggplot(
    pred_grid,
    aes(x = !!sym(focal_var), y = linear_zero, color = life.form.integrated)
  ) +
    geom_line(size = 1.5) +
    scale_color_manual(values = lifeform_colors, name = "Life forms:") +
    labs(
      x = focal_label,
      y = "Linear predictor\n(zero component)"
    ) +
    scale_x_continuous(limits = c(-1.2, 4.5), breaks = seq(-1, 4, by = 1)) +
    scale_y_continuous(limits = ylim_zero, breaks = ybreaks_zero) +
    theme_publication +
    theme(
      axis.title   = element_text(size = 24),
      axis.text    = element_text(size = 24),
      legend.title = element_text(size = 24),
      legend.text  = element_text(size = 24),
      legend.position = "none"
    )
  
  # --- Plot: Count component ---
  plot_count <- ggplot(
    pred_grid,
    aes(x = !!sym(focal_var), y = linear_count, color = life.form.integrated)
  ) +
    geom_line(size = 1.5) +
    scale_color_manual(values = lifeform_colors, name = "Life forms:") +
    labs(
      x = focal_label,
      y = "Linear predictor\n(count component)"
    ) +
    scale_x_continuous(limits = c(-1.2, 4.5), breaks = seq(-1, 4, by = 1)) +
    scale_y_continuous(limits = ylim_count, breaks = ybreaks_count) +
    theme_publication +
    theme(
      axis.title   = element_text(size = 24),
      axis.text    = element_text(size = 24),
      legend.title = element_text(size = 24),
      legend.text  = element_text(size = 24),
      legend.position = "none"
    )
  
  list(zero = plot_zero, count = plot_count)
}

# 7. Generate Plots for Each Predictor -----------------------------------------
# Hypothesis visualization: Separate effects of native range, cultivation inside
#                           China, and cultivation outside China on naturalization

# --- Native range size ---
plots_native <- create_linear_predictor_plots(
  focal_var     = "scale_native_range",
  focal_label   = "Log(native range size + 1) (scaled)",
  ylim_zero     = c(-6.5, 4),
  ylim_count    = c(-2.5, 3),
  ybreaks_zero  = seq(-6, 4, by = 2),
  ybreaks_count = seq(-2, 3, by = 1)
)

# --- Cultivation outside China ---
plots_worldcup <- create_linear_predictor_plots(
  focal_var     = "scale_WorldCuP",
  focal_label   = "Log(cultivation outside China + 1) (scaled)",
  ylim_zero     = c(-6.5, 4),
  ylim_count    = c(-4, 6),
  ybreaks_zero  = seq(-6, 4, by = 2),
  ybreaks_count = seq(-4, 6, by = 2)
)

# --- Cultivation within China ---
plots_china <- create_linear_predictor_plots(
  focal_var     = "scale_China_cultiv",
  focal_label   = "Log(cultivation within China + 1) (scaled)",
  ylim_zero     = c(-6, -1),
  ylim_count    = c(-3, 2),
  ybreaks_zero  = seq(-6, -1, by = 1),
  ybreaks_count = seq(-3, 2, by = 1)
)

# 8. Compose Multi-Panel Figure ------------------------------------------------
# Rationale: 6-panel layout (3 predictors × 2 components) with shared legend
#            facilitates comparison of effect patterns across variables

figure_s3 <- (
  plots_native$zero + plots_native$count +
    plots_china$zero + plots_china$count +
    plots_worldcup$zero + plots_worldcup$count +
    plot_layout(ncol = 2, guides = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 24))
) &
  theme(legend.position = "top") &
  scale_color_manual(
    values = c(
      "annual herb"    = "#E31A1C",
      "perennial herb" = "#1F78B4",
      "woody"          = "#F8AF66"
    ),
    name   = "Life forms:",
    labels = c(
      "annual herb"    = "Annual herb.",
      "perennial herb" = "Perennial herb.",
      "woody"          = "Woody"
    ),
    breaks = c("annual herb", "perennial herb", "woody")
  )

# 9. Export --------------------------------------------------------------------
# Rationale: High-resolution PNG for supplementary materials
ggexport(
  figure_s3,
  filename  = here("results", "Figure_S3.png"),
  width     = 4200,
  height    = 6000,
  pointsize = 12,
  res       = 300
)


# ==============================================================================
# End of Script
# ==============================================================================