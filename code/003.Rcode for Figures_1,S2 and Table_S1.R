# ==============================================================================
# Manuscript:   Cultivation facilitates global naturalization
# Description:  Hurdle models (binomial + zero-truncated negative binomial) 
#               analyzing how native range size, cultivation extent (within/
#               outside China), and life form predict plant naturalization 
#               incidence and extent in China.
# Author:       YW and DBC
# Date:         2026-01-16
# ==============================================================================

# 1. Environment Setup ---------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  pscl,
  here,        # Relative path management
  tidyverse,   # Data manipulation and visualization (includes dplyr, ggplot2, readr)
  patchwork,   # Combine multiple plots
  ggpubr,      # Export publication-quality plots
  VGAM,        # Zero-truncated negative binomial models
  broom        # Tidy model outputs
)

# Clear workspace to ensure reproducibility
rm(list = ls())
cat("\014")

# 2. Data Loading & Preprocessing ----------------------------------------------
# Rationale: Load matched native-cultivated plant dataset with standardized
#            TDWG Level 3 geographic units and trait data
native_flora <- read_csv(
  here("data", "20250530.native_plant.matched.cultivated_plant_DO11.csv"),
  show_col_types = FALSE
)

# Rationale: Select key variables for analysis (naturalization extent, native/
# cultivation ranges, life form) and filter out missing life form records
native_flora_clean <- native_flora %>%
  select(
    taxon_name,
    nat.extent,                  # Response: naturalization extent (count of TDWG3 regions)
    planting.China.tdwg3,        # Predictor: cultivation extent within China
    WorldCuP.n.tdwg3,            # Predictor: cultivation extent outside China
    native.global.tdwg3,         # Predictor: native range size
    life.form.integrated         # Moderator: life form category
  ) %>%
  filter(!is.na(life.form.integrated)) %>%
  mutate(
    # Rationale: Set "woody" as reference level to compare herbs against woody plants
    life.form.integrated = factor(
      life.form.integrated,
      levels = c("woody", "annual herb", "perennial herb")
    )
  )

# Rationale: Standardize predictors using log-transformation (to handle skewness)
# and z-score scaling (to enable effect size comparison across variables)
native_flora_scaled <- native_flora_clean %>%
  mutate(
    scale.WorldCuP.n.tdwg3       = as.numeric(scale(log(WorldCuP.n.tdwg3 + 1))),
    scale.native.global.tdwg3    = as.numeric(scale(log(native.global.tdwg3 + 1))),
    scale.planting.China.tdwg3   = as.numeric(scale(log(planting.China.tdwg3 + 1)))
  )

# Inspect data structure
str(native_flora_scaled)

# 3. Statistical Modeling ------------------------------------------------------

# --- 3.1 Full Hurdle Model ---
# Rationale: Hurdle models address zero-inflation by separately modeling:
#   (1) Binomial part: probability of naturalization occurrence (yes/no)
#   (2) Count part: naturalization extent (number of regions, given occurrence)
# Hypothesis: Native range size and cultivation extent interact with life form
# to determine both naturalization probability and extent.

hurdle_full <- hurdle(
  nat.extent ~ 
    scale.native.global.tdwg3 +
    scale.WorldCuP.n.tdwg3 +
    scale.planting.China.tdwg3 +
    life.form.integrated +
    scale.native.global.tdwg3:life.form.integrated +
    scale.WorldCuP.n.tdwg3:life.form.integrated +
    scale.planting.China.tdwg3:life.form.integrated |
    scale.native.global.tdwg3 +
    scale.WorldCuP.n.tdwg3 +
    scale.planting.China.tdwg3 +
    life.form.integrated +
    scale.native.global.tdwg3:life.form.integrated +
    scale.WorldCuP.n.tdwg3:life.form.integrated +
    scale.planting.China.tdwg3:life.form.integrated,
  dist = "negbin",
  data = native_flora_scaled
)

summary(hurdle_full)

# --- Extract and Save Hurdle Model Coefficients ---
result_count <- as.data.frame(summary(hurdle_full)$coefficients$count) %>%
  tibble::rownames_to_column(var = "term") %>%
  mutate(part = "count")

result_zero <- as.data.frame(summary(hurdle_full)$coefficients$zero) %>%
  tibble::rownames_to_column(var = "term") %>%
  mutate(part = "zero")

result_hurdle <- bind_rows(result_count, result_zero)
write_csv(result_hurdle, here("results", "hurdle_full_coefficients.csv"))

# --- 3.2 Separate Binomial and Zero-Truncated Models ---
# Rationale: Fit component models separately for detailed inspection and plotting

model_formula <- as.formula(
  nat.extent ~ 
    scale.native.global.tdwg3 + 
    scale.WorldCuP.n.tdwg3 + 
    scale.planting.China.tdwg3 + 
    life.form.integrated +
    scale.native.global.tdwg3:life.form.integrated +
    scale.WorldCuP.n.tdwg3:life.form.integrated +
    scale.planting.China.tdwg3:life.form.integrated
)

# Binomial model: Probability of naturalization (yes/no)
model_binomial <- glm(
  model_formula,
  data = native_flora_scaled %>% mutate(nat.extent = as.numeric(nat.extent > 0)),
  family = binomial
)
summary(model_binomial)

# Zero-truncated negative binomial: Extent of naturalization (given occurrence)
data_nonzero <- subset(native_flora_scaled, nat.extent > 0)
model_truncated_nb <- vglm(
  model_formula,
  data = data_nonzero,
  family = posnegbinomial
)
summary(model_truncated_nb)

# --- Extract and Combine Coefficients from Component Models ---
tidy_binomial <- tidy(model_binomial) %>%
  mutate(model = "Binomial (zero part)")

coef_mat_truncated <- coef(summary(model_truncated_nb))
tidy_truncated <- as.data.frame(coef_mat_truncated) %>%
  tibble::rownames_to_column("term") %>%
  rename(
    estimate   = Estimate,
    std.error  = `Std. Error`,
    statistic  = `z value`,
    p.value    = `Pr(>|z|)`
  ) %>%
  mutate(model = "Zero-truncated NB (count part)")

tidy_combined <- bind_rows(tidy_binomial, tidy_truncated) %>%
  select(model, term, estimate, std.error, statistic, p.value)

write_csv(tidy_combined, here("results", "component_models_coefficients.csv"))

# 4. Model Predictions for Visualization --------------------------------------

# --- 4.1 Define Prediction Sequences for Focal Variables ---
native_seq <- seq(
  min(native_flora_scaled$scale.native.global.tdwg3, na.rm = TRUE),
  max(native_flora_scaled$scale.native.global.tdwg3, na.rm = TRUE),
  length.out = 100
)

worldcup_seq <- seq(
  min(native_flora_scaled$scale.WorldCuP.n.tdwg3, na.rm = TRUE),
  max(native_flora_scaled$scale.WorldCuP.n.tdwg3, na.rm = TRUE),
  length.out = 100
)

china_seq <- seq(
  min(native_flora_scaled$scale.planting.China.tdwg3, na.rm = TRUE),
  max(native_flora_scaled$scale.planting.China.tdwg3, na.rm = TRUE),
  length.out = 100
)

life_levels <- levels(native_flora_scaled$life.form.integrated)

# --- 4.2 Prediction Function for Hurdle Model ---
# Rationale: Generate predicted naturalization extent across focal variable gradients
predict_hurdle_effect <- function(model, var_name, seq_values, data, life_levels) {
  
  pred_data <- expand.grid(
    life.form.integrated = life_levels,
    x_seq = seq_values
  )
  
  # Hold non-focal predictors at their mean values
  pred_data$scale.native.global.tdwg3  <- mean(data$scale.native.global.tdwg3, na.rm = TRUE)
  pred_data$scale.WorldCuP.n.tdwg3     <- mean(data$scale.WorldCuP.n.tdwg3, na.rm = TRUE)
  pred_data$scale.planting.China.tdwg3 <- mean(data$scale.planting.China.tdwg3, na.rm = TRUE)
  
  # Replace focal variable with sequence
  pred_data[[var_name]] <- pred_data$x_seq
  pred_data$x_seq <- NULL
  
  # Predict expected naturalization extent
  pred_data$response <- predict(model, newdata = pred_data, type = "response")
  
  pred_data <- pred_data %>%
    mutate(
      life.form.integrated = factor(
        life.form.integrated,
        levels = c("annual herb", "perennial herb", "woody")
      )
    )
  
  return(pred_data)
}

# Generate predictions for hurdle model
pred_hurdle_native   <- predict_hurdle_effect(hurdle_full, "scale.native.global.tdwg3", native_seq, native_flora_scaled, life_levels)
pred_hurdle_worldcup <- predict_hurdle_effect(hurdle_full, "scale.WorldCuP.n.tdwg3", worldcup_seq, native_flora_scaled, life_levels)
pred_hurdle_china    <- predict_hurdle_effect(hurdle_full, "scale.planting.China.tdwg3", china_seq, native_flora_scaled, life_levels)

# --- 4.3 Prediction Function for Binomial Model ---
predict_binomial_effect <- function(model, var_name, seq_values, data, life_levels) {
  
  pred_data <- expand.grid(
    life.form.integrated = life_levels,
    x_seq = seq_values
  )
  
  pred_data$scale.native.global.tdwg3  <- mean(data$scale.native.global.tdwg3, na.rm = TRUE)
  pred_data$scale.WorldCuP.n.tdwg3     <- mean(data$scale.WorldCuP.n.tdwg3, na.rm = TRUE)
  pred_data$scale.planting.China.tdwg3 <- mean(data$scale.planting.China.tdwg3, na.rm = TRUE)
  
  pred_data[[var_name]] <- pred_data$x_seq
  pred_data$x_seq <- NULL
  
  pred_data$prob_naturalized <- predict(model, newdata = pred_data, type = "response")
  
  pred_data <- pred_data %>%
    mutate(
      life.form.integrated = factor(
        life.form.integrated,
        levels = c("annual herb", "perennial herb", "woody")
      )
    )
  
  return(pred_data)
}

pred_binom_native   <- predict_binomial_effect(model_binomial, "scale.native.global.tdwg3", native_seq, native_flora_scaled, life_levels)
pred_binom_worldcup <- predict_binomial_effect(model_binomial, "scale.WorldCuP.n.tdwg3", worldcup_seq, native_flora_scaled, life_levels)
pred_binom_china    <- predict_binomial_effect(model_binomial, "scale.planting.China.tdwg3", china_seq, native_flora_scaled, life_levels)

# --- 4.4 Prediction Function for Zero-Truncated Negative Binomial Model ---
predict_truncnegbin_effect <- function(model, var_name, seq_values, data, life_levels) {
  
  pred_data <- expand.grid(
    life.form.integrated = life_levels,
    x_seq = seq_values
  )
  
  pred_data$scale.native.global.tdwg3  <- mean(data$scale.native.global.tdwg3, na.rm = TRUE)
  pred_data$scale.WorldCuP.n.tdwg3     <- mean(data$scale.WorldCuP.n.tdwg3, na.rm = TRUE)
  pred_data$scale.planting.China.tdwg3 <- mean(data$scale.planting.China.tdwg3, na.rm = TRUE)
  
  pred_data[[var_name]] <- pred_data$x_seq
  pred_data$x_seq <- NULL
  
  pred_data$expected_count <- predict(model, newdata = pred_data, type = "response")
  
  pred_data <- pred_data %>%
    mutate(
      life.form.integrated = factor(
        life.form.integrated,
        levels = c("annual herb", "perennial herb", "woody")
      )
    )
  
  return(pred_data)
}

pred_trunc_native   <- predict_truncnegbin_effect(model_truncated_nb, "scale.native.global.tdwg3", native_seq, data_nonzero, life_levels)
pred_trunc_worldcup <- predict_truncnegbin_effect(model_truncated_nb, "scale.WorldCuP.n.tdwg3", worldcup_seq, data_nonzero, life_levels)
pred_trunc_china    <- predict_truncnegbin_effect(model_truncated_nb, "scale.planting.China.tdwg3", china_seq, data_nonzero, life_levels)

# --- 4.5 Prepare Observed Data for Plotting ---
observed_data <- native_flora_scaled %>%
  select(
    scale.native.global.tdwg3,
    scale.WorldCuP.n.tdwg3,
    scale.planting.China.tdwg3,
    life.form.integrated,
    nat.extent
  ) %>%
  mutate(
    observed_count  = nat.extent,
    observed_binary = as.numeric(nat.extent > 0),
    life.form.integrated = factor(
      life.form.integrated,
      levels = c("annual herb", "perennial herb", "woody")
    )
  )

# 5. Visualization Setup -------------------------------------------------------

# --- 5.1 Define Global Plot Theme ---
# Rationale: Consistent publication-quality aesthetics across all figures
theme_standard <- theme_classic(base_size = 14, base_family = "serif") +
  theme(
    panel.background  = element_rect(fill = "white", colour = NA),
    plot.background   = element_rect(fill = "white", colour = NA),
    panel.grid        = element_blank(),
    axis.line         = element_line(colour = "black", linewidth = 1.5),
    axis.ticks        = element_line(colour = "black", linewidth = 1.5),
    axis.text         = element_text(colour = "black", size = 14),
    axis.title        = element_text(colour = "black", size = 14),
    legend.title      = element_text(size = 14),
    legend.text       = element_text(size = 14),
    legend.position   = "top"
  )

# Define color palette for life forms
color_palette <- c(
  "woody"           = "#F8AF66",
  "perennial herb"  = "#1F78B4",
  "annual herb"     = "#E31A1C"
)

# --- 5.2 Plotting Functions ---

# Hurdle model (full response): Naturalization extent
plot_hurdle_response <- function(pred_data, observed_data, x_var, x_label) {
  
  ggplot(pred_data, aes(x = .data[[x_var]], y = response, color = life.form.integrated)) +
    geom_point(
      data  = observed_data,
      aes(x = .data[[x_var]], y = observed_count),
      alpha = 0.4,
      size  = 4
    ) +
    geom_line(size = 1.5) +
    scale_color_manual(values = color_palette, name = "Life form:") +
    labs(x = x_label, y = "Naturalization success\n(no. of regions)") +
    scale_x_continuous(limits = c(-1.2, 4.5), breaks = seq(-1, 4, by = 1)) +
    scale_y_continuous(limits = c(-5, 200), breaks = seq(0, 200, by = 50)) +
    theme_standard +
    theme(
      axis.title.x  = element_text(size = 26),
      axis.title.y  = element_text(size = 26),
      axis.text     = element_text(size = 26),
      legend.title  = element_text(size = 26),
      legend.text   = element_text(size = 26)
    )
}

# Binomial model: Naturalization incidence (probability)
plot_binomial_response <- function(pred_data, observed_data, x_var, x_label) {
  
  ggplot(pred_data, aes(x = .data[[x_var]], y = prob_naturalized, color = life.form.integrated)) +
    geom_point(
      data     = observed_data,
      aes(x = .data[[x_var]], y = observed_binary),
      alpha    = 0.4,
      size     = 4,
      position = position_jitter(height = 0.05, width = 0, seed = 123)
    ) +
    geom_line(size = 1.5) +
    scale_color_manual(values = color_palette, name = "Life form:") +
    labs(x = x_label, y = "Naturalization incidence\n(no, yes)") +
    scale_x_continuous(limits = c(-1.2, 4.5), breaks = seq(-1, 4, by = 1)) +
    scale_y_continuous(limits = c(-0.05, 1.05), breaks = seq(0, 1), labels = c("0", "1")) +
    theme_standard +
    theme(
      axis.title.x  = element_text(size = 26),
      axis.title.y  = element_text(size = 26),
      axis.text     = element_text(size = 26),
      legend.title  = element_text(size = 26),
      legend.text   = element_text(size = 26)
    )
}

# Zero-truncated negative binomial: Naturalization extent (given occurrence)
plot_truncnegbin_response <- function(pred_data, observed_data, x_var, x_label) {
  
  ggplot(pred_data, aes(x = .data[[x_var]], y = expected_count, color = life.form.integrated)) +
    geom_point(
      data  = filter(observed_data, observed_count > 0),
      aes(x = .data[[x_var]], y = observed_count),
      alpha = 0.4,
      size  = 4
    ) +
    geom_line(size = 1.5) +
    scale_color_manual(values = color_palette, name = "Life form:") +
    labs(x = x_label, y = "Naturalization extent\n(no. of regions)") +
    scale_x_continuous(limits = c(-1.2, 4.5), breaks = seq(-1, 4, by = 1)) +
    scale_y_continuous(limits = c(-5, 200), breaks = seq(0, 200, by = 50)) +
    theme_standard +
    theme(
      axis.title.x  = element_text(size = 26),
      axis.title.y  = element_text(size = 26),
      axis.text     = element_text(size = 26),
      legend.title  = element_text(size = 26),
      legend.text   = element_text(size = 26)
    )
}

# Inset zoom plot for hurdle model
plot_hurdle_inset <- function(pred_data, x_var, y_max = 10, y_breaks = 2) {
  
  ggplot(pred_data, aes(x = .data[[x_var]], y = response, color = life.form.integrated)) +
    geom_line(linewidth = 1.5) +
    scale_color_manual(values = color_palette) +
    scale_x_continuous(limits = c(-1, 4.5), breaks = seq(-1, 4, by = 1)) +
    scale_y_continuous(limits = c(0, y_max), breaks = seq(0, y_max, by = y_breaks)) +
    labs(x = "", y = "", tag = "") +
    theme_standard +
    theme(
      legend.position = "none",
      axis.text       = element_text(size = 14),
      title           = element_text(size = 14),
      panel.border    = element_rect(color = "black", fill = NA, linewidth = 1.5)
    )
}

# Inset zoom plot for zero-truncated negative binomial model
plot_truncnegbin_inset <- function(pred_data, x_var, y_max = 10, y_breaks = 2) {
  
  ggplot(pred_data, aes(x = .data[[x_var]], y = expected_count, color = life.form.integrated)) +
    geom_line(linewidth = 1.5) +
    scale_color_manual(values = color_palette) +
    scale_x_continuous(limits = c(-1, 4.5), breaks = seq(-1, 4, by = 1)) +
    scale_y_continuous(limits = c(0, y_max), breaks = seq(0, y_max, by = y_breaks)) +
    labs(x = "", y = "", tag = "") +
    theme_standard +
    theme(
      legend.position = "none",
      axis.text       = element_text(size = 14),
      title           = element_text(size = 14),
      panel.border    = element_rect(color = "black", fill = NA, linewidth = 1.5)
    )
}

# Inset zoom plot for binomial model
plot_binomial_inset <- function(pred_data, x_var, y_max = 0.1, y_breaks = 0.1) {
  
  ggplot(pred_data, aes(x = .data[[x_var]], y = prob_naturalized, color = life.form.integrated)) +
    geom_line(size = 1.5) +
    scale_color_manual(values = color_palette) +
    scale_x_continuous(limits = c(-1, 4.5), breaks = seq(-1, 4, by = 1)) +
    scale_y_continuous(limits = c(0, y_max), breaks = c(0, y_max), labels = c("0", as.character(y_max))) +
    labs(x = "", y = "", tag = "") +
    theme_standard +
    theme(
      legend.position = "none",
      axis.text       = element_text(size = 14),
      title           = element_text(size = 14),
      panel.border    = element_rect(color = "black", fill = NA, linewidth = 1.5)
    )
}

# 6. Generate Publication Figures ----------------------------------------------

# --- 6.1 Hurdle Model: Full Response Plots with Insets ---

plot_hurdle_native <- plot_hurdle_response(
  pred_hurdle_native,
  observed_data,
  "scale.native.global.tdwg3",
  "Log(native range size + 1) (scaled)"
)

plot_hurdle_china <- plot_hurdle_response(
  pred_hurdle_china,
  observed_data,
  "scale.planting.China.tdwg3",
  "Log(cultivation within China + 1) (scaled)"
)

plot_hurdle_worldcup <- plot_hurdle_response(
  pred_hurdle_worldcup,
  observed_data,
  "scale.WorldCuP.n.tdwg3",
  "Log(cultivation outside China + 1) (scaled)"
)

# Create inset plots
inset_hurdle_native   <- ggplotGrob(plot_hurdle_inset(pred_hurdle_native, "scale.native.global.tdwg3"))
inset_hurdle_china    <- ggplotGrob(plot_hurdle_inset(pred_hurdle_china, "scale.planting.China.tdwg3", y_max = 0.5, y_breaks = 0.1))
inset_hurdle_worldcup <- ggplotGrob(plot_hurdle_inset(pred_hurdle_worldcup, "scale.WorldCuP.n.tdwg3"))

# Add insets to main plots
plot_hurdle_native_final <- plot_hurdle_native +
  annotation_custom(grob = inset_hurdle_native, xmin = -1, xmax = 1, ymin = 120, ymax = 220)

plot_hurdle_china_final <- plot_hurdle_china +
  annotation_custom(grob = inset_hurdle_china, xmin = -1, xmax = 1, ymin = 120, ymax = 220)

plot_hurdle_worldcup_final <- plot_hurdle_worldcup +
  annotation_custom(grob = inset_hurdle_worldcup, xmin = -1, xmax = 1, ymin = 120, ymax = 220)

# Combine hurdle plots
figure_hurdle <- (
  plot_hurdle_native_final +
    plot_hurdle_china_final +
    plot_hurdle_worldcup_final
) +
  plot_layout(ncol = 1, guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "top") &
  scale_color_manual(
    values = c(
      "annual herb"     = "#E31A1C",
      "perennial herb"  = "#1F78B4",
      "woody"           = "#F8AF66"
    ),
    name = "Life forms:",
    labels = c(
      "annual herb"    = "Annual herb.",
      "perennial herb" = "Perennial herb.",
      "woody"          = "Woody"
    )
  )

# Save hurdle model figure
ggexport(
  figure_hurdle,
  filename  = here("results", "figure_hurdle_full_response.png"),
  width     = 2500,
  height    = 6000,
  pointsize = 12,
  res       = 300
)

# --- 6.2 Component Models: Binomial and Zero-Truncated NB with Insets ---

# Binomial plots
plot_binom_native <- plot_binomial_response(
  pred_binom_native,
  observed_data,
  "scale.native.global.tdwg3",
  "Log(native range size + 1) (scaled)"
)

plot_binom_china <- plot_binomial_response(
  pred_binom_china,
  observed_data,
  "scale.planting.China.tdwg3",
  "Log(cultivation within China + 1) (scaled)"
)

plot_binom_worldcup <- plot_binomial_response(
  pred_binom_worldcup,
  observed_data,
  "scale.WorldCuP.n.tdwg3",
  "Log(cultivation outside China + 1) (scaled)"
)

# Truncated negative binomial plots
plot_trunc_native <- plot_truncnegbin_response(
  pred_trunc_native,
  observed_data,
  "scale.native.global.tdwg3",
  "Log(native range size + 1) (scaled)"
)

plot_trunc_china <- plot_truncnegbin_response(
  pred_trunc_china,
  observed_data,
  "scale.planting.China.tdwg3",
  "Log(cultivation within China + 1) (scaled)"
)

plot_trunc_worldcup <- plot_truncnegbin_response(
  pred_trunc_worldcup,
  observed_data,
  "scale.WorldCuP.n.tdwg3",
  "Log(cultivation outside China + 1) (scaled)"
)

# Create inset plots for binomial model
inset_binom_china <- ggplotGrob(
  plot_binomial_inset(pred_binom_china, "scale.planting.China.tdwg3", y_max = 0.1)
)

plot_binom_china_final <- plot_binom_china +
  annotation_custom(grob = inset_binom_china, xmin = -1.5, xmax = 0.6, ymin = 0.4, ymax = 0.95)

# Create inset plots for truncated negative binomial model
inset_trunc_native   <- ggplotGrob(plot_truncnegbin_inset(pred_trunc_native, "scale.native.global.tdwg3", y_max = 20, y_breaks = 4))
inset_trunc_worldcup <- ggplotGrob(plot_truncnegbin_inset(pred_trunc_worldcup, "scale.WorldCuP.n.tdwg3", y_max = 20, y_breaks = 4))
inset_trunc_china    <- ggplotGrob(plot_truncnegbin_inset(pred_trunc_china, "scale.planting.China.tdwg3", y_max = 30, y_breaks = 6))

plot_trunc_native_final <- plot_trunc_native +
  annotation_custom(grob = inset_trunc_native, xmin = -1, xmax = 1, ymin = 120, ymax = 220)

plot_trunc_china_final <- plot_trunc_china +
  annotation_custom(grob = inset_trunc_china, xmin = -1, xmax = 1, ymin = 120, ymax = 220)

plot_trunc_worldcup_final <- plot_trunc_worldcup +
  annotation_custom(grob = inset_trunc_worldcup, xmin = -1, xmax = 1, ymin = 120, ymax = 220)

# Combine component model plots
figure_components <- (
    plot_binom_native +
    plot_trunc_native_final +
    plot_binom_china_final +
    plot_trunc_china_final +
    plot_binom_worldcup +
    plot_trunc_worldcup_final
) +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "top") &
  scale_color_manual(
    values = c(
      "annual herb"     = "#E31A1C",
      "perennial herb"  = "#1F78B4",
      "woody"           = "#F8AF66"
    ),
    name = "Life forms:",
    labels = c(
      "annual herb"    = "Annual herb.",
      "perennial herb" = "Perennial herb.",
      "woody"          = "Woody"
    )
  )

# Save component models figure
ggexport(
  figure_components,
  filename  = here("results", "figure_component_models.png"),
  width     = 4200,
  height    = 5500,
  pointsize = 12,
  res       = 300
)

# ==============================================================================
# 7. Create Combined 9-Panel Figure (3×3 Grid) --------------------------------
# ==============================================================================

# Rationale: Generate a comprehensive 9-panel figure showing all three models
# (Hurdle full response, Binomial incidence, Zero-truncated extent) across
# all three focal predictors (Native range, Cultivation in China, Cultivation
# outside China) for direct visual comparison in a single publication figure.

cat("\n--- Generating 9-panel combined figure ---\n")

# --- 7.1 Prepare Individual Panels Without Legends ---
# Remove legends from individual plots to avoid redundancy

# Row 1: Hurdle model (Full response - expected naturalization extent)
panel_hurdle_native <- plot_hurdle_response(
  pred_hurdle_native,
  observed_data,
  "scale.native.global.tdwg3",
  "Log(native range size + 1) (scaled)"
) + theme(legend.position = "none")

panel_hurdle_china <- plot_hurdle_response(
  pred_hurdle_china,
  observed_data,
  "scale.planting.China.tdwg3",
  "Log(cultivation within China + 1) (scaled)"
) + theme(legend.position = "none")

panel_hurdle_worldcup <- plot_hurdle_response(
  pred_hurdle_worldcup,
  observed_data,
  "scale.WorldCuP.n.tdwg3",
  "Log(cultivation outside China + 1) (scaled)"
) + theme(legend.position = "none")

# Row 2: Binomial model (Naturalization incidence - probability)
panel_binom_native <- plot_binomial_response(
  pred_binom_native,
  observed_data,
  "scale.native.global.tdwg3",
  "Log(native range size + 1) (scaled)"
) + theme(legend.position = "none")

panel_binom_china <- plot_binomial_response(
  pred_binom_china,
  observed_data,
  "scale.planting.China.tdwg3",
  "Log(cultivation within China + 1) (scaled)"
) + theme(legend.position = "none")

panel_binom_worldcup <- plot_binomial_response(
  pred_binom_worldcup,
  observed_data,
  "scale.WorldCuP.n.tdwg3",
  "Log(cultivation outside China + 1) (scaled)"
) + theme(legend.position = "none")

# Row 3: Zero-truncated NB (Naturalization extent - conditional on occurrence)
panel_trunc_native <- plot_truncnegbin_response(
  pred_trunc_native,
  observed_data,
  "scale.native.global.tdwg3",
  "Log(native range size + 1) (scaled)"
) + theme(legend.position = "none")

panel_trunc_china <- plot_truncnegbin_response(
  pred_trunc_china,
  observed_data,
  "scale.planting.China.tdwg3",
  "Log(cultivation within China + 1) (scaled)"
) + theme(legend.position = "none")

panel_trunc_worldcup <- plot_truncnegbin_response(
  pred_trunc_worldcup,
  observed_data,
  "scale.WorldCuP.n.tdwg3",
  "Log(cultivation outside China + 1) (scaled)"
) + theme(legend.position = "none")

# --- 7.2 Add Insets to Relevant Panels ---

# Row 1 insets (Hurdle model)
inset_9p_hurdle_native   <- ggplotGrob(plot_hurdle_inset(pred_hurdle_native, "scale.native.global.tdwg3", y_max = 10, y_breaks = 2))
inset_9p_hurdle_china    <- ggplotGrob(plot_hurdle_inset(pred_hurdle_china, "scale.planting.China.tdwg3", y_max = 0.5, y_breaks = 0.1))
inset_9p_hurdle_worldcup <- ggplotGrob(plot_hurdle_inset(pred_hurdle_worldcup, "scale.WorldCuP.n.tdwg3", y_max = 10, y_breaks = 2))

panel_hurdle_native_9p <- panel_hurdle_native +
  annotation_custom(grob = inset_9p_hurdle_native, xmin = -1, xmax = 1, ymin = 120, ymax = 220)

panel_hurdle_china_9p <- panel_hurdle_china +
  annotation_custom(grob = inset_9p_hurdle_china, xmin = -1, xmax = 1, ymin = 120, ymax = 220)

panel_hurdle_worldcup_9p <- panel_hurdle_worldcup +
  annotation_custom(grob = inset_9p_hurdle_worldcup, xmin = -1, xmax = 1, ymin = 120, ymax = 220)

# Row 2 insets (Binomial model - only China needs inset)
inset_9p_binom_china <- ggplotGrob(
  plot_binomial_inset(pred_binom_china, "scale.planting.China.tdwg3", y_max = 0.1, y_breaks = 0.1)
)

panel_binom_china_9p <- panel_binom_china +
  annotation_custom(grob = inset_9p_binom_china, xmin = -1.5, xmax = 0.6, ymin = 0.4, ymax = 0.95)

# Row 3 insets (Zero-truncated NB)
inset_9p_trunc_native   <- ggplotGrob(plot_truncnegbin_inset(pred_trunc_native, "scale.native.global.tdwg3", y_max = 20, y_breaks = 4))
inset_9p_trunc_china    <- ggplotGrob(plot_truncnegbin_inset(pred_trunc_china, "scale.planting.China.tdwg3", y_max = 30, y_breaks = 6))
inset_9p_trunc_worldcup <- ggplotGrob(plot_truncnegbin_inset(pred_trunc_worldcup, "scale.WorldCuP.n.tdwg3", y_max = 20, y_breaks = 4))

panel_trunc_native_9p <- panel_trunc_native +
  annotation_custom(grob = inset_9p_trunc_native, xmin = -1, xmax = 1, ymin = 120, ymax = 220)

panel_trunc_china_9p <- panel_trunc_china +
  annotation_custom(grob = inset_9p_trunc_china, xmin = -1, xmax = 1, ymin = 120, ymax = 220)

panel_trunc_worldcup_9p <- panel_trunc_worldcup +
  annotation_custom(grob = inset_9p_trunc_worldcup, xmin = -1, xmax = 1, ymin = 120, ymax = 220)

# --- 7.3 Combine into 3×3 Grid Layout ---
# Layout:
#   Column 1: Native range | Column 2: China cultivation | Column 3: WorldCuP cultivation
#   Row 1: Hurdle (full response)
#   Row 2: Binomial (incidence)
#   Row 3: Zero-truncated NB (extent conditional on occurrence)

figure_9panel <- (
  # Row 1: Hurdle model
  panel_hurdle_native_9p + panel_binom_native + panel_trunc_native_9p + 
  # Row 2: Binomial model
  panel_hurdle_china_9p + panel_binom_china_9p + panel_trunc_china_9p + 
  # Row 3: Zero-truncated NB
  panel_hurdle_worldcup_9p + panel_binom_worldcup + panel_trunc_worldcup_9p
) +
  plot_layout(
    ncol   = 3,
    nrow   = 3,
    guides = "collect"
  ) +
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(plot.tag = element_text(size = 24)) &
  theme(legend.position = "top") &
  scale_color_manual(
    values = c(
      "annual herb"     = "#E31A1C",
      "perennial herb"  = "#1F78B4",
      "woody"           = "#F8AF66"
    ),
    name = "Life forms:",
    labels = c(
      "annual herb"    = "Annual herb.",
      "perennial herb" = "Perennial herb.",
      "woody"          = "Woody"
    )
  )

# --- 7.4 Save 9-Panel Combined Figure ---

# Standard resolution PNG
ggexport(
  figure_9panel,
  filename  = here("results", "figure_9panel_combined.png"),
  width     = 7500,
  height    = 7500,
  pointsize = 12,
  res       = 300
)

# ==============================================================================
# END OF SCRIPT
# ==============================================================================