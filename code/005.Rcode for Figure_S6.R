# ==============================================================================
# Manuscript:   Cultivation facilitates global naturalization
# Description:  Correlation analysis of native range size, global cultivation,
#               and domestic cultivation patterns for Chinese native flora
# Author:       YW and DBC
# Date:         2026-01-16
# ==============================================================================

# 1. Environment Setup ---------------------------------------------------------

# Load required packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  here,      # Relative path management
  dplyr,     # Data manipulation
  ggplot2,   # Visualization
  ggpubr,    # Export publication-quality plots
  tidyverse, # Data science workflow (includes dplyr, ggplot2)
  Hmisc,     # Statistical tools
  GGally     # Pairwise correlation matrices
)

# Clear workspace to ensure reproducibility
cat("\014")
rm(list = ls())

# 2. Data Loading --------------------------------------------------------------
# Rationale: Load dataset of native Chinese plant species with cultivation 
# records to examine relationship between native range and cultivation extent
native_species_raw <- read_csv(
  here("data", "20250530.native_plant.matched.cultivated_plant_DO11.csv"), 
  show_col_types = FALSE
)
str(native_species_raw)

# 3. Data Processing -----------------------------------------------------------
# Extract variables of interest for biogeographic analysis
native_species <- native_species_raw %>%
  select(
    taxon_name,                # Species identity
    nat.extent,                # Native extent metric
    planting.China.tdwg3,      # Cultivation provinces in China (TDWG level 3)
    WorldCuP.n.tdwg3,          # Global cultivation regions outside China
    native.global.tdwg3,       # Native range size (number of TDWG3 regions)
    life.form.integrated       # Growth form classification
  )

# --- Data Validation ---
# Inspect raw data ranges before transformation
range(native_species$WorldCuP.n.tdwg3, na.rm = TRUE)      # Global cultivation
range(native_species$native.global.tdwg3, na.rm = TRUE)   # Native range
range(native_species$planting.China.tdwg3, na.rm = TRUE)  # Domestic cultivation

# --- Normalization for Statistical Modeling ---
# Rationale: Log-transform count data to reduce skewness, then standardize 
# (mean = 0, SD = 1) to make regression coefficients comparable
native_analysis <- native_species %>%
  mutate(
    scale.native.global.tdwg3    = scale(log(native.global.tdwg3 + 1)),
    scale.WorldCuP.n.tdwg3       = scale(log(WorldCuP.n.tdwg3 + 1)),
    scale.planting.China.tdwg3   = scale(log(planting.China.tdwg3 + 1))
  )

# Verify transformation results
str(native_analysis)
glimpse(native_analysis)
range(native_analysis$scale.native.global.tdwg3, na.rm = TRUE)
range(native_analysis$scale.WorldCuP.n.tdwg3, na.rm = TRUE)
range(native_analysis$scale.planting.China.tdwg3, na.rm = TRUE)

# 4. Correlation Analysis ------------------------------------------------------
# Hypothesis: Native range size positively correlates with both global and 
# domestic cultivation extent due to increased propagule pressure and 
# horticultural accessibility

# Prepare correlation matrix dataset
correlation_data <- native_analysis %>%
  select(
    scale.native.global.tdwg3, 
    scale.WorldCuP.n.tdwg3, 
    scale.planting.China.tdwg3
  ) %>%
  mutate_all(~as.numeric(.))

# Check missing data
colSums(is.na(correlation_data))

# --- Custom ggpairs Lower Panel Function ---
# Rationale: Display raw data points with regression line to visualize 
# strength and direction of bivariate relationships
lowerFn <- function(data, mapping, method = "lm", ...) {
  p <- ggplot(data = data, mapping = mapping) +
    geom_point(colour = "grey", alpha = 0.5, size = 3) +
    geom_smooth(method = method, color = "#8B1E3F", ...)
  p
}

# 5. Visualization: Pairwise Correlation Matrix --------------------------------
# Create correlation plot matrix with custom aesthetics
# Reference: https://allisonhorst.github.io/palmerpenguins/articles/intro.html
correlation_plot <- ggpairs(
  correlation_data,
  columns = c(
    "scale.native.global.tdwg3", 
    "scale.WorldCuP.n.tdwg3", 
    "scale.planting.China.tdwg3"
  ),
  upper        = list(continuous = wrap("cor", size = 4, color = "black")),
  lower        = list(continuous = wrap(lowerFn, method = "lm")),
  diag         = list(continuous = wrap("barDiag", 
                                        bins  = 30, 
                                        fill  = "#3F688C",  # NCS S 4040-R90B
                                        color = "#1B3B6F")),
  columnLabels = c(
    "Native range size",
    "Cultivation outside China",
    "Cultivation within China"
  )
)

# --- Apply Publication-Grade Theme ---
correlation_plot_final <- correlation_plot +
  theme_bw() +
  theme(
    panel.background = element_rect(fill = NA),
    axis.text        = element_text(family = "serif", colour = "black", size = 12),
    strip.text       = element_text(family = "serif", colour = "black", size = 12),
    axis.title       = element_text(family = "serif", colour = "black", size = 12)
  ) 

correlation_plot_final

# 6. Export Results ------------------------------------------------------------
# Save high-resolution figure for manuscript supplementary materials
if (!dir.exists(here("results"))) dir.create(here("results"))
ggexport(
  correlation_plot_final, 
  filename  = here("results", "Figure_S6.png"),
  width     = 2000,
  height    = 2000,
  pointsize = 12,
  res       = 300
)

# ==============================================================================
# End of Script
# ==============================================================================