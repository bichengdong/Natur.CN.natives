# ==============================================================================
# Manuscript:   Cultivation facilitates global naturalization
# Description:  Network visualization of SEM pathways linking cultivation history,
#               native range, and naturalization extent across plant life forms
# Author:       YW and DBC
# Date:         2026-01-16
# ==============================================================================

# 1. Environment Setup ---------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  here,        # Project-relative file paths
  tidyverse,   # Data manipulation (includes stringr, readr)
  igraph,      # Graph data structures
  ggraph,      # Network visualization
  tidygraph,   # Tidy interface for graphs
  purrr,       # Functional programming tools
  patchwork,   # Multi-panel plot composition
  ggpubr       # Publication-ready plot export
)

# Clear workspace to ensure reproducibility
cat("\014")
rm(list = ls())

# 2. Custom Functions ----------------------------------------------------------
# Purpose: Generate SEM network diagrams showing causal pathways
# Rationale: Visualize direct/indirect effects of cultivation and native range
#            on plant naturalization success, stratified by life form
visualize_sem_network <- function(data, model_name, response_var) {
  
  # Filter to specific model and prepare categorical color variable
  data <- data %>%
    filter(model == model_name) %>%
    mutate(color = factor(color))
  
  color_palette <- levels(data$color)
  
  # Convert edge list to tbl_graph object
  # Node renaming: Map model terms to human-readable labels for publication
  data.g <- as_tbl_graph(data) %>%
    activate(nodes) %>%
    mutate(name = case_when(
      name == "culChinatdwg3scaled"      ~ "cul.China\n(tdwg3)",
      name == "culforeigntdwg3scaled"    ~ "cul.foreign\n(tdwg3)",
      name == "nativeglobaltdwg3scaled"  ~ "native.global\n(tdwg3)",
      name == "natextentscaled"          ~ "Naturalization extent\n(tdwg3)",
      TRUE                               ~ name
    ))
  
  # Extract nodes and get count for dynamic coordinate assignment
  nodes_base <- data.g %>%
    activate(nodes) %>%
    as_tibble()
  
  n_nodes <- nrow(nodes_base)
  
  # Define node spatial coordinates for consistent layout across panels
  # Rationale: Fixed positions enable visual comparison of path coefficients
  #            between life form groups
  # Layout: Response variable at bottom, predictors arranged above
  if (n_nodes == 4) {
    # Standard 4-node layout: 2 predictors, 1 native range, 1 response
    nodes_df <- tibble(
      name = nodes_base$name,
      x    = c(0, 1, -1, 0),   # Horizontal positions
      y    = c(0.5, 0, 0, -1)  # Vertical positions (response at bottom)
    )
  } else {
    # Fallback for different node counts: arrange in circle
    angles <- seq(0, 2 * pi, length.out = n_nodes + 1)[1:n_nodes]
    nodes_df <- tibble(
      name = nodes_base$name,
      x    = cos(angles),
      y    = sin(angles)
    )
  }
  
  # Construct SEM network diagram
  sem_plot <- ggraph(data.g, layout = nodes_df, x = x, y = y) +
    geom_edge_link(
      aes(
        label      = label,              # Standardized coefficients + significance
        color      = factor(color),      # Red = positive, Blue = negative, Grey = NS
        width      = width,               # Edge thickness proportional to effect size
        start_cap  = circle(20, 'mm'),   # Avoid overlap with node labels
        end_cap    = circle(20, 'mm')
      ),
      angle_calc   = "along",            # Rotate labels to follow edge direction
      label_dodge  = unit(8, "mm"),      # Offset labels from edge
      label_size   = 4,
      arrow        = arrow(length = unit(3, "mm"), type = "closed")
    ) +
    geom_node_point(size = 0, color = "lightblue") +  # Invisible nodes (labels only)
    geom_node_text(
      aes(label = name),
      size   = 5,
      color  = 'black'
    ) +
    scale_edge_color_manual(values = color_palette) +
    theme_graph() +
    theme(
      plot.title      = element_text(hjust = 0.5, size = 16),
      legend.position = "none"
    ) +
    labs(title = model_name)
  
  return(sem_plot)
}

# 3. Data Loading --------------------------------------------------------------
# Source: Supplementary Table S2 from manuscript (Bayesian SEM results)
# Contains: Standardized path coefficients, posterior probabilities, model IDs
naturalization_data <- read_csv(here("results", "Table_S2_hypothesis_tests.csv"))

# 4. Data Processing -----------------------------------------------------------
# Hypothesis: Cultivation history (China vs foreign) and native range extent
#             differentially affect naturalization success across life forms

# Extract fixed effects (exclude random effects and intercepts)
# Rationale: Focus on predictor-response relationships, not baseline levels
naturalization_data <- naturalization_data %>%
  filter(effect == "fixed" & term != "(Intercept)") %>%
  mutate(term = str_replace_all(term, "[.]", "")) %>%  # Remove dots from variable names
  mutate(
    estimate  = round(estimate, 3),                     # Standardized coefficients
    Post.Prob = round(Post.Prob, 4)                     # Posterior probability (Bayesian p-value)
  ) %>%
  # Assign significance stars based on Bayesian posterior probabilities
  # Thresholds: 95% (*), 99% (**), 99.9% (***)
  mutate(Star = case_when(
    Post.Prob > 0.999 ~ "***",
    Post.Prob > 0.99  ~ "**",
    Post.Prob > 0.95  ~ "*",
    Post.Prob < 0.95  ~ ""
  )) %>%
  # Color-code edges by effect direction and significance
  # Rationale: Red = facilitation, Blue = inhibition, Grey = uncertain
  mutate(color = case_when(
    estimate > 0 & Post.Prob > 0.95 ~ "red",
    estimate < 0 & Post.Prob > 0.95 ~ "blue",
    Post.Prob < 0.95                ~ "grey"
  ))

# Prepare edge list for network construction
# Edge width scaled by absolute effect size for visual emphasis
network_edges <- naturalization_data %>%
  transmute(
    from      = term,                              # Predictor variable
    to        = response,                          # Response variable
    model     = Model,                             # Model stratification (life form)
    estimate  = estimate,                          # Path coefficient
    Post.Prob = Post.Prob,                         # Bayesian confidence
    width     = abs(estimate) * 5,                 # Visual scaling factor
    color     = color,                             # Edge color (effect direction)
    label     = paste(estimate, Star)              # Label: "0.42 **"
  )

# 5. Visualization -------------------------------------------------------------
# Generate separate SEM diagrams for each life form group
# Rationale: Different life history strategies may respond differently to
#            propagule pressure and native range predictors

# Model names MUST match the "Model" column in your data exactly
model_names    <- c("All species", "Annual herbs", 
                    "Perennial herbs", "Woody plants")
response_vars  <- c("All species", "Annual herbs", 
                    "Perennial herbs", "Woody plants")

sem_plots_naturalization <- purrr::map2(
  model_names, 
  response_vars, 
  ~ visualize_sem_network(network_edges, .x, .y)
)
sem_plots_naturalization <- setNames(sem_plots_naturalization, model_names)

# Extract individual plots for flexible composition
sem_plot_all_species    <- sem_plots_naturalization$`All species`
sem_plot_annual_herb    <- sem_plots_naturalization$`Annual herbs`
sem_plot_perennial_herb <- sem_plots_naturalization$`Perennial herbs`
sem_plot_woody          <- sem_plots_naturalization$`Woody plants`

# 6. Multi-Panel Figure Assembly -----------------------------------------------
# Compose Figure 2: Life form-specific SEM pathways
# Layout: Three panels (excluding all-species model for main text)
figure_2_sem_life_forms <-
  sem_plot_annual_herb +
  sem_plot_perennial_herb +
  sem_plot_woody +
  plot_layout(ncol = 3) +
  plot_annotation(tag_levels = "A")

print(figure_2_sem_life_forms)

# Export publication-quality figure (300 dpi, PNG format)
ggexport(
  figure_2_sem_life_forms, 
  filename  = here("results", "Figure.2.png"),
  width     = 8000,
  height    = 3000,
  pointsize = 12,
  res       = 300
)

# Note: Figure 2 was manually enhanced in draw.io for final presentation

# ==============================================================================
# END OF SCRIPT
# ==============================================================================