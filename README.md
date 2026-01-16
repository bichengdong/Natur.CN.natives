## README: Data and R Code for "Cultivation Interacts with Life Form to Mediate the Effect of Native Range Size on Plant Naturalization"

---

### Data

**1. 20250330.native_plant.matched.cultivated_plant_DO11.csv**

**Description**: This dataset constitutes the core data for analyzing the global naturalization patterns of native Chinese plants. It integrates species naturalization information, cultivation information, life forms, native range size, and cultivation both within and outside China.

*Taxonomic identifiers*

- `accepted_plant_name_id`: The unique identifier for the accepted name by WCVP
- `taxon_name`: The accepted name (Latin name) of the species
- `taxon_authors`: The authors who named the species

*Life forms, cultivation status, and naturalization status*

- `life.form.integrated`: Integrated life form of the species
- `cultivation.lin2019`: Cultivation status determined according to *Catalogue of Cultivated Plants in China* (yes or no)
- `nat.extent`: Naturalization extent (number of regions where plants are naturalized)
- `nat.incidence`: Naturalization incidence (likelihood of becoming naturalized)

*Native range size and cultivation extent*

- `native.global.tdwg3`: Native range size (that is, number of regions where native)
- `planting.China.tdwg3`: Cultivation within China, measured as the number of TDWG level-3 regions within China where the species is recorded as cultivated
- `WorldCuP.n.tdwg3`: Cultivation outside China, measured as the number of TDWG level-3 regions outside China where the species is recorded as cultivated

**2. NativeChinesePlantsCultivationExtentWorldCuP.csv**

**Description**: This dataset quantifies the cultivation extent of Chinese native plants outside of China. Out of the 25,808 native taxa analyzed, 9,670 are identified as being cultivated internationally.

- `taxon_name`: The accepted name (Latin name) of the species
- `taxon_authors`: The authors who named the species
- `accepted_plant_name_id`: The unique identifier for the accepted name
- `number_records`: The number of unique WorldCuP lists in which each taxon is present
- `number_tdwg3`: The total number of unique TDWG level-3 regions outside of China where the taxon is recorded as cultivated according to WorldCuP

**3. References_WorldCuPSources.csv**

**Description**: This file contains the comprehensive list of 365 bibliographic references and data sources used to compile the cultivation extent information for Chinese native plants within the World Checklist of Cultivated Plants (WorldCuP) database.

---

### Code

#### 1. Rcode for Figure_S4 and Table_S2.R

**Purpose**: Bayesian structural equation modeling (BSEM) to quantify direct and indirect effects of native range size on naturalization success via cultivation pathways.

**Key analyses**:

- Specify three-step mediation model:
  1. Native range size → Cultivation within China
  2. Native range size → Cultivation outside China (direct + mediated via China)
  3. Native range size → Naturalization extent (direct + dual-mediated pathways)
- Fit separate BSEMs for annual herbs, perennial herbs, woody plants, and all species combined
- Extract posterior distributions of pathway-specific effects (direct, indirect via China, indirect via foreign cultivation, indirect via both)
- Conduct Bayesian hypothesis testing for directional effects (β > 0 or β < 0)

**Outputs**:

- Figure S4: Ridge density plots of posterior effect distributions stratified by life form
- Table S2: Hypothesis test results with Bayesian credible intervals and evidence ratios

---

#### 2. Rcode for Figure_2.R

**Purpose**: Network visualization of structural equation model (SEM) pathways linking cultivation history, native range size, and naturalization extent across plant life forms.

**Key analyses**:

- Construct directed acyclic graphs (DAGs) representing causal pathways from Bayesian SEM results
- Visualize standardized path coefficients with effect direction (facilitation vs. inhibition) and statistical significance
- Generate life form-specific network diagrams (annual herbs, perennial herbs, woody plants)

**Output**: Figure 2 (3-panel SEM pathway comparison across life forms)

---

#### 3. Rcode for Figures_1, S2 and Table_S1.R

**Purpose**: Comprehensive hurdle model analysis examining how native range size, cultivation extent (within and outside China), and life form predict plant naturalization incidence and extent.

**Key analyses**:

- Fit full hurdle model with life form interactions (combines binomial and zero-truncated negative binomial components)
- Fit component models separately for detailed inspection:
  - Binomial model: probability of naturalization (yes/no)
  - Zero-truncated negative binomial: naturalization extent (number of TDWG3 regions, conditional on occurrence)
- Generate model predictions across focal variable gradients
- Create visualization panels with inset zoom plots for low-range details

**Outputs**:

- Figure 1: Hurdle model predictions (naturalization extent vs. focal predictors)
- Figure S2: Component model predictions (binomial incidence + truncated count model)
- Table S1: Model coefficients with standard errors and p-values
- Additional output: 9-panel combined figure (3 predictors × 3 models)

---

#### 4. Rcode for Figure_S3.R

**Purpose**: Visualization of linearized relationships from two-part hurdle models examining predictors of naturalization extent.

**Key analyses**:

- Fit negative binomial hurdle models with life form interactions
- Zero component: model probability of naturalization occurrence (binomial regression)
- Count component: model naturalization breadth conditional on occurrence (zero-truncated negative binomial)
- Generate linear predictor plots across focal variables (native range size, cultivation within China, and cultivation outside China)

**Output**: Figure S3 (6-panel plot: 3 predictors × 2 model components)

---

#### 5. Rcode for Figure_S6.R

**Purpose**: Pairwise correlation analysis of native range size, cultivation within China, and cultivation outside China.

**Key analyses**:

- Compute Pearson correlations among log-transformed, z-standardized predictors
- Visualize correlation matrix with scatterplots (lower triangle), correlation coefficients (upper triangle), and marginal distributions (diagonal)
- Test hypothesis that native range size positively correlates with both cultivation within China and cultivation outside China

**Output**: Figure S6 (correlation matrix with regression lines and distribution histograms)
