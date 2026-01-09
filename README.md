# Paper I - Objective carcass grading for bovine animals based on carcass length

### Statistical Validation of Objective Classification Models

### 📄 Overview
This repository contains the **R analysis scripts** used for *Paper I: Objective carcass grading for bovine animals based on carcass length*.

The project involved the statistical validation of a new objective prediction model intended to replace subjective human classification in the Norwegian market. The core objective was to quantify the precision and bias of the prediction models relative to the existing human classifier standard.

### 🛠️ Statistical Methodology
The analysis utilizes **Variance Component Analysis (Random Effects Models & ICC)** to decompose the total error into its constituent sources. By modeling the classification process hierarchically, we isolated the systematic signal from random noise.

**Key Statistical Concepts Implemented:**
*   **Bias Estimation ($\mu$):** Quantifying the systematic deviation that resulted from the prediction models.
*   **Variance Decomposition:**
    *   **Prediction Model Error ($\sigma^2_t$):** The shared variance component attributable to the model's inability to perfectly predict conformation.
    *   **Classifier Error ($\sigma^2_u$):** The independent random noise introduced by individual human classifiers.
*   **Intraclass Correlation (ICC):** Calculated to benchmark the reliability of the objective method against human classifiers.

### 📂 Repository Contents
The analysis is split into execution code and helper functions to ensure modularity.

*   **`OBKLAS_Rcode.R`**: The main execution script. Performs data import, model fitting (`lme4`), residual analysis, and parameter estimation.
*   **`OBKLAS_Rfunctions.R`**: Custom R functions developed for the project used in **`OBKLAS_Rcode.R`**.
*   **`LICENSE`**: CC0 1.0. (see repository link under "References")

### 💻 Technologies
*   **Language:** R
*   **Key Packages:** `lme4` (Linear Mixed Models), `tidyverse` (Data Manipulation).

---
### 🔗 References
*   **Full Article:** [Taylor & Francis Online](https://www.tandfonline.com/doi/full/10.1080/09064702.2021.1906940)
*   **Data Repository:** [DataverseNO](https://dataverse.no/dataset.xhtml?persistentId=doi:10.18710/T9SXVF)
