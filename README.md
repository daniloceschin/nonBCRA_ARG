# Non-BRCA Pathogenic Variants in Córdoba, Argentina

## Study Overview

This repository contains the statistical analysis code and supplementary materials for the research paper:

**"PALB2 Founder Variant and High Non-BRCA Pathogenic Variant Prevalence with Clinical Implications in Córdoba, Argentina"**

### Authors
- Claudia Alejandra Martin¹,²
- Catalina Bono¹  
- Juliana Mandrile¹
- Christine Susana Mercedes Kunst¹
- Adriana Elizabeth Borello¹
- Rodolfo Ávila¹
- Florencia Celeste Salazar¹
- Virginia Miretti¹,²
- Florencia Pabletich¹,²
- Verónica Andreoli¹*
- Danilo Guillermo Ceschin¹,²,³*

*Co-senior authors

### Affiliations
1. Hospital Privado Universitario de Córdoba, Córdoba, Argentina
2. Instituto Universitario de Ciencias Biomédicas de Córdoba (IUCBC), Córdoba, Argentina  
3. Centro de Investigación en Medicina Traslacional "Severo R. Amuchástegui" (CIMETSA), CONICET, Córdoba, Argentina

## Study Description

This retrospective study analyzed 283 patients with personal and/or family history of breast or ovarian cancer evaluated at Hospital Privado Universitario de Córdoba between 2016-2022. The study identified a high prevalence of pathogenic variants in non-BRCA genes (11.3%), with PALB2 being the most frequent (4.59%), including evidence of a potential founder effect for the c.1653T>A variant.

### Key Findings
- **Total pathogenic variants identified:** 85/283 patients (30.0%, 95% CI: 24.7-35.8%)
- **Non-BRCA pathogenic variants:** 32/283 patients (11.3%, 95% CI: 8.0-15.5%)
- **PALB2 prevalence:** 13/283 patients (4.59%, 95% CI: 2.6-7.7%)
- **PALB2 c.1653T>A founder variant:** 5/283 patients (1.77%, 95% CI: 0.7-4.1%)
- **Statistical power:** >99% for key comparisons with European populations

## Repository Contents

### Main Analysis Script
- **`NonBCRA_statistics.R`** - Complete statistical analysis script including:
  - Confidence interval calculations using exact binomial methods
  - Statistical power analysis for population comparisons
  - Fisher's exact tests comparing Córdoba cohort with European populations
  - Summary statistics and data export functions

### Generated Output Files
- **`confidence_intervals.csv`** - Frequency estimates with 95% confidence intervals
- **`population_comparisons.csv`** - Statistical comparisons with European cohorts
- **`power_analysis.csv`** - Statistical power calculations for key findings

## Requirements

### R Version
- R version 4.4.1 or higher

### Required Packages
```r
install.packages(c("binom", "stats"))
```

## Usage Instructions

1. **Clone the repository:**
```bash
git clone https://github.com/daniloceschin/nonBCRA_ARG.git
cd nonBCRA_ARG
```

2. **Open R/RStudio and set working directory:**
```r
setwd("path/to/nonBCRA_ARG")
```

3. **Install required packages:**
```r
install.packages(c("binom", "stats"))
```

4. **Run the analysis:**
```r
source("NonBCRA_statistics.R")
```

## Key Statistical Methods

### Confidence Intervals
- **Method:** Exact binomial confidence intervals (Clopper-Pearson method)
- **Confidence Level:** 95%
- **Implementation:** `binom.confint()` function with `methods="exact"`

### Statistical Power Analysis
- **Test Type:** Two-proportions z-test
- **Significance Level:** α = 0.05 (two-tailed)
- **Power Calculation:** Normal approximation method

### Population Comparisons  
- **Statistical Test:** Fisher's exact test
- **Comparison Groups:** Córdoba cohort (n=283) vs European literature (n≈1000)
- **Significance Level:** α = 0.05

## Results Summary

| Gene | Córdoba Frequency | European Frequency | Statistical Significance |
|------|-------------------|-------------------|-------------------------|
| PALB2 | 4.59% (2.6-7.7%) | ~1.0% | p<0.001 |
| CHEK2 | 2.47% (1.2-5.0%) | ~2.5% | p>0.05 |
| TP53 | 1.06% (0.4-3.0%) | ~0.5% | p>0.05 |
| Non-BRCA Total | 11.3% (8.0-15.5%) | ~6.5% | p<0.05 |

## Clinical Implications

The study provides evidence-based recommendations for:
- PALB2 carrier management (surveillance protocols, risk-reducing surgery)
- CHEK2 and TP53 carrier management  
- Population-specific genetic testing strategies for Argentina
- Healthcare system implementation considerations
- Cascade testing protocols for founder variants

## Data Availability

The raw data supporting the conclusions of this article are available from the corresponding authors upon reasonable request, subject to appropriate ethical approvals and data sharing agreements.

## Ethics and Compliance

This study was conducted in accordance with the Declaration of Helsinki and approved by the local institutional ethics committee. All participants provided informed consent for genetic testing and research use of de-identified data.

## Contact Information

**Corresponding Authors:**
- Dr. Danilo Guillermo Ceschin: danilo.ceschin@iucbc.edu.ar
- Dr. Verónica Andreoli: veronica.andreoli@hospitalprivadosa.com.ar

## Citation

If you use this code or data in your research, please cite:

......

## License

This work is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

We thank all patients and families who participated in this study, and the healthcare professionals involved in patient care and sample collection at Hospital Privado Universitario de Córdoba.

---

**Repository Last Updated:** January 2025  
**Manuscript Status:** Submitted to Journal of Medical Genetics  
**Data Collection Period:** January 2016 - December 2022
