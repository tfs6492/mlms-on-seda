# Charter and Traditional Public Schools in a Multilevel Frame

### Variance components, response surfaces, and measurement error in SEDA 6.0

Author: Thomas Szymanski Status: 🛠️ In progress

------------------------------------------------------------------------

## Overview

This project uses the Stanford Education Data Archive (SEDA 6.0) to study how school **learning rates** (and, as a companion outcome, **average achievement**) vary across states, districts, and schools in the United States. It also asks how that variation **differs between charter schools and traditional public schools (TPS)**.

The analysis is descriptive. Sector is treated as a grouping that may reshape the entire structure of the multilevel model:

- **how much variance** sits at the state, district, and school levels;
- **the response surface**, meaning the functional form linking the outcome to school composition and district context;
- **how precisely** each school's outcome is measured, since charter schools tend to be smaller and their estimates noisier.

A separate, secondary analysis models the **ecological selection of schools into the charter sector**. "Ecological" means it uses school- and district-level aggregates, not student-level data.

The emphasis throughout is on choosing, fitting, and defending an appropriate multilevel model: specification, heteroscedastic variance components, known sampling error, missing data, and diagnostics. All results are interpreted as **conditional associations**, not causal effects.

------------------------------------------------------------------------

## Data

**Source:** Reardon, S. F., Ho, A. D., Shear, B. R., Fahle, E. M., Kalogrides, D., Saliba, J., et al. *Stanford Education Data Archive (Version 6.0).* <https://purl.stanford.edu/xh833nn4025>

- U.S. public schools, grades 3–8, math and reading/ELA state tests, school years 2008–09 through 2018–19.
- All states linked to a **common NAEP-based scale**, so schools are comparable across states.
- Metric: **GCS (grade cohort scale)**. 4 equals the national average NAEP 4th-grade score, and one unit is roughly one grade level. A learning rate of 1 means one grade level of growth per grade.
- Every estimate comes with a **standard error**.

### Files (CSV, "pool" form: pooled across years, grades, and subjects)

| File | Unit | Role |
|----|----|----|
| `seda_school_pool_gcs_6.0` | school | Outcomes and their standard errors |
| `seda_cov_school_pool_6.0` | school | School covariates; contains `sedalea` (geographic district) and `stateabb`, which define the nesting |
| `seda_cov_geodist_pool_6.0` | geographic district | District covariates, linked on `sedalea` |
| `seda_crosswalk_6.0` | school × year | Year-by-year charter status, virtual-school flag, latitude/longitude |

### Key variables

| Variable | File | Description |
|----|----|----|
| `gcs_mn_lrn_ol`, `gcs_mn_lrn_ol_se` | outcomes | **Learning rate** (growth per grade) and its SE (primary outcome) |
| `gcs_mn_avg_ol`, `gcs_mn_avg_ol_se` | outcomes | **Average achievement** and its SE (companion outcome) |
| `charter` | school covariates / crosswalk | Charter status (CCD 2018; yearly in crosswalk) |
| `perecd` | school covariates | Share of economically disadvantaged students |
| `lep` | school covariates | Share of limited-English-proficient students (CRDC average) |
| `perwht`, `perblk`, `perhsp`, `perasn`, … | school covariates | Race/ethnicity shares |
| `totenrl`, `mingrd`, `maxgrd`, `level` | school covariates | Enrollment and grade span |
| `sesavgall` | district covariates | ACS socioeconomic composite |
| `perell`, `urbanicity`, `hsecdnec` | district covariates | \% English learners, locale, economic segregation |

### Data access

SEDA is free to download after accepting a data use agreement, which does not allow redistribution. **This repository contains code only.** To reproduce the analysis, download the files above from the [SEDA archive](https://edopportunity.org/get-the-data/seda-archive-downloads) and place them in `data/raw/`.

------------------------------------------------------------------------

## Design

Cross-sectional, three-level nesting:

```         
state (50 + DC)
 └── geographic district (~12,400; key: sedalea)
      └── school (~64,000; key: sedasch)   ← sector: charter vs. TPS
```

Structural features that shape the modeling:

- **Sector varies within districts.** About 10% of districts contain charter schools. After dropping schools with missing learning rates, about 930 districts contain both sectors (roughly 2,900 charters and 18,700 TPS). Only these districts directly inform the *within-district* contrast and the district-level covariance between sectors.
- **"District" means something different for each sector.** SEDA 6.0 places charter schools in the geographic district where they physically sit. For a charter, the district effect therefore reflects local context (neighborhood, labor market, peer schools), not governance by that district.
- **Many single-school districts**, mostly rural. They identify little about district-level variance.
- **States with few or no charter schools** contribute TPS information only, which limits how well state-level charter variance can be estimated.
- **Non-geographic placeholder districts** (705 schools in 69 state-level placeholders) have no district covariates.

------------------------------------------------------------------------

## Research questions

### RQ1 (primary): How do the variance components and the response surface differ between charter and traditional public schools, controlling for district and state?

- **1a. Variance partition by sector.** How is the variance in school learning rates split across states, districts, and schools, and is that split different for charters and TPS? Are charter schools more heterogeneous within a district than TPS, and are a district's charter and TPS outcomes correlated?
- **1b. Response surface by sector.** Does the functional form of the outcome given school composition (poverty, English learners, race/ethnicity, enrollment, grade span) and district context (SES) differ between the two sectors? This includes slopes, nonlinearity, and within- vs. between-district relationships.
- **1c. Signal vs. noise.** How much of the apparent between-school variance, and of any sector difference in that variance, is sampling error rather than true variation?

### RQ2 (secondary, modeled separately): Which school and district characteristics are associated with a school being a charter?

This is an **ecological** analysis. School composition is partly a *consequence* of charter status (who enrolls), so the model describes how the sectors differ and where charters are located. It does not explain why families choose them.

------------------------------------------------------------------------

## Data issues

- **Missing learning rates.** Learning rates are missing far more often than average achievement, and much more often for charter schools than for TPS. Some are structurally undefined (too few tested grades to estimate growth); others are suppressed (small cells or imprecision).
- **Heterogeneous sampling error.** Standard errors vary widely across schools. Smaller schools, and schools testing fewer grades, have noisier estimates.
- **Raw vs. empirical-Bayes estimates.** SEDA provides both raw (`_ol`) and already-shrunken empirical-Bayes (`_eb`) estimates. The raw estimates are used to avoid shrinking twice inside a multilevel model.
- **Changes in sector over time.** Some schools converted to or from charter status during 2009–2019, so a single sector label is ambiguous for them.
- **Virtual schools.** Virtual schools have no meaningful physical location and are flagged separately.
- **Co-located schools.** Some charter schools share a building with a TPS, which shows up as identical coordinates in the crosswalk.
- **Non-geographic placeholder districts.** A small number of schools sit in state-level placeholder districts with no district covariates.
- **Single-school districts.** Many rural districts contain only one school.
- **Special-education schools.** Schools serving mainly special-education populations are not comparable to general schools.
- **Free/reduced-price lunch measures.** The Community Eligibility Provision distorts these, so the economically-disadvantaged share (`perecd`) is used instead.
