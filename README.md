# Charter Schools, Learning Rates, and the Geography of Achievement: Multilevel Models on SEDA

A multilevel analysis of U.S. public school achievement using the Stanford Education Data Archive (SEDA 6.0). The project has two parts:

1. **Primary (causal):** Do charter schools show higher learning rates than comparable traditional public schools in the same district?
2. **Secondary (descriptive):** How is the variance in school achievement split across states, districts, and schools?

> **Status:** In progress. The data have been obtained and profiled; modeling has not yet started. This README will be updated as the project develops.

## Data

This project uses the SEDA (Stanford Education Data Archive) dataset, version 6.0. SEDA contains average scores on mathematics and ELA (English Language Arts) state tests for grades 3–8 across most public schools in the country between 2009 and 2019. It places all states on a common scale (linked to NAEP), which allows direct comparisons, and it reports a standard error with each estimate.

The dataset offers two outcomes of interest: school average achievement and school learning rate, i.e., growth across grades. Predictors are available at more than one level. At the school level, there are predictors like the percent of economically disadvantaged students, demographics, enrollment, and locale. At the district level, there are predictors such as an SES composite and the percent of English learners.

### Files

Four files are used, each pooled across years, grades, and subjects:

| File | Contents | Role |
|---|---|---|
| `seda_school_pool_gcs_6.0.csv` | One row per school (64,415 schools): average achievement (`gcs_mn_avg_ol`) and learning rate (`gcs_mn_lrn_ol`), with standard errors, in grade-level units (4 = the national average NAEP 4th-grade score; one unit ≈ one grade level) | Outcomes |
| `seda_cov_school_pool_6.0.csv` | School-level predictors (demographics, `perecd`, `lep`, enrollment, locale, 2018 charter status) and each school's geographic district ID (`sedalea`) | School predictors; defines the school-in-district nesting |
| `seda_cov_geodist_pool_6.0.csv` | District-level predictors (ACS SES composite `sesavgall`, `perell`, segregation indices) | District predictors, linked on `sedalea` |
| `seda_crosswalk_6.0.csv` | School-by-year records (2009–2019), including charter status in every year and school latitude/longitude | Needed only for the charter analysis |

Merge keys: outcomes ↔ school covariates on `sedasch`; school ↔ district covariates on `sedalea`. Variable definitions are in the SEDA codebooks (`codebooks/`).

### Getting the data

SEDA is free to download from the [SEDA archive](https://edopportunity.org/get-the-data/seda-archive-downloads) after accepting its data use agreement. Place the four CSVs above in `data/`.

Citation: Reardon, S. F., Ho, A. D., Shear, B. R., Fahle, E. M., Kalogrides, D., Saliba, J., et al. *Stanford Education Data Archive (Version 6.0).* https://purl.stanford.edu/xh833nn4025

## Design

The data are cross-sectional with three levels: about 64 thousand schools nested in about 12.4 thousand geographic districts nested in 50 states plus Washington DC. A further 705 schools sit in 69 non-geographic, state-level placeholder districts that have no district covariates.

**Primary research question (causal).** Do charter schools show higher learning rates than comparable traditional public schools in the same district? The outcome is learning rate, and the treatment is attending a charter school; treatment assignment (charter status) is in turn the outcome of the propensity score model. Relevant predictors may include the share of economically disadvantaged students (`perecd`), the percent of limited English proficient students (`lep`), race, enrollment, and the district SES composite.

**Secondary research question (variance decomposition).** How is the variance in school achievement split across states, districts, and schools, with a possible extension to a geospatial component? This uses a more standard `lmer` approach with state and district information. Unlike the primary question, the outcome is average achievement, not learning rate.

## Analysis plan

**Primary question.** An MLM on top of propensity score modeling. The propensity score is itself estimated with a multilevel logistic model, for example:

```r
glmer(charter ~ ecd + lep + race + enrollment + ses_dist + (1 | district),
      family = binomial)
```

The outcome model is then an MLM of learning rate, fit within propensity score strata or matched sets, e.g.:

```r
lmer(lrn ~ charter + covariates + (1 | state/district))
```

**Secondary question.** The model is built as an extension of:

```r
lmer(avg ~ 1 + (1 | state/district))
```

**Diagnostics.**

- *Propensity score:* overlap plots of propensity scores; covariate balance after matching/subclassification; a sensitivity analysis for unmeasured confounding.
- *Multilevel models (both questions):* likelihood-ratio tests on the random effects (conservative, since variances are tested on the boundary); normal QQ plots of district- and state-level random effects.

## Anticipated data issues

- **Missing learning rates.** Missing data is minimal for average achievement and for most covariates (under 2% of schools), but not for learning rate. It is missing for 14.5% of schools, mostly those with few tested grades, and for 31% of charter schools compared with 13% of traditional public schools. This gap matters for the primary question and needs to be accounted for.
- **Charter status over time.** Charter status is recorded for 2018 in the school covariates file and year-by-year in the crosswalk. 423 schools in the sample changed charter status at some point between 2009 and 2019.
- **Co-location.** Some independently operated charter schools share a building with district public schools. This creates potential spillover: if charter status is the treatment, co-location would violate the no-interference assumption. Co-location, especially on vulnerable campuses, has been known to create pedagogical and practical problems for the public school involved. The plan is to find co-location arrangements and handle them uniformly, for instance by dropping them. The crosswalk's coordinates make a first pass possible: about 580 schools in the outcome file share exact coordinates with a school of the other type in the same year, mostly in CA, NY, and WI. Shared coordinates are only a proxy, though.
- **Finding comparable schools within districts.** Comparable traditional public schools may be hard to find for some charter schools within the same district.
- **Narrowing to districts with both school types.** Only about 10% of districts contain any charter school. After dropping schools with missing learning rates, about 930 districts with both school types remain, with roughly 2,900 charters and 18,700 traditional public schools. If matching within districts proves not feasible, one option is matching across districts under both geospatial and modeling-adjacent constraints (e.g., using estimated BLUPs as a similarity metric).

## Repository structure

```
codebooks/   SEDA 6.0 codebooks (.xlsx)
data/        SEDA CSVs (see "Getting the data")
```
