# Robustness to fat-tailed school residuals: a Student-t (t4) likelihood in `brms` vs. weighting in `lme4::lmer()`

*Theory and feasibility notes for the robustness section of the SEDA multilevel analysis. No code. The goal is to understand what each method does, why it would help, where it would not, and what to compare at the end.*

------------------------------------------------------------------------

## 0. Where we stand: what the diagnostics actually say

Everything here concerns the current reference models in `mlms_on_seda.Rmd`:

- `fit.full.rndm`: sector, Mundlak term, school composition, district context, a state random intercept, and correlated sector-specific district random effects `(0 + factor(charter) | sedalea)`.
- `fit.full.rndm.int`: the same model plus sector × covariate interactions (the RQ1b response-surface test).

Both are fitted on `seda_an`: N = 54,074 schools (2,998 charters), 11,569 districts, 51 states.

I refit both models from `data/seda_an` to get the numbers below. They match the figures in `figs/`.

### 0.1 Level 1 (school residuals): clearly heavy-tailed

`figs/resid_l1_school.png` and `figs/resid_int_l1_school.png` look the same: a symmetric, sharply peaked center and long tails on both sides. The lower tail is longer.

| Quantity (interaction model) | Value |
|----|----|
| Residual SD σ̂ | 0.203 GCS |
| Most extreme scaled residuals | −15.1 σ and +8.0 σ |
| Excess kurtosis of residuals (0 if normal) | ≈ 4.2 |
| Schools with \|residual\| \> 3σ | 406 (a normal model expects ≈ 146) |
| Schools with \|residual\| \> 4σ | 98 (normal expects ≈ 3) |
| Schools with \|residual\| \> 6σ | 19 (normal expects ≈ 0) |
| Between ±2σ | Close to normal (2,392 vs. ≈ 2,460 expected beyond ±2σ) |

So the bulk of the distribution is close to normal and the problem is concentrated in the tails. That pattern is exactly what a Student-t likelihood is built to describe.

### 0.2 Levels 2 and 3

- **TPS district BLUPs** (`figs/resid_l2_tps.png`): also fat-tailed, but milder than at level 1.
- **Charter district BLUPs** (`resid_l2_charter.png`): close to normal, with a slight upper tail.
- **State BLUPs** (`resid_l3_state.png`): fine.

This matters later: a t-likelihood at level 1 does **not** address the TPS district tails by itself (§2.6).

### 0.3 Who are the outliers?

- **Not a charter phenomenon.** Of the 98 schools beyond ±4σ, 93 are TPS and 5 are charters (5.1%, vs. a 5.5% charter base rate). The *within-sector* residual distributions differ, though: TPS residuals have excess kurtosis ≈ 4.5, while charter residuals have ≈ 1.0. Robust methods will therefore mostly change how TPS schools inform the model.
- **Noisier, but not explained by noise.** Extreme schools have larger standard errors (mean SE 0.091 vs. 0.062 overall; see `hist_se_outliers.png` vs. `hist_se_all.png`), and the correlation between \|residual\| and SE is 0.23. But the largest possible SE in the data is ≈ 0.15, which is under 1σ. Sampling error alone cannot produce a residual of −3 GCS (a school "losing" two grade levels per grade). Averaged over schools, sampling variance (mean SE² ≈ 0.0045) is only **≈ 11% of the residual variance σ̂² ≈ 0.041**. Most school-level residual variation is real between-school variation. The extreme tail is probably a mix of genuinely unusual schools and data or measurement problems such as few tested grades, restructuring, sector switches, or suppressed cells.
- `se_vs_resid.png` shows this directly. The outliers are not lined up at high SE. They occur across the SE range.

**Implication for both methods.** "Downweight imprecise schools" (inverse-variance weighting) and "downweight surprising schools" (robust likelihood or robust weights) are **two different operations aimed at two different problems**. The data have both problems, but the fat tails are mainly the second. Most of what follows rests on keeping these two apart.

------------------------------------------------------------------------

## 1. Foundations of Bayesian inference

### 1.1 The core idea

In the frequentist (REML/ML) framework used so far, parameters such as the charter gap β or the variance components are **fixed unknown constants**. We find the values that make the observed data most likely. Uncertainty is described by how estimates would vary over hypothetical repeated samples (standard errors, confidence intervals, p-values).

In the Bayesian framework, unknown parameters are treated as **random quantities whose uncertainty we describe with probability distributions**. Inference is the process of updating those distributions with data, using Bayes' theorem:

$$\underbrace{p(\theta \mid y)}_{\text{posterior}} \;\propto\; \underbrace{p(y \mid \theta)}_{\text{likelihood}} \times \underbrace{p(\theta)}_{\text{prior}}$$

- **Likelihood** $p(y \mid \theta)$: the data model. It answers "if the parameters were θ, how probable are the data we saw?" This is **the same object REML/ML maximizes**. The choice of normal vs. Student-t errors lives here. That is why the t4 idea is mainly a change to the *likelihood*, not something exclusively Bayesian.
- **Prior** $p(\theta)$: what we assume about plausible parameter values before seeing these data.
- **Posterior** $p(\theta \mid y)$: the updated distribution. Everything we report comes from it.

### 1.2 How to read Bayesian output

- **Point estimate.** Usually the posterior mean or median, analogous to the REML estimate.
- **Uncertainty.** A **95% credible interval** is the range containing 95% of the posterior. Unlike a confidence interval, it can be read directly as "given the model and the data, there is a 95% probability the parameter lies here".
- **No p-values by default.** Instead of asking whether β ≠ 0, you look at the posterior directly, for example the probability that the charter gap is positive, or that a sector × covariate interaction exceeds a substantively meaningful size.

### 1.3 Why priors matter (and when they don't)

With 54,074 schools, the likelihood for the fixed effects and the level-1 scale is very informative. Any reasonable, weakly informative prior will be overwhelmed and the posterior will be driven by the data. Priors matter more where information is thin:

- **State variance**: only 51 states.
- **Charter district variance and the TPS–charter correlation**: informed only by the ≈ 930–990 districts that contain charters.
- **The degrees-of-freedom parameter ν**, if estimated (§2.4).

"Weakly informative" means priors that rule out absurd values (for example, a district SD of 3 GCS when the outcome SD is 0.23) without pushing estimates toward any particular answer. For a robustness check, the principle is: **choose priors that are defensible on the GCS scale, report them, and confirm results don't hinge on them.**

### 1.4 How the posterior is computed: MCMC

For multilevel models the posterior has no closed form. `brms` writes the model in **Stan**, which draws thousands of samples from the posterior using **Hamiltonian Monte Carlo (HMC)**. Every quantity you report (means, intervals, variance ratios, ICCs) is computed from those draws. Derived quantities such as an ICC or a pseudo-R² get full uncertainty intervals for free, which is awkward to obtain from `lmer`.

Because it is a simulation, you have to check that it worked:

| Diagnostic | What it checks | Rule of thumb |
|----|----|----|
| **R-hat** | Whether the independent chains agree | ≤ 1.01 |
| **Effective sample size (ESS)** | How many independent draws the correlated samples are worth | Bulk and tail ESS in the hundreds or more |
| **Divergent transitions** | Stan struggling with the geometry of the posterior, common with small variance components | Ideally zero |
| **Posterior predictive checks** | Whether data simulated from the fitted model look like the real data | The Bayesian analog of the Q-Q plot |

The **posterior predictive check** is the key one here. You simulate outcome distributions from the Gaussian model and from the t4 model and overlay them on the observed outcome. The Gaussian model should visibly fail to produce schools as extreme as the observed ones; the t4 model should not.

### 1.5 `brms` vs. `lmer`: what transfers

`brms` uses the same formula syntax as `lme4`. The current model structure, including the correlated sector-specific district effects and the state intercept, carries over one-to-one. What is new is that you must specify (a) a **family** (Gaussian or Student-t), (b) **priors**, and (c) **sampler settings**.

------------------------------------------------------------------------

## 2. The Student-t (t4) error structure

### 2.1 What changes in the model

Currently, at level 1:

$$y_{ijk} = \mathbf{x}_{ijk}'\boldsymbol\beta + u_{jk} + v_k + \varepsilon_{ijk}, \qquad \varepsilon_{ijk} \sim \mathcal N(0, \sigma^2)$$

The t4 model changes only the last part:

$$\varepsilon_{ijk} \sim t_{\nu}(0, \sigma), \qquad \nu = 4$$

The fixed effects, random effects, and their covariance structure are untouched. The Student-t distribution looks like a normal distribution but puts far more probability in the tails. ν (the "degrees of freedom") controls how heavy those tails are: ν = 1 is the extremely heavy-tailed Cauchy, and ν → ∞ returns the normal. ν = 4 is a conventional "robust but well-behaved" choice: the variance is still finite, but the tails are much heavier than normal.

### 2.2 Why it limits the influence of outliers

**Under normal errors.** The penalty for a residual grows with its *square*. A school 15σ from its prediction contributes 225 times the "pull" of a school 1σ away. The estimator moves β, σ, and the BLUPs substantially to reduce that one enormous discrepancy. Specifically:

- **σ̂ inflates** to accommodate the extremes. That makes every school look relatively more precise than it is and distorts the variance partition (RQ1a).
- **Fixed effects can be pulled** toward the outliers, especially when outliers sit at unusual covariate values (high leverage).
- **District BLUPs absorb outliers.** One wild school in a 3-school district drags that district's effect.

**Under t errors.** The penalty grows only *logarithmically* for large residuals. The model effectively says "an observation this far out is more plausibly an unusual case than a reason to move the whole fit".

A useful way to see this is the **scale-mixture representation**. A t-distributed error is equivalent to a normal error whose variance is individually scaled for each school:

$$\varepsilon_i \mid \lambda_i \sim \mathcal N(0, \lambda_i\sigma^2), \qquad \lambda_i \sim \text{Inverse-Gamma}(\nu/2, \nu/2)$$

Each school gets its own latent "noisiness" factor λᵢ, which the model infers from how surprising that school is. The resulting effective weight each school receives in the fit is approximately

$$w_i = \frac{\nu + 1}{\nu + (r_i/\sigma)^2}$$

where $r_i/\sigma$ is the school's standardized residual. With ν = 4:

| Standardized residual | Effective weight (t4) | Weight under normal |
|-----------------------|-----------------------|---------------------|
| 0                     | 1.25                  | 1                   |
| 1                     | 1.00                  | 1                   |
| 2                     | 0.63                  | 1                   |
| 3                     | 0.38                  | 1                   |
| 6                     | 0.13                  | 1                   |
| 15 (our worst school) | 0.02                  | 1                   |

Two points stand out:

1.  **The model decides the weights itself, and jointly with everything else.** You do not choose a cutoff. Uncertainty about which schools are outliers is propagated into the uncertainty of β and the variance components. This is the main conceptual advantage over manual weighting (§3).
2.  **Typical schools are essentially unaffected.** Since 99% of residuals sit within ±3σ, the t4 model behaves like the Gaussian model for the bulk of the data. It only stops the ≈ 100–400 extreme schools from steering the fit.

### 2.3 What it would help with, issue by issue

| Issue | How t4 helps |
|----|----|
| Q-Q tails at level 1 | The likelihood now *expects* tails of roughly this shape, so the model is better specified. Diagnostics should use posterior predictive checks or t-quantile Q-Q plots, not normal Q-Q plots. |
| Inflated residual variance | σ is estimated from the bulk of the data, so σ is less dominated by extremes. |
| Fixed effects pulled by extreme schools | Bounded influence. Most relevant to effects estimated from fewer schools: the sector × covariate interactions (the charter × `perecd`, × `n_grades`, × `sesavgall` terms currently significant), `perasn`, and the Mundlak term. |
| District BLUPs contaminated by one wild school | The wild school's downweighting is carried into the district effect, so district-level variance is less inflated by level-1 outliers. |
| Standard errors | Intervals reflect uncertainty about which observations are outliers. Note that `lmer` SEs under misspecified normality are not robust to fat tails. |

### 2.4 Fixing ν = 4 vs. estimating ν

- **Fixing ν = 4** is a deliberate robustness choice. It says "assume tails at least this heavy" and is easy to explain and reproduce.
- **Estimating ν** lets the data choose the tail weight. `brms` puts a prior on ν by default. A crude moment-based estimate from our residuals' kurtosis (excess kurtosis = 6/(ν−4) for a t-distribution) gives ν ≈ 5.4, so the data are roughly in t4–t6 territory. The kurtosis estimate is itself dominated by the outliers, though, so treat it as a ballpark figure only.

Technical note: a t4 distribution has *infinite* theoretical kurtosis (kurtosis is only finite for ν \> 4). So t4 is slightly heavier-tailed than what the sample moments suggest, which is fine for a conservative robustness check.

**Suggestion:** fit ν = 4 as the headline robustness model and an estimated-ν model as a secondary check. If both agree with each other and differ from the Gaussian model, the conclusion is about the tails, not about the specific ν.

### 2.5 Interpretation pitfall: σ is no longer the residual SD

In a t-distribution, σ is a **scale** parameter, not the standard deviation. The SD is $\sigma\sqrt{\nu/(\nu-2)}$, which for ν = 4 is **√2 · σ ≈ 1.41σ**.

This matters directly for RQ1a (variance partition) and the pseudo-R² calculations in the `.Rmd`. Comparing the t4 model's σ² with `lmer`'s σ² compares different things. Either convert to the implied variance (2σ² for ν = 4) or compare the *fixed effects and level-2/3 variances*, which keep their meaning, and treat the level-1 term separately.

Also, "share of variance at level 1" is less meaningful when level 1 is heavy-tailed, because a few extreme schools dominate any variance. That is a substantive point worth stating in the write-up.

### 2.6 What t4 does *not* fix

1.  **TPS district BLUP tails.** `brms`' standard Student-t family changes only the *level-1* distribution. Random effects remain normal. `brms` does allow Student-t random effects through the `gr()` term, but whether it applies to our correlated two-dimensional sector-specific structure should be checked in the documentation before relying on it. Otherwise, after the t4 fit, re-examine the district effects: some of the TPS district tail may disappear once level-1 outliers stop contaminating district means.
2.  **Bad data.** If a learning rate of −2 GCS per grade comes from a data artifact, downweighting it is better than letting it drive the fit, but the right response is still to look at those schools: n_grades, sector switches, enrollment, and whether the school existed for all years. Robust estimation and outlier investigation complement each other.
3.  **Known sampling error.** t4 treats outliers symmetrically regardless of their SE. Using the SEs is a separate step (§2.8, §3).
4.  **Heteroscedasticity by sector.** Charter and TPS residuals differ in shape and spread (SD 0.211 vs. 0.196). A single σ and a single ν cannot represent that. See §2.8.

### 2.7 Feasibility

**Software.** `brms` and `cmdstanr` are **not currently installed** on this machine; `rstan` is. `brms` also needs a working C++ toolchain to compile each model. This is routine setup on macOS, but budget time for it.

**Computation is the real constraint.**

- The model has ≈ 23,000 district random effects (two per district), 51 state effects, \~25 fixed effects, the covariance parameters, and ν.
- HMC on 54k observations with this structure is feasible, but expect **hours rather than minutes** per model with 4 chains.
- Heavy-tailed likelihoods sometimes sample less efficiently. With ν estimated, ν and σ can trade off, which slows mixing.

Practical levers to discuss before fitting:

- Use `cmdstanr` as the backend; it is faster and more memory-efficient than `rstan`.
- Use within-chain parallelization (`brms` supports threading). This helps a lot with large N.
- **Pilot on a subsample first.** The natural one is the ≈ 930 districts containing both sectors (≈ 21–22k schools), because those districts identify the within-district charter–TPS contrast and the sector covariance anyway. Confirm the pipeline, priors, and diagnostics there, then scale up.
- Save fitted objects to `model-objects/` as you already do for the ecology fits. These fits are too expensive to rerun at every knit.

**Separating two changes at once.** If you go straight from `lmer` (Gaussian, REML) to `brms` (Student-t, Bayesian), any difference could come from the likelihood *or* from the estimation framework and priors. The clean design is:

1.  `brms` with a **Gaussian** family. This should closely reproduce `lmer`; if it doesn't, the priors or sampling need attention.
2.  `brms` with **Student-t, ν = 4**.
3.  Optionally, Student-t with ν estimated.

Differences between (1) and (2) are then attributable to the fat-tail correction alone.

**Model comparison.** Bayesian models are compared with **LOO-CV** (approximate leave-one-out cross-validation, via the `loo` package that `brms` uses). Its **Pareto-k diagnostic** flags individual observations that are highly influential under each model. Under the Gaussian model you would expect many high-k schools among our outliers, and far fewer under t4. That is a concrete, reportable demonstration that the robust likelihood resolves the influence problem.

### 2.8 Extensions available in `brms` that `lmer` can't do

Not needed for the core robustness check, but directly relevant to the research questions:

- **Known measurement error.** `brms` can take each school's SE as known measurement error, so the level-1 variance becomes "true between-school variance + that school's known sampling variance". This is the principled version of RQ1c (signal vs. noise) and of inverse-variance weighting (§3.2), and it can be combined with the t likelihood.
- **Distributional regression.** `brms` can model σ itself (and even ν) as a function of predictors, for example letting the level-1 scale differ by sector. That would let RQ1a's question "are charters more heterogeneous within a district?" be answered at level 1 as well as level 2, within one model.

------------------------------------------------------------------------

## 3. Weighting in `lme4::lmer()`

### 3.1 What the `weights` argument actually does

In `lmer`, `weights` are **precision weights**: school i's residual variance becomes

$$\operatorname{Var}(\varepsilon_i) = \sigma^2 / w_i$$

A school with weight 0.5 is treated as twice as noisy as a school with weight 1, so it pulls less on the fixed effects, on its district's BLUP, and on the variance components. Things to know:

- **They are not sampling or frequency weights.** They do not mean "this school represents w schools in the population". Survey-style weighting in multilevel models is a different, harder problem.
- **The reported σ is for a school with weight 1.** If weights are not normalized (for example to mean 1), σ̂ changes meaning and is not comparable to the unweighted model's σ̂. The fixed effects are unaffected by rescaling the weights, but variance reporting is.
- **`lmer` treats the weights as fixed and known.** Any uncertainty about the weights (because they were estimated, or derived from residuals) is ignored. The resulting SEs and t-values are too optimistic to an unknown degree.
- **Weighted and unweighted likelihoods are not comparable.** `anova()`, AIC, and BIC between a weighted and an unweighted model, or between models with different weights, are not meaningful. Within a set of models sharing the *same* weights, comparisons are fine.

There are two different reasons to weight, and they answer different questions.

### 3.2 Rationale A: inverse-variance weights from the known SEs (precision)

**The idea.** SEDA gives a standard error for every school's learning rate. Schools whose learning rate is measured imprecisely should count less.

**The subtlety specific to multilevel models.** The naive weight $w_i = 1/SE_i^2$ is the right weight only if *all* level-1 variance were sampling error. Here it is about 11%. The actual variance of a school's residual is

$$\operatorname{Var}(\varepsilon_i) = \underbrace{\tau^2}_{\text{true between-school variance}} + \underbrace{SE_i^2}_{\text{sampling variance}}$$

so the appropriate weight is proportional to $1/(\tau^2 + SE_i^2)$. Because τ² ≈ 0.036 is much larger than a typical SE² ≈ 0.0045, the correct weights vary *modestly* across schools: roughly by a factor of 1.5, from the most to the least precise school. The naive 1/SE² weights vary by a factor of **≈ 196**, which massively overstates how much more informative a precisely measured school is.

Using 1/SE² makes a few large, precisely measured schools dominate. That is not a robustness correction; it changes what the estimate represents.

A practical approach is **two-step (feasible GLS)**:

1.  Estimate τ² from an unweighted fit (σ̂² minus the average SE²).
2.  Form $w_i \propto 1/(\hat\tau^2 + SE_i^2)$.
3.  Refit with those weights and normalize them to mean 1.

`brms`' measurement-error term (§2.8) does the same thing in one model and properly propagates τ² uncertainty.

**What it helps with.**

- RQ1c: it begins to separate signal from noise.
- **Heteroscedasticity.** Charter schools tend to be smaller with wider SE distributions (charter SE spread is wider: IQR 0.038–0.083 vs. 0.044–0.075 for TPS), so precision weighting matters for sector comparisons.
- It **modestly** reduces the outliers' influence, because they have above-average SEs.

**What it doesn't help with.** It does **not** fix the fat tails. Outliers at −15σ are not explained by an SE of 0.1 to 0.14; even the noisiest school's sampling SD is below 1σ. Precision weighting cuts their weight by a factor of at most \~1.5 under the correct weights, compared with \~50 under t4. **Precision weighting and robustness to outliers are different corrections.**

**Quick illustration (exploratory, not a result).** In a one-off refit of the interaction model with naive 1/SE² weights, several quantities moved noticeably:

- The charter district SD rose from 0.073 to 0.104.
- The TPS–charter district correlation fell from 0.44 to 0.17.
- The charter × `perecd` interaction shrank from 0.115 to 0.071.
- The charter × `perecd_dist` interaction (−0.078) went to ≈ 0.

The main charter gap barely moved (0.055 → 0.058). This shows the variance-structure findings (RQ1a) and the interaction findings (RQ1b) are sensitive to *how* schools are weighted, and it is exactly why the choice of weight needs to be principled (τ² + SE², not SE² alone).

### 3.3 Rationale B: robustness weights based on residual size (outliers)

**The idea.** Downweight schools *because their residuals are extreme*. This is the manual counterpart of what t4 does automatically. The procedure is **iteratively reweighted** estimation:

1.  Fit the model.
2.  Compute each school's standardized level-1 residual.
3.  Assign weights with a robust weight function.
4.  Refit, recompute residuals, and repeat until estimates stabilize.

Common weight functions:

| Weight function | Behavior | Character |
|----|----|----|
| **Huber** | Weight 1 inside a cutoff (often ≈ 1.345σ), then decreasing like (cutoff / \|z\|) | Mild; no school is ever discarded |
| **Tukey bisquare** | Smoothly decreasing to exactly 0 beyond a cutoff (often ≈ 4.685σ) | Aggressive; very extreme schools are effectively removed |
| **t-weights** | $(\nu+1)/(\nu + z^2)$, §2.2 | Iterating this *is* fitting a t-likelihood, approximately |

**The deep connection:** iterating `lmer` with t-weights approximately recovers the t4 point estimates. The frequentist weighting approach and the Bayesian t4 approach are two routes to the same target. They differ in what is treated as known and how uncertainty is handled.

**What it helps with.**

- It directly addresses the fat level-1 tails visible in the Q-Q plots.
- It stops the ≈ 100 extreme schools from inflating σ̂, pulling coefficients, and contaminating their districts' BLUPs.
- It stays entirely inside `lmer`, so no new software or long runtimes are needed. This is a clear feasibility advantage.

**Limitations to state.**

- **Circularity.** The weights come from residuals of a model that was itself distorted by the outliers. Iteration mitigates this, but the result depends on the starting fit and on the cutoff, which is an analyst choice.
- **Understated uncertainty.** `lmer` treats the final weights as known, so SEs ignore that the weights were estimated from the same data. Bootstrapping or a dedicated robust method is needed for honest intervals. Dedicated robust mixed-model methods exist (for example the `robustlmm` package, not installed here) that downweight extreme values at *every* level, including random effects, and provide adjusted inference.
- **Level-1 only.** As with t4, weights in `lmer` act on level-1 observations. TPS district-level tails are only affected indirectly.
- **Variance components change meaning.** With residual-based weights, σ̂ describes "well-behaved" schools. Report it that way; it is not the variance of all schools.

### 3.4 Combining rationales A and B

They can be multiplied: precision weight × robustness weight. That gives one model addressing both measurement precision and outliers. The more weighting layers are stacked, though, the more the result depends on analyst choices. For a robustness section it is clearer to show each correction separately and then together.

------------------------------------------------------------------------

## 4. Side-by-side: how the two approaches compare

|   | Student-t (t4) in `brms` | Weights in `lmer` |
|----|----|----|
| What addresses the outliers | The likelihood: tails are part of the model | External weights you construct |
| Who decides how much each school counts | The model, jointly with all parameters | You: weight function, cutoff, iteration |
| Uncertainty about which schools are outliers | Propagated into all intervals | Ignored (weights treated as known) |
| Can use the known SEs | Yes: measurement-error term, can be combined with t | Yes: precision weights $1/(\tau^2+SE^2)$ |
| Effect on level-2/3 tails | Indirect (unless t random effects are used) | Indirect |
| Interpretation of σ | Scale, not SD (SD = √2·σ for ν = 4) | Residual SD for a weight-1 school |
| Model comparison | LOO-CV with Pareto-k influence diagnostics | Only among models with identical weights |
| New learning required | Bayesian workflow, priors, MCMC diagnostics | Little |
| Software and compute | New installs; hours per fit | Already available; seconds to minutes |
| Main risk | Long runtimes, sampler issues, prior sensitivity in the thin parts | Arbitrary cutoffs; overconfident SEs; naive 1/SE² weights distort the estimand |

**They complement each other.** If a residual-weighted `lmer` and a t4 `brms` model agree with each other and both differ from the Gaussian baseline in the same direction, that is strong evidence that the fat tails were affecting the original findings. If neither moves the key estimates, the original findings are robust to the tails. Either outcome is a useful result.

------------------------------------------------------------------------

## 5. What to compare in the robustness section

Fix the quantities that carry the substantive conclusions, then put them in one table across specifications.

**Quantities to track:**

- Charter–TPS within-district gap (`factor(charter)1`), currently ≈ 0.043 (main model) / 0.055 (interaction model)
- Mundlak coefficient (`charter_share_dist_c`), ≈ −0.035
- The sector × covariate interactions and the joint Wald test (RQ1b): especially charter × `perecd` (≈ 0.115), × `n_grades` (≈ −0.011), × `sesavgall` (≈ −0.029)
- Variance components: TPS district SD (≈ 0.069), charter district SD (≈ 0.073), their correlation (≈ 0.44), state SD (≈ 0.056), level-1 scale (≈ 0.203). Keep the t-scale caveat (§2.5) in mind.
- The variance partition and pseudo-R² values reported for RQ1a

**Specifications:**

1.  Baseline `lmer` (current)
2.  `lmer` + precision weights $1/(\hat\tau^2 + SE^2)$
3.  `lmer` + iterated robust weights (Huber as the mild version; optionally bisquare as the aggressive one)
4.  `brms` Gaussian (framework check)
5.  `brms` Student-t, ν = 4
6.  `brms` Student-t, ν estimated (optional)
7.  *(Transparency)* baseline refit excluding the schools beyond ±4σ, after inspecting them

**Criteria for "robust":** the same sign, similar magnitude (judged against substantive meaning in GCS units, not only statistical significance), and overlapping intervals.

Given the diagnostics in §0, expect:

- **The main charter gap** to be fairly stable. Outliers are proportionally spread across sectors, and the precision-weighted quick check barely moved it.
- **The interactions and the sector-specific variance components** to be the most sensitive. They rest on fewer schools and already moved under weighting.

That is where the robustness discussion will be most informative.

------------------------------------------------------------------------

## 6. Summary

- The level-1 residuals are **normal in the bulk and much too heavy in the tails** (98 schools beyond ±4σ vs. ≈ 3 expected; excess kurtosis ≈ 4.2). That is the textbook case for a Student-t likelihood.
- **t4 in `brms`** addresses this by assuming heavy tails in the model. Each school's influence shrinks automatically as it becomes more surprising (weight ≈ 0.38 at 3σ, ≈ 0.02 at 15σ), and that uncertainty is propagated into all estimates. It is feasible, but it requires installing `brms`/`cmdstanr`, learning the Bayesian workflow (priors, convergence diagnostics, posterior predictive checks, LOO), and accepting long runtimes. Pilot on the two-sector districts, and pair it with a Gaussian `brms` fit so the effect of the tails can be separated from the effect of switching frameworks.
- **Weights in `lmer`** can serve two distinct purposes:
  - **Precision weights** from the SEs address measurement noise and RQ1c. They must be $1/(\tau^2 + SE^2)$, not $1/SE^2$, because sampling error is only ≈ 11% of level-1 variance. They do **not** fix the fat tails.
  - **Residual-based robust weights** do target the fat tails and are cheap, but they rely on analyst-chosen cutoffs and give overconfident SEs.
- **Neither approach fixes the heavy-tailed TPS district effects directly**, and neither replaces inspecting the extreme schools. Several have learning rates that look implausible enough to suggest data problems.
- **Run them as a set** (baseline, precision-weighted, robust-weighted, Gaussian `brms`, t4 `brms`) and compare a fixed list of key quantities. The main charter gap will probably hold; the interactions and the sector-specific variance structure are where robustness is genuinely in question.
