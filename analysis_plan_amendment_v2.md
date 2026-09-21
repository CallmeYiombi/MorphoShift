# MorphoShift Analysis Plan Amendment v2

**Project:** MorphoShift  
**Amendment date:** 2026-09-22  
**Project root:** `/home/genetics/yiombi/Morpho`  
**Status:** To be frozen before inspection of remaining cold-benchmark performance results  
**Revision basis:** Prior-art / novelty / confounder assessment supplied in `Claude.pdf`  
**Scope:** Analysis, baseline, uncertainty, ceiling, matched-control, data-provenance, prior-art claim-boundary, and model-selection rules added after the frozen split and phenotype-preprocessing definitions

---

## 1. Purpose

This amendment supplements the already frozen MorphoShift split and phenotype-preprocessing definitions.

The existing split definitions, target preprocessing policy, and benchmark membership are **not changed** by this amendment.

The purpose of this amendment is to predefine:

1. phenotype predictability ceiling analyses;
2. uncertainty estimation and permutation nulls;
3. simple confounder and relational baselines;
4. matched warm/cold and population-matched controls;
5. secondary evaluation metrics;
6. continuous biological novelty analyses;
7. cold-benchmark hyperparameter-selection rules;
8. PPI / graph-access rules;
9. reproducibility and result-versioning requirements;
10. confirmatory decision rules;
11. prior-art-informed novelty claim boundaries;
12. JUMP ORF/CRISPR provenance and technical-confound audits.

The central principle is:

> No primary split, phenotype-preprocessing, metric definition, model-selection rule, or confirmatory decision rule will be changed after inspecting the remaining cold-benchmark performance results.


---

## 1A. Prior-art-informed positioning and claim boundaries

The external prior-art / novelty assessment supplied before remaining cold-result unblinding changes
**how the study is positioned and audited**, but does not change the frozen benchmark membership or
target preprocessing.

The literature assessment identifies several adjacent research directions:

- morphology/image generation from perturbation or transcriptomic conditioning
  (`MorphoDiff`, `MorphDiff`, `IMPA`, `CellFlux`, `FORM`, related image-generation work);
- gene-embedding benchmarks and gene representation learning
  (`GenePert`, `MorphoHELM`, gene-embedding benchmark studies);
- prior-knowledge models for unseen perturbation prediction
  (`GEARS`, `TxPert`, `PRESAGE`, `Stable-Shift`, `PerturbGraph`);
- leakage-aware split design in molecular ML (`DataSAIL` and related leakage-aware benchmarks);
- JUMP morphology used as an input/node attribute for other tasks (`MOTIVE`).

The closest methodological analogue identified in the assessment is `GenePert`, which maps a gene
embedding through a ridge-style predictor to an unseen-gene transcriptomic response. It is therefore
treated as an important architectural comparator in manuscript framing, while the MorphoShift target
is a JUMP Cell Painting morphology profile rather than transcriptomic response.

`MOTIVE` is treated as a directionally different use of JUMP morphology because it uses morphology
profiles as node attributes for inductive link prediction rather than predicting morphology from a
protein/gene representation.

### Allowed novelty framing

The primary contribution may be described conservatively as the **combination** of:

1. protein-language-model representation -> JUMP Cell Painting profile regression for held-out genes;
2. six predefined biological generalization families spanning random, sequence, family, pathway, and
   network novelty;
3. explicit pairwise leakage audits and pathway/network neighborhood purging;
4. construct-level and replicate-level predictability-ceiling analysis;
5. confounder-controlled tests of whether ESM contributes beyond length/amino-acid composition;
6. comparison of sequence-only, relational-shortcut, and later multimodal models under the same
   frozen evaluation framework.

Preferred manuscript wording is:

> "To our knowledge, we systematically evaluate protein-language-model-to-JUMP-morphology profile
> prediction under multi-axis, leakage-aware biological holdouts and explicit phenotype-ceiling
> analysis."

This wording must remain qualified by "to our knowledge" until primary literature verification is
complete.

### Claims that are prohibited without stronger evidence

Do **not** claim:

- first unseen-gene morphology prediction;
- first gene-embedding-to-morphology model;
- first use of biological prior knowledge for unseen perturbation prediction;
- first leakage-aware molecular split;
- state-of-the-art morphology prediction;
- that the proposed multimodal model "beats deep learning" as a class;
- that a cold-benchmark score difference alone proves biological novelty causes the performance drop.

### Literature-verification rule

The supplied novelty assessment contains a mixture of peer-reviewed papers, preprints, repository
material, and search-derived details. Before manuscript submission:

- every novelty-sensitive claim must be checked against the primary paper/repository;
- bibliographic metadata and dataset/version claims must be verified;
- absence-of-prior-art claims must remain conservative even after search;
- this amendment uses the assessment to define **analysis safeguards and claim boundaries**, not as a
  substitute for final literature verification.


---

## 2. Frozen design inherited from the previous benchmark specification

The following definitions remain frozen and unchanged.

### Datasets

- ORF and CRISPR are evaluated as separate prediction tasks.
- ORF target space: standardized 722-dimensional morphology.
- CRISPR target space: standardized 259-dimensional morphology.
- ESM-2 input representation: 1,280-dimensional embedding.
- Primary phenotype transform: TRAIN-only `StandardScaler`.
- Primary PCA target transform: not used.

### Benchmark families

1. `random_gene`
2. `homology50`
3. `homology30`
4. `panther_family`
5. `reactome_pathway`
6. `string_network700`

### Cold-start semantics

- `homology50`: no qualifying TEST–TRAIN/VAL pair at sequence identity >= 50% and bidirectional coverage >= 80%.
- `homology30`: no qualifying TEST–TRAIN/VAL pair at sequence identity >= 30% and bidirectional coverage >= 80%.
- `panther_family`: PANTHER root families do not cross partitions.
- `reactome_pathway`: held-out pathway with pathway-neighborhood purge.
- `string_network700`: held-out STRING >=700 Louvain community with direct-neighbor purge.

### Seed semantics

- `random_gene`, `homology50`, `homology30`, `panther_family`: seeds 42 / 123 / 2024 alter TRAIN/VAL/TEST partitions.
- `reactome_pathway`, `string_network700`: TEST biological challenge is fixed across seeds; seeds alter only TRAIN/VAL allocation in the safe pool.

### Primary baseline metric

The primary metric remains:

\[
R^2_{\mathrm{trainmean}}
=
1 -
\frac{\sum_i \|y_i-\hat{y}_i\|^2}
     {\sum_i \|y_i\|^2}
\]

because standardized target zero corresponds exactly to the TRAIN phenotype mean.

Therefore:

- global-mean baseline: \(R^2_{\mathrm{trainmean}} = 0\);
- positive values improve on the TRAIN-mean predictor;
- negative values are worse than the TRAIN-mean predictor.

---

## 3. Results inspected before this amendment

The following performance results were inspected before freezing this amendment.

### `random_gene`

For both ORF and CRISPR:

- seeds 42 / 123 / 2024;
- global-mean baseline;
- ESM-Ridge baseline;
- Ridge alpha-selection results;
- TEST `R²_trainmean`;
- RMSE;
- gene-wise Pearson;
- feature-wise Pearson.

The Ridge alpha search was expanded and ultimately fixed to:

`0.01, 0.1, 1, 10, 100, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8`

Observed qualitative pattern before this amendment:

- ORF `random_gene`: small positive ESM-Ridge signal across all three seeds.
- CRISPR `random_gene`: approximately TRAIN-mean-level ESM-Ridge performance.

### `homology50`

For both ORF and CRISPR:

- seeds 42 / 123 / 2024;
- global-mean baseline;
- ESM-Ridge baseline;
- Ridge alpha-selection results;
- TEST `R²_trainmean`;
- RMSE;
- gene-wise Pearson;
- feature-wise Pearson.

Observed qualitative pattern before this amendment:

- ORF `homology50`: small positive ESM-Ridge signal remained across three seeds.
- CRISPR `homology50`: ESM-Ridge remained approximately near the TRAIN-mean baseline.


### Label-only / reproducibility observations known before this amendment

The following descriptive ORF observations were also known before the formal ceiling analysis:

- median within-construct replicate-well profile Pearson was approximately `0.107`;
- median same-gene cross-construct profile Pearson was approximately `0.037`.

These values are **not** treated as the formal model ceiling because they use a different descriptive
metric and may be affected by phenotype strength, replicate structure, plate position, batch, and
construct composition.

They serve only as motivation for the metric-aligned construct-shared and split-half ceiling analyses
predeclared below.

### Not inspected before this amendment

Performance results for the following must not be inspected until this amendment is frozen:

- `homology30`
- `panther_family`
- `reactome_pathway`
- `string_network700`

Also not inspected:

- nonlinear ESM models;
- multimodal models;
- PPI models;
- ceiling-adjusted analyses;
- permutation-null analyses;
- confounder-adjusted incremental ESM analyses;
- matched warm/cold controls;
- population-matched random controls.

---

## 4. Immediate freeze rule

This document must be committed and tagged **before opening remaining cold-benchmark performance outputs**.

Recommended Git actions:

```bash
git add analysis_plan_amendment_v2.md
git commit -m "Freeze MorphoShift analysis plan amendment v2"
git tag morphoshift-analysis-plan-v2
```

The commit hash and tag must be copied into downstream result metadata.

---

## 5. Blind QC before performance unblinding

After the current baseline job finishes, the first QC pass must not inspect performance metrics.

### Blind QC may inspect

- expected file existence;
- expected number of preprocessing/model instances;
- expected TRAIN/VAL/TEST counts;
- prediction array shapes;
- target dimensionality;
- ESM dimensionality;
- NaN / Inf counts;
- missing rows;
- duplicated identifiers;
- split/scaler/hash consistency;
- selected alpha identity;
- whether selected alpha lies at the upper/lower grid boundary;
- runtime errors and incomplete jobs.

### Blind QC must not inspect

- `R²_trainmean`;
- RMSE / MAE values;
- Pearson values;
- cosine values;
- retrieval values;
- benchmark-to-benchmark performance differences.

Performance unblinding occurs only after this amendment is committed/tagged and blind QC passes.

---

# PART I — PHENOTYPE CEILING AND RELIABILITY

## 6. ORF construct-shared ceiling proxy

### Source

Do **not** use a gene-aggregated all-construct table if construct identity has already been collapsed.

Use reagent/construct-level consensus data, expected source:

`processed/orf_reagent_consensus.parquet`

together with the reagent/JCP-to-gene and reagent-to-UniProt/Prot_Match mapping used during preprocessing.


### Mandatory `Prot_Match` provenance check

The supplied prior-art assessment notes that `Prot_Match` may be a JUMP metadata field rather than a
variable explicitly defined in the associated morphology paper.

Before applying the `Prot_Match >= 99` ceiling sensitivity rule:

1. identify the exact source metadata file and column;
2. record its data type, range, missingness, and definition if available;
3. store the source-file SHA256;
4. verify that the value refers to the specific ORF construct/reagent being paired;
5. do not infer a biological meaning beyond the metadata definition.

If the provenance/meaning cannot be verified, the `Prot_Match >= 99` analysis is retained only as a
metadata-stratified sensitivity analysis and not described as a sequence-verified ground truth filter.

### Primary construct ceiling population

Include genes with:

- one primary construct defined by the original metadata-only primary-selection rule;
- at least one alternate construct;
- both primary and alternate construct satisfying `Prot_Match >= 99`.

No predictive-performance information may be used to select constructs.

### Pair handling

For a gene with multiple alternates:

- compare the primary construct with every eligible alternate;
- ensure every gene contributes equal total weight regardless of number of alternate constructs.

Sensitivity analysis:

- all eligible pairwise construct combinations;
- again normalized so each gene contributes equal total weight.

### Scaling

For each benchmark instance used for comparison:

- use that instance's frozen TRAIN phenotype scaler;
- apply the scaler unchanged to construct-level phenotypes;
- do not refit a scaler on multi-construct genes.

### Ceiling estimand

For phenotype feature \(j\):

\[
C_j =
\frac{\operatorname{Cov}(Y^{primary}_j, Y^{alternate}_j)}
     {\operatorname{Var}(Y^{primary}_j)}
\]

Report:

1. uniform feature average:
   \[
   C_{\mathrm{feature}} = \frac{1}{p}\sum_j C_j
   \]

2. pooled estimate:
   \[
   C_{\mathrm{pooled}}
   =
   \frac{\sum_j \operatorname{Cov}(Y^{primary}_j,Y^{alternate}_j)}
        {\sum_j \operatorname{Var}(Y^{primary}_j)}
   \]

The primary ceiling proxy for the primary uniform-feature interpretation is `C_feature`.

This quantity is referred to as a:

> construct-shared predictable-variance ceiling proxy

and not as an unconditional mathematical upper bound.

### Bootstrap

Bootstrap unit: gene.

Use 1,000 bootstrap resamples.

Report:

- point estimate;
- percentile 95% CI;
- number of genes;
- number of eligible primary–alternate comparisons.

### Matched model comparison

When comparing Ridge performance with the construct ceiling:

- restrict Ridge evaluation to the same multi-construct genes present in the relevant TEST set;
- recompute Ridge `R²_trainmean` on that subset;
- do not divide a whole-test-set R² by a ceiling estimated from a different subset.

### Ceiling audit trigger

If any model's observed R² exceeds the **upper 95% CI** of the corresponding predeclared ceiling proxy:

- do not interpret this automatically as superior biological prediction;
- trigger mandatory leakage / batch / construct / technical-confounding audit.

---

## 7. ORF and CRISPR split-half technical ceiling

### Purpose

Estimate technical reproducibility independently of construct-level biological reproducibility.

This analysis is mandatory for both ORF and CRISPR.

For CRISPR, where an independent alternate reagent may not be available for every gene, split-half reproducibility is treated as a **technical reliability estimate**, not a gene-invariant biological ceiling.

### Source

Use the clean replicate-well-level profile data underlying the frozen phenotype construction.

The exact input files and columns used must be recorded in analysis metadata.

### Procedure

For every gene/reagent with sufficient replicate wells:

1. randomly divide replicate wells into two approximately equal halves;
2. aggregate each half using the same consensus statistic used by the phenotype pipeline;
3. transform both half-consensus profiles using the appropriate TRAIN-only target scaler;
4. compute the same covariance/variance reliability estimand;
5. repeat random split-half assignment 100 times where possible;
6. summarize the mean split-half reliability;
7. apply Spearman–Brown correction:

\[
\rho_{SB} = \frac{2\rho_{half}}{1+\rho_{half}}
\]

when the underlying reliability statistic lies in the valid correlation interpretation range.

Report both:

- uncorrected split-half reliability;
- Spearman–Brown-corrected technical reliability.

Bootstrap unit: gene/reagent.

Bootstrap resamples: 1,000.

---

## 8. Plate and well-position diagnostics

These analyses accompany the ceiling analyses.

### Required audits

1. **Replicate-position sharing**
   - quantify how often replicate wells for the same reagent occupy the same well position across plates;
   - test whether same-position replicate pairs are more similar than different-position replicate pairs.

2. **Construct plate overlap**
   - quantify whether primary and alternate ORF constructs are profiled on the same plates/batches;
   - compare same-plate and different-plate construct similarity where possible.

3. **Related-gene plate enrichment**
   - test whether genes belonging to the same PANTHER family / close sequence cluster / paralog group are overrepresented on the same plate.

Permutation-based enrichment tests should be used where feasible.

Plate effects are treated as potential technical confounding, not biological signal.

### Additional spatial/batch diagnostics

Also record, when metadata permit:

- row effect;
- column effect;
- edge-vs-center well effect;
- plate ID;
- batch ID;
- plate-map ID.

Evaluate whether phenotype norm and the first major phenotype components differ systematically across
these technical variables.

---

## 8A. Target-feature provenance audit

The prior-art assessment indicates that the public JUMP morphology resource may involve a larger
image-feature space followed by transformations such as sphering/Harmony before the final profile
dimensions used by downstream analyses.

Before interpreting the 722-dimensional ORF and 259-dimensional CRISPR targets as "raw morphology
features", verify the exact provenance of the frozen target tables:

1. source dataset/release;
2. source feature space;
3. aggregation/consensus procedure;
4. sphering, Harmony, PCA, feature selection, or other transformations already applied upstream;
5. plate/batch corrections already applied;
6. whether feature names are original CellProfiler features or transformed coordinates.

The frozen targets are **not changed** by this audit.

Manuscript terminology must match the verified provenance. If the coordinates are transformed,
describe them as processed morphology-profile dimensions rather than raw CellProfiler features.

---

## 8B. ORF library ascertainment and construct-length audit

The prior-art assessment highlights a lentiviral ORF packaging/infection-efficiency issue for long
constructs and reports that the source study removed low-infection-efficiency ORF reagents.

Therefore sequence length may affect both:

- whether an ORF construct is represented in the usable morphology dataset;
- perturbation strength / profile magnitude among retained constructs.

This creates a potential ascertainment/confounding pathway:

`construct length -> packaging/infection/expression -> morphology magnitude / dataset inclusion`.

Before nonlinear model interpretation:

1. recover actual ORF insert length where available;
2. compare insert length with selected UniProt protein length;
3. quantify the relation of length to:
   - availability in the frozen model population;
   - replicate count;
   - `Prot_Match`;
   - phenotype profile norm;
   - replicate reproducibility;
4. record whether long constructs are underrepresented;
5. use construct insert length rather than canonical protein length for the primary ORF length
   confounder whenever construct length is available.

This audit is descriptive and does not alter frozen model-population membership.

---

# PART II — BASELINE UNCERTAINTY AND NULL MODELS

## 9. Resampling units by benchmark

Bootstrap units must reflect split-generating dependencies.

### `random_gene`

Primary resampling unit:

- exact selected UniProt accession group.

Genes sharing the same selected UniProt accession must stay together within a bootstrap draw.

### `homology50` / `homology30`

Primary resampling unit:

- MMseqs split cluster/group used by the frozen partition.

### `panther_family`

Primary resampling unit:

- PANTHER root family.

### `reactome_pathway`

Primary panel summary:

- pathway-macro summary over the 12 fixed held-out pathways.

For a pathway-specific estimate:

- bootstrap TEST genes within that pathway for a pathway-specific CI.

For panel-level variability:

- pathway-level resampling may be reported as a descriptive panel-variability analysis;
- do not interpret this as population-level uncertainty over all possible Reactome pathways.

Primary panel report:

- macro mean across pathways.

Secondary report:

- micro pooled statistic across all pathway test instances, clearly labelled.

### `string_network700`

Because only five fixed communities are included:

- report each community separately;
- use gene-level bootstrap CI within each community;
- panel-level macro mean / median are descriptive summaries only;
- do not use five-community bootstrap as the primary inferential CI.

---

## 10. Bootstrap count

Default:

- 1,000 bootstrap resamples.

All random seeds used for bootstrap must be explicitly stored.

---

## 11. Model comparison bootstrap

Because global mean corresponds exactly to `R²_trainmean = 0`, a separate mean-vs-Ridge paired bootstrap is not necessary.

Paired bootstrap is reserved for comparisons on identical TEST examples, especially:

- `length + AA` vs `length + AA + ESM`;
- ESM-Ridge vs ESM-kNN;
- ESM-Ridge vs ESM-MLP;
- warm vs cold matched controls;
- other models evaluated on the same test set.

Benchmark-to-benchmark comparisons with different TEST genes are **not** treated as paired causal comparisons.

---

## 12. Permutation null

### Null mechanism

Keep ESM/features fixed.

Permute phenotype rows relative to genes.

For each permutation:

- TRAIN phenotype rows are randomly permuted;
- VAL phenotype rows are independently permuted;
- TEST labels remain unchanged and are used only for final evaluation.

### Hyperparameter selection under the null

For each permutation:

1. fit all Ridge alphas on permuted TRAIN labels;
2. select alpha on independently permuted VAL labels;
3. evaluate the selected model on the original TEST labels.

This preserves the full training and alpha-selection procedure under the null.

### Number of permutations

Primary:

- 1,000 permutations.

Permutation p-value:

\[
p =
\frac{1+\#(T_{null}\ge T_{obs})}
     {1001}
\]

### Predeclared permutation scope

Run 1,000 permutations for all instances of:

- `random_gene`;
- `homology50`;
- `homology30`;
- `panther_family`;

for both ORF and CRISPR and all three seeds.

For the fixed biological panels:

- `reactome_pathway`: run the permutation null for every held-out pathway using seed 42 as the primary
  inferential run; seeds 123/2024 are sensitivity runs if computationally feasible;
- `string_network700`: run the permutation null for every held-out community using seed 42 as the
  primary inferential run; seeds 123/2024 are sensitivity runs if computationally feasible.

This scope is fixed before remaining cold-result inspection. If computational limitations prevent a
declared run, the omission and reason are reported; instances are not selected retrospectively based
on observed performance.

---

# PART III — SIMPLE CONFOUNDER AND SEQUENCE BASELINES

## 13. Baseline suite

The following baseline classes are declared before nonlinear modeling.

### A. Non-informative baseline

1. Global TRAIN-mean phenotype.

### B. Simple sequence/statistical covariates

2. Protein/construct length.
3. Amino-acid composition (20 dimensions).
4. Length + amino-acid composition.

### C. Rich sequence representation

5. ESM-Ridge.
6. Length + AA + ESM Ridge.
7. ESM-kNN.

### D. Relational shortcut controls

8. PANTHER-family phenotype mean, where available from TRAIN.
9. Reactome shared-pathway phenotype mean, where available from TRAIN.
10. STRING direct-neighbor phenotype mean / label-propagation baseline, where allowed by the benchmark.
11. STRING degree-only baseline.

### E. Phenotypic-activity / protein-class confounder controls

12. Protein-class indicators where a pre-existing annotation can be frozen independently of performance.
13. Phenotypic-activity / profile-norm predictor or stratification analysis.

These are used to test whether sequence models primarily recover broad protein-class or perturbation-strength
differences rather than detailed morphology-profile structure.

### F. CRISPR-specific contextual controls

14. Chromosome-arm phenotype mean baseline.
15. U2OS expression baseline.
16. Essentiality baseline.

The exact U2OS expression and essentiality source/version must be fixed before those baseline results are inspected.

---

## 14. Length definition

Preferred for ORF:

- actual construct insert length, if recoverable from frozen reagent metadata.

Fallback:

- selected UniProt protein sequence length.

For CRISPR:

- selected UniProt protein sequence length unless a more relevant construct-level quantity is explicitly defined before analysis.

### Length functional form

Primary length model:

- `log(length)` represented using a cubic B-spline;
- 4 degrees of freedom;
- spline basis fit on TRAIN only and applied unchanged to VAL/TEST.

A simple linear `log(length)` model may be reported as secondary.

---

## 15. Amino-acid composition

Use the 20 standard amino-acid fractions.

Composition is computed from the same protein sequence used for the ESM embedding.

No phenotype information is used in feature construction.

---

## 16. Incremental ESM test

The confirmatory sequence-representation test is not ESM vs mean alone.

Primary incremental contrast:

\[
\Delta R^2_{\mathrm{ESM}}
=
R^2_{\mathrm{length+AA+ESM}}
-
R^2_{\mathrm{length+AA}}
\]

Evaluate on identical TEST genes.

Use paired bootstrap with the benchmark-appropriate resampling unit.

### Decision rule

If the 95% paired-bootstrap CI for \(\Delta R^2_{\mathrm{ESM}}\) includes 0:

> ESM-specific predictive information is not distinguished from protein-length and amino-acid-composition effects.

If the lower 95% CI is > 0:

> ESM provides incremental predictive information beyond length and amino-acid composition.

---

## 17. Phenotype magnitude confounding

For every dataset:

- compute target profile norm per gene;
- evaluate association between profile norm and:
  - protein/construct length;
  - amino-acid composition summary variables;
  - selected technical covariates if available.

Report both:

- correlation with phenotype magnitude;
- direction/profile-shape prediction metrics.

This analysis is required before interpreting a nonlinear model as learning rich morphology rather than magnitude effects.

### Protein-class / phenotypic-activity audit

The prior-art assessment notes that broad protein class can be associated with phenotypic activity in
JUMP-style perturbation profiles.

Accordingly:

1. choose one frozen protein-class annotation source before model-performance stratification;
2. report profile norm / reproducibility by class;
3. report class distribution across benchmark TEST sets;
4. test whether ESM prediction performance is concentrated in high-activity classes;
5. where sample size permits, compare ESM incremental \(\Delta R^2\) after controlling for the
   predeclared class indicators.

This analysis is secondary unless the annotation source and complete model form are frozen before
execution.

---

## 18. ESM-kNN

### Representation

- frozen ESM-2 1,280-d embeddings.

### Distance

Primary:

- cosine distance.

### Hyperparameter grid

Predeclared:

`k = {1, 3, 5, 10, 20, 50}`

Select k on VAL only.

Prediction:

- unweighted mean phenotype of the k nearest TRAIN genes.

A distance-weighted sensitivity analysis may be reported separately, but is exploratory unless frozen before execution.

---

# PART IV — MATCHED CONTROLS

## 19. Population-matched random controls

Population-matched random controls are constructed without changing the primary frozen benchmarks.

### Reactome

For every Reactome held-out instance:

- restrict to the same Reactome-mapped population;
- match `n_train`, `n_val`, and `n_test` to the corresponding cold instance;
- construct random gene-disjoint controls.

### STRING

For every STRING community instance:

- restrict to the same STRING-connected population;
- match `n_train`, `n_val`, and `n_test`;
- construct random gene-disjoint controls.

These controls estimate the effect of biological coldness separately from population composition and sample size.

---

## 20. Same-test warm/cold controls

### Reactome / STRING

Keep exactly the same TEST genes.

#### Cold

- use the frozen safe pool;
- exclude purged genes.

#### Warm matched control

- candidate training pool = `safe ∪ purged`;
- retain exactly the same TEST genes;
- sample the same `n_train` and `n_val` as the cold instance.

Use 20 matched warm resamples per instance.

Primary contrast:

\[
\Delta R^2_{\mathrm{warm-cold}}
=
R^2_{\mathrm{warm}}
-
R^2_{\mathrm{cold}}
\]

evaluated on identical TEST genes.

### Homology / PANTHER same-group exposure control

These benchmarks do not have an explicit purge set.

Secondary mechanistic control:

- for each eligible held-out test group with at least two genes:
  - split the group into `test-core` and `warm-exposure` subsets;
  - keep `test-core` fixed for evaluation;
  - allow `warm-exposure` genes only in the warm training condition;
  - randomly remove an equal number of non-group training genes to keep training size fixed.

Singleton groups are excluded from this control.

This is called:

> same-group exposure control

and is not treated as a replacement primary benchmark.

---

# PART V — CONTINUOUS NOVELTY ANALYSIS

## 21. Per-gene novelty variables

For every TEST gene, store where available:

1. maximum sequence identity to any TRAIN protein;
2. alignment query coverage;
3. alignment target coverage;
4. shared Pfam/InterPro domain indicator/count with TRAIN;
5. PANTHER family relation to TRAIN;
6. number of direct STRING >=700 TRAIN neighbors;
7. maximum STRING edge score to TRAIN;
8. STRING degree;
9. annotation density / number of Reactome pathways if used as an exploratory research-bias covariate.

Novelty analyses are secondary/exploratory unless explicitly included in the confirmatory hypotheses below.

Performance-vs-novelty relationships must be assessed at the per-gene level using stored predictions.

---

# PART VI — SECONDARY METRICS

## 22. Primary metric remains unchanged

Primary:

- `R²_trainmean`.

It is not replaced after examining model results.

---

## 23. Secondary metrics

Report:

1. RMSE.
2. MAE.
3. gene-wise profile Pearson.
4. feature-wise across-gene Pearson.
5. gene-wise cosine similarity.
6. normalized retrieval percentile rank.
7. phenotype magnitude error.
8. phenotype direction similarity.
9. TRAIN-PCA-space R².
10. performance stratified by phenotype reproducibility.

---

## 24. Retrieval normalization

Raw MRR or fixed Recall@k are not used as the sole benchmark-comparison metric because candidate-set size varies substantially.

Primary retrieval summary:

\[
\text{normalized rank percentile}
=
\frac{\mathrm{rank}-1}{N-1}
\]

where 0 is best and 1 is worst.

If Recall@k is reported, it is secondary and always accompanied by candidate-set size.

---

## 25. PCA-space sensitivity metric

PCA is not used for primary target training.

It is used only as a secondary evaluation projection.

For each instance:

- fit PCA on TRAIN standardized phenotypes only;
- project VAL/TEST unchanged.

Predeclared component counts:

- ORF: first 309 PCs;
- CRISPR: first 204 PCs.

These counts correspond approximately to the previously audited 90% TRAIN-variance scale and are fixed here for sensitivity analysis.

---

## 26. Reproducibility strata

Training population is not filtered by phenotype reproducibility.

For secondary reporting, TEST genes may be stratified by a predeclared reproducibility measure.

Primary strata, where a reproducibility statistic in [0,1] is available:

- `< 0.10`
- `0.10–0.30`
- `>= 0.30`

If the available reproducibility statistic has a different scale or definition, the strata must be frozen before viewing model performance by stratum.

---

# PART VII — CRISPR-SPECIFIC AUDITS

## 27. CRISPR proximity-bias audit

Before interpreting CRISPR sequence-model performance:

1. determine whether the supplied CRISPR phenotype profiles were corrected for chromosome-arm /
   proximity effects;
2. document the exact correction, source, and processing version if present;
3. specifically verify whether the profile file received a geometric/chromosome-arm correction of the
   type described in prior CRISPR morphology work, rather than assuming generic PCA/Harmony correction
   removes proximity bias;
4. if not clearly corrected, test:
   - phenotype similarity vs genomic distance;
   - same-chromosome-arm similarity;
   - between-arm vs within-arm similarity;
   - chromosome-arm mean phenotype baseline;
5. store chromosome and arm annotations used by the audit.

No assumption is made that a file labelled "PCA corrected", "batch corrected", sphered, or Harmony
corrected is necessarily chromosome-arm corrected without explicit documentation.

The prior-art assessment reports strong within-arm similarity for cpg0016-like CRISPR morphology data;
this motivates the audit but is not imported as a correction factor into MorphoShift.

---

## 28. CRISPR cell-context variables

The following may be introduced as later multimodal inputs / baselines:

- U2OS gene expression;
- gene essentiality in a relevant cell context.

Before use, freeze:

- source;
- release/version;
- preprocessing;
- missing-data handling;
- TRAIN-only transforms;
- model form.

These variables are not retrospectively chosen based on which source yields higher test performance.

---

# PART VIII — GRAPH / PPI RULES

## 29. STRING graph-access policy

Before PPI/GNN training, the primary graph setting is:

> strict inductive cold evaluation.

For a given STRING cold instance:

- TEST genes may not pass messages through PURGED genes to TRAIN genes during training or TEST inference;
- PURGED genes are excluded from model fitting;
- TEST-to-PURGED-to-TRAIN two-hop paths are not usable by the primary model;
- no TEST phenotype labels are available;
- no graph transformation may use TEST phenotype information.

Primary graph for message passing:

- induced graph on active TRAIN/VAL nodes during training;
- TEST nodes are attached only according to a predeclared inductive inference rule that does not introduce paths through PURGED nodes.

A separate transductive external-graph sensitivity analysis may later be reported, but it must be labelled secondary and must not replace the strict inductive primary result.

---

# PART IX — MODEL SELECTION

## 30. Ridge hyperparameter rule

Frozen alpha grid:

`0.01, 0.1, 1, 10, 100, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8`

Selection:

- minimum VAL MSE.

If the optimum is `1e8`:

- record it as an upper-boundary / strongly shrunk solution;
- do not expand the grid post hoc;
- interpret it as approaching the mean predictor.

If multiple alphas are effectively tied within numerical tolerance:

- choose the larger alpha.

Numerical tie tolerance:

- absolute VAL MSE difference <= `1e-8`.

Ridge solver:

- `cholesky`.

---

## 31. Nonlinear ESM model development

Architecture and hyperparameter development occur using `random_gene` only.

Cold TEST results must not be used for architecture selection.

Before the first MLP cold-benchmark run, freeze:

- architecture;
- hidden dimensions;
- dropout grid;
- weight decay grid;
- learning-rate grid;
- batch size;
- optimizer;
- maximum epochs;
- early-stopping rule;
- model-selection metric.

### Primary cold-benchmark selection rule

Preferred strict rule:

- choose architecture/hyperparameters using `random_gene`;
- freeze them for all cold benchmarks;
- cold-benchmark VAL sets may be used only for a predeclared early-stopping rule, not for architecture/hyperparameter search.

A stricter sensitivity analysis may freeze the training epoch count from `random_gene` and avoid cold-VAL-driven early stopping.

Reactome/STRING warm validation must not be used to redesign the model after observing cold TEST results.

### Prior-art comparator interpretation

The nonlinear model must be compared against strong simple baselines because prior perturbation-
prediction literature has repeatedly shown that mean/additive/linear baselines can be difficult to
beat under well-calibrated evaluation.

Therefore:

- a deep model is not considered successful merely because it beats the global mean;
- its primary sequence-model comparison is ESM-Ridge and the incremental length+AA controls;
- its relational comparison is the corresponding family/pathway/STRING shortcut baseline;
- improvement must be evaluated on identical TEST genes with paired uncertainty.

---

# PART X — CONFIRMATORY HYPOTHESES

## 32. Confirmatory hypothesis family

The following are confirmatory.

### H0 — ESM association above the label-permutation null

For each predeclared permutation instance:

\[
H_0:
T_{\mathrm{observed}}
\text{ is exchangeable with the phenotype-permutation null}
\]

where the primary statistic is TEST `R²_trainmean`.

A positive ESM-associated signal is not described as statistically detectable unless the permutation
test rejects the null at the declared multiplicity-adjusted level for the confirmatory family being
interpreted.

This test does not by itself establish ESM-specific information beyond length/AA composition; that is
addressed by H1.

### H1 — Incremental ESM information

For ORF and CRISPR separately:

\[
\Delta R^2 =
R^2_{\mathrm{length+AA+ESM}}
-
R^2_{\mathrm{length+AA}}
\]

Hypothesis:

- ESM contributes predictive information beyond length and amino-acid composition.

Decision:

- support only if the 95% paired-bootstrap CI lower bound is > 0.

### H2 — Nonlinear improvement over linear ESM

\[
\Delta R^2 =
R^2_{\mathrm{ESM-MLP}}
-
R^2_{\mathrm{ESM-Ridge}}
\]

Decision:

- support only if the paired-bootstrap 95% CI lower bound is > 0 on the same TEST set.

### H3 — Biological coldness effect under matched TEST composition

For matched warm/cold controls:

\[
\Delta R^2 =
R^2_{\mathrm{warm}}
-
R^2_{\mathrm{cold}}
\]

Decision:

- a cold-generalization penalty is supported only if the paired-bootstrap 95% CI lower bound is > 0.

### H4 — Relational model improvement beyond relational shortcut baseline

For PPI/pathway models:

- compare learned model against the corresponding declared family/pathway/neighbor baseline on identical TEST genes.

Decision:

- improvement is supported only when the paired-bootstrap 95% CI for model minus shortcut baseline is > 0.

### Multiple-comparison correction

Apply Holm correction within the confirmatory hypothesis family where multiple parallel tests are interpreted jointly.

All other analyses are labelled secondary or exploratory.

---

# PART XI — INTERPRETATION RULES

## 32A. Prior-art-sensitive manuscript claim rules

The following distinctions must be preserved in all abstracts, titles, discussions, and figure
captions.

### Task claim

Allowed:

> protein-LM / sequence representations are evaluated for predicting held-out-gene JUMP morphology
> profiles under frozen leakage-aware splits.

Not allowed:

> first unseen-gene morphology prediction.

Image-generation and perturbation-conditioned morphology-generation work is adjacent prior art even
when its output object differs from profile-level regression.

### Representation claim

Allowed:

> protein-LM-to-profile transfer under the MorphoShift benchmark.

Not allowed:

> first gene-embedding-to-morphology model.

The novelty assessment identified prior gene-embedding/morphology work using non-ESM gene
representations.

### Prior-knowledge claim

Allowed:

> joint evaluation of sequence, family, pathway, and PPI information under the same morphology-profile
> benchmark.

Not allowed:

> first prior-knowledge model for unseen perturbation response.

### Benchmark claim

The strongest novelty claim is the **combined evaluation design**, not any single split primitive:

- random gene;
- homology50;
- homology30;
- PANTHER family;
- Reactome pathway holdout + neighborhood purge;
- STRING community holdout + 1-hop purge;
- explicit pairwise sequence leakage audit.

### Ceiling claim

Allowed after successful audit:

> construct-shared and replicate-level reliability analyses quantify how much predictable phenotype
> signal is available to a sequence-conditioned model.

Avoid "true upper bound" unless the assumptions required for an upper-bound interpretation are
explicitly justified.

---

## 33. Language restrictions

Do not describe a small positive score as:

- "strong";
- "meaningful";
- "reproducible biological signal";
- "generalizes biologically";
- "state of the art";
- "first unseen perturbation predictor";

unless supported by the corresponding predefined uncertainty/null/control analyses.

Allowed descriptive wording before inferential support:

> A small positive predictive association was observed.

---

## 34. Benchmark-difference interpretation

Differences such as:

`random_gene` vs `homology50`

cannot by themselves be interpreted causally as the effect of sequence novelty because TEST composition differs.

Causal-style interpretation of a novelty penalty must primarily rely on:

- same-test matched warm/cold controls;
- population-matched controls;
- continuous novelty analysis.

---

## 35. Ceiling interpretation

If a model exceeds a construct-shared or technical ceiling proxy:

- treat this as an audit signal;
- inspect batch, plate, position, split, metadata, and target leakage;
- do not automatically claim super-ceiling biology.

---

## 36. Metric disagreement

Metric disagreement is expected.

In particular:

- strong shrinkage can produce near-zero `R²_trainmean` while preserving positive cosine/profile-direction agreement;
- magnitude-sensitive and direction-sensitive metrics answer different questions.

No metric will be promoted to primary status after inspecting results.

---

# PART XII — REPRODUCIBILITY AND RESULT VERSIONING

## 37. Never overwrite historical results

Existing results must remain unchanged.

Do not use `--overwrite` to replace the original baseline outputs when generating prediction-preserving reruns.

Use a new result tree, for example:

```text
results/
├── baselines_v1_metrics_only/
├── baselines_v2_predictions/
├── ceiling_v1/
├── null_v1/
├── controls_v1/
└── mlp_v1/
```

---

## 38. Prediction-preserving rerun

After unblinding and blind-QC completion, rerun deterministic Ridge baselines in a new output directory with TEST predictions saved.

The rerun must be compared with the original metric-only results.

Expected comparison:

- identical or numerically equivalent within declared tolerance;
- selected-alpha changes caused only by exact numerical ties must be explicitly recorded.

---

## 39. Required provenance metadata

Every result instance must record, where applicable:

- dataset;
- benchmark;
- seed;
- held-out pathway/community/group ID;
- model;
- git commit hash;
- analysis-plan tag;
- split-file SHA256;
- target-scaler SHA256;
- reagent/construct metadata SHA256 where applicable;
- raw/clean replicate-profile manifest hash where applicable;
- verified upstream morphology-processing provenance/version;
- ESM embedding SHA256 or file-manifest hash;
- code version;
- Python version;
- NumPy version;
- pandas version;
- scikit-learn version;
- PyTorch version for neural models;
- CUDA version for GPU models;
- GPU model;
- random seeds;
- hyperparameters;
- timestamp.

---

## 40. Numerical reproducibility tolerance

CPU Ridge/Cholesky reruns are expected to be effectively deterministic.

For GPU neural models:

- exact bitwise equality is not required;
- report seed-level variation;
- deterministic flags/settings must be stored;
- performance agreement is judged by predefined numerical tolerance and seed distribution, not exact checkpoint bytes.

---

# PART XIII — EXECUTION ORDER FROM THIS AMENDMENT

## 41. Locked execution order

### Phase A — before remaining cold-result inspection

1. Save this amendment.
2. Commit and tag it.
3. Do not inspect remaining cold performance outputs.
4. Freeze a literature/provenance snapshot listing the prior-art papers and dataset-processing claims used by this amendment.
5. Start in parallel:
   - ORF construct-shared ceiling;
   - ORF split-half technical ceiling;
   - CRISPR split-half technical ceiling;
   - plate/row/column/well-position audits;
   - `Prot_Match` provenance verification;
   - target-feature provenance audit;
   - ORF construct-length / ascertainment audit;
   - CRISPR proximity-correction status audit.

### Phase B — after current Ridge jobs finish

6. Run blind QC only.
7. Confirm expected instance/file counts and alpha-boundary metadata.
8. After amendment commit/tag and blind QC, unblind remaining Ridge performance.

### Phase C — baseline inference

9. Rerun Ridge with stored predictions in a new output directory.
10. Verify reproducibility against original results.
11. Compute split-aware bootstrap CIs.
12. Run the predeclared 1,000-permutation nulls.
13. Run:
    - length spline;
    - AA composition;
    - length + AA;
    - length + AA + ESM;
    - ESM-kNN.
14. Evaluate incremental ESM contrast and protein-class/activity sensitivity.

### Phase D — benchmark controls

15. Population-matched random controls.
16. Same-test warm/cold controls.
17. Relational shortcut baselines.
18. Continuous novelty analyses.
19. Secondary metrics.

### Phase E — nonlinear / multimodal modeling

20. Freeze MLP architecture/hyperparameter policy.
21. Train ESM-MLP.
22. Compare against ESM-Ridge and length+AA controls with paired inference.
23. Introduce pathway/PPI/context models only after graph/data-access rules are frozen.
24. Evaluate learned relational models against declared shortcut baselines.

---

## 42. Changes requiring a future amendment

Any later change to the following requires a new versioned amendment:

- primary metric;
- benchmark membership;
- target preprocessing;
- bootstrap unit;
- permutation mechanism;
- ceiling estimand;
- confirmatory hypothesis;
- multiple-testing procedure;
- hyperparameter grid after result inspection;
- graph-access rule;
- matched-control definition;
- MLP model-selection rule;
- novelty-sensitive manuscript claim beyond the boundaries declared here;
- introduction of a new external expression/essentiality/protein-class source after performance inspection.

Future amendments must explicitly list which results had already been inspected at the time of change.

---

## 43. Final status

At the time of this amendment:

- split construction: frozen;
- split QC: passed;
- phenotype preprocessing: frozen;
- phenotype preprocessing QC: passed;
- global-mean and ESM-Ridge baseline implementation: established;
- `random_gene` baseline results: inspected;
- `homology50` baseline results: inspected;
- remaining cold-baseline performance: not to be inspected until this amendment is committed/tagged;
- nonlinear/multimodal model development: not yet started.

**Analysis Plan Amendment v2 is intended to be frozen before remaining cold-result unblinding.**

---

# Appendix A — Changes introduced in v2 from the supplied prior-art / novelty assessment

This revision adds or strengthens the following items before remaining cold-result unblinding:

1. **Prior-art claim boundaries**
   - separates profile regression from image generation;
   - recognizes GenePert as a close embedding->response architectural comparator;
   - recognizes MOTIVE as a reverse-direction JUMP-profile use case;
   - prevents unsupported "first" and SOTA claims.

2. **Construct-ceiling provenance**
   - requires reagent-level source data;
   - requires explicit verification of the `Prot_Match` field;
   - records the previously observed descriptive 0.107 within-construct vs 0.037 cross-construct
     Pearson medians only as motivation, not as the formal ceiling.

3. **JUMP technical-confound audits**
   - plate / row / column / well position;
   - target-feature processing provenance;
   - ORF construct length and library ascertainment;
   - CRISPR chromosome-arm proximity correction status.

4. **Stronger baseline interpretation**
   - deep models must beat appropriate linear/confounder/relational baselines, not only the global mean;
   - ESM-specific evidence is based on incremental `length+AA+ESM` vs `length+AA`.

5. **Protein-class / phenotypic-activity audit**
   - tests whether broad activity/class effects explain apparent sequence predictability.

6. **Permutation scope**
   - declares which benchmark instances receive the 1,000-permutation null before their performance is
     inspected.

7. **Manuscript interpretation safeguards**
   - benchmark novelty is framed as the integrated six-axis leakage-aware evaluation;
   - construct/replicate "ceiling" is described as a ceiling proxy unless stronger assumptions are
     justified.

---

# Appendix B — Prior-art items to verify before manuscript submission

The supplied assessment identifies the following as especially relevant to final literature
positioning. Exact titles, versions, claims, and bibliographic metadata must be verified from primary
sources before citation:

- GenePert;
- MOTIVE;
- MorphoHELM;
- gene-embedding benchmark studies;
- MorphoDiff / MorphDiff / IMPA / CellFlux / FORM and related morphology-generation work;
- GEARS / TxPert / PRESAGE / Stable-Shift / PerturbGraph;
- DataSAIL and related leakage-aware split work;
- Ahlmann-Eltze et al. on simple baselines for perturbation prediction;
- Systema on systematic variation / calibration in perturbation-response evaluation;
- CRISPR proximity-bias work in image-based perturbation datasets;
- JUMP/morphmap data-processing and reagent/QC documentation.

This appendix is a literature-verification checklist, not a statement that every search-derived detail
has already been independently confirmed.

