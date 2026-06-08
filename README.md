# afterglow

`afterglow` is an R package for half-life-aware interpretation of expression changes. It uses gene-specific mRNA half-lives to infer model-specific expression-equivalent fold-changes from standard RNA-seq, microarray, or other expression matrices.

The package is local-only: it does not transmit user data.

## Install

Install `afterglow` from the local source directory:

```r
install.packages("path/to/afterglow",
                 repos = NULL, type = "source")
```

Or install from a source tarball:

```r
install.packages("afterglow_0.1.0.tar.gz", repos = NULL, type = "source")
```

After the public repository is created, reviewers can install from GitHub:

```r
devtools::install_github("JinluHuyan/afterglow")
```

The core kinetic bracketing and confidence functions use base R dependencies. Publication-style DE workflows using `limma`, `edgeR`, or `DESeq2` require users to install those optional Bioconductor packages in their own R environment; they are not uploaded to or bundled with this repository.

```r
install.packages("BiocManager")
BiocManager::install(c("limma", "edgeR", "DESeq2"))
```

## Input Files

The main workflow needs:

1. An expression matrix with genes as rows and samples as columns.
2. A group table with columns `sample` and `group`.
3. A gene-level half-life table with columns `gene_id` and `half_life_hours`.
4. The elapsed time between stimulus and sampling, in hours.
5. The kinetic model to use: `pulse` or `synthesis_rate`.

Read expression and group CSV files:

```r
library(afterglow)

expr <- read_expression_csv("expr.csv", gene_col = "gene_id")
group <- read_group_csv("group.csv", sample_col = "sample", group_col = "group")
```

You can also provide matrices and data frames directly:

```r
expr <- as.matrix(read.csv("expr.csv", row.names = 1, check.names = FALSE))
group <- data.frame(
  sample = colnames(expr),
  group = c("Control", "Control", "Stim", "Stim")
)
```

## TTDB Half-Life Reference

Download TTDB CSV files from `https://sysbio.gzzoc.com/ttdb/download.html` and keep them together in one directory. The full TTDB CSV snapshot is intentionally not bundled in the public source tarball; configure a local TTDB directory or set `AFTERGLOW_TTDB_DIR` before building a half-life reference. For manuscript reproduction, use the unified input manifest in `写作/resources/manuscript_input_manifest.tsv` to check expected file names and source accessions.

Configure the TTDB directory in R:

```r
set_ttdb_dir("D:/data/ttdb")
ttdb <- read_ttdb_downloads()
```

Or set an environment variable before starting R.

Windows PowerShell:

```powershell
$env:AFTERGLOW_TTDB_DIR = "D:/data/ttdb"
```

macOS/Linux shell:

```bash
export AFTERGLOW_TTDB_DIR=/path/to/ttdb
```

Build a one-row-per-gene half-life table:

```r
half_life_ref <- make_half_life_reference(
  stability_data = ttdb,
  species = "Human",
  condition = "WT",
  aggregate = "median",
  min_r_squared = 0.8
)
```

If your expression matrix uses gene symbols but TTDB has mixed identifier columns, choose the identifier order explicitly:

```r
half_life_ref <- make_half_life_reference(
  stability_data = ttdb,
  species = "Human",
  gene_col_preference = c("gene_name_x", "gene_name_y", "ensembl_id"),
  aggregate = "median"
)
```

You can also supply your own half-life table:

```r
half_life_ref <- data.frame(
  gene_id = c("TNF", "NFKBIA", "IL6"),
  half_life_hours = c(0.5, 1.2, 2.4),
  n_records = c(3, 2, 1)
)
```

## Choosing A Kinetic Model

Use `model = "pulse"` when the design is best interpreted as an instantaneous or transient net RNA response that may decay before sampling. This model is most unstable for short-half-life genes sampled late.

Use `model = "synthesis_rate"` when the design is best interpreted as sustained stimulation or a new ongoing RNA production rate. This model is most unstable for long-half-life genes sampled early.

If a single RNA-seq time point cannot distinguish the temporal response shape, run both models and treat them as a bracketing sensitivity analysis.

```r
pulse_fit <- run_afterglow_de(..., model = "pulse")
synth_fit <- run_afterglow_de(..., model = "synthesis_rate")
```

## Differential Expression Methods

`afterglow` supports common two-group DE methods through `differential_expression()` and `run_afterglow_de()`.

Use `limma` for microarray, chip, qPCR-like, log2 TPM/CPM, or other log-scale expression:

```r
fit <- run_afterglow_de(
  expr = expr,
  group = group,
  half_life_ref = half_life_ref,
  time_hours = 2,
  control_group = "Control",
  treatment_group = "Stim",
  expression_scale = "log2",
  model = "pulse",
  de_method = "limma",
  assay_type = "microarray"
)
```

Use `limma_voom` for RNA-seq raw counts or count-like inferred values:

```r
fit_voom <- run_afterglow_de(
  expr = count_matrix,
  group = group,
  half_life_ref = half_life_ref,
  time_hours = 2,
  control_group = "Control",
  treatment_group = "Stim",
  expression_scale = "linear",
  model = "pulse",
  de_method = "limma_voom",
  assay_type = "rna_seq"
)
```

Use `edgeR` for RNA-seq quasi-likelihood testing:

```r
fit_edgeR <- run_afterglow_de(
  expr = count_matrix,
  group = group,
  half_life_ref = half_life_ref,
  time_hours = 2,
  control_group = "Control",
  treatment_group = "Stim",
  expression_scale = "linear",
  model = "pulse",
  de_method = "edgeR",
  assay_type = "rna_seq"
)
```

Use `deseq2` for integer RNA-seq counts:

```r
fit_deseq2 <- run_afterglow_de(
  expr = count_matrix,
  group = group,
  half_life_ref = half_life_ref,
  time_hours = 2,
  control_group = "Control",
  treatment_group = "Stim",
  expression_scale = "linear",
  model = "pulse",
  de_method = "deseq2",
  assay_type = "rna_seq",
  deseq_round_counts = TRUE
)
```

DESeq2 requires integer counts. If turnover correction creates non-integer count-like values, `deseq_round_counts = TRUE` rounds them and labels the method as `DESeq2_rounded_turnover_counts`; report this approximation. More generally, p-values from DE tests on afterglow-transformed count-like matrices are screening statistics on the model-transformed expression scale. They should not be described as fully calibrated likelihood tests on the original RNA-seq count-generating process, because the per-gene kinetic inversion changes the error structure.

Use `welch` only for quick checks or environments without Bioconductor:

```r
fit_quick <- run_afterglow_de(..., de_method = "welch")
```

## Confidence Scores

Confidence scoring uses three components:

- Model-specific amplification factor, weight 0.60.
- Signal-to-noise ratio from the observed log2 fold-change, weight 0.25.
- Half-life evidence from `n_records`, weight 0.15.

The amplification factor is `exp(lambda * time)` for the pulse model and `1 / (1 - exp(-lambda * time))` for the synthesis-rate model.

Estimate measurement CV from replicates:

```r
cv <- estimate_cv(expr_linear, group)
cv["pooled"]
```

Add confidence to an afterglow result:

```r
ag <- infer_turnover_expression(
  expr = count_matrix,
  group = group,
  half_life_ref = half_life_ref,
  time_hours = 2,
  control_group = "Control",
  treatment_group = "Stim",
  expression_scale = "linear",
  model = "pulse"
)

ag <- add_confidence(ag, cv_measurement = cv["pooled"])
head(ag$confidence)
```

Confidence tiers:

- `High`: amplification factor <= 2.
- `Good`: > 2 and <= 5.
- `Moderate`: > 5 and <= 10.
- `Low`: > 10 and <= 50.
- `Unreliable`: > 50.

For candidate gene prioritization, prefer genes with significant adjusted p-values, meaningful corrected log2 fold-change, `High` or `Good` confidence, and no extreme non-positive-value diagnostics.

## Half-Life Perturbation Sensitivity

Use `half_life_sensitivity()` to test whether conclusions are stable when half-life estimates are perturbed.

```r
sens <- half_life_sensitivity(
  expr = count_matrix,
  group = group,
  half_life_ref = half_life_ref,
  time_hours = 2,
  control_group = "Control",
  treatment_group = "Stim",
  expression_scale = "linear",
  model = "pulse",
  perturbation_sd = c(0, 0.05, 0.10, 0.20),
  n_iter = 50,
  seed = 1
)

sens$summary
write.csv(sens$summary, "afterglow_half_life_sensitivity_summary.csv",
          row.names = FALSE)
write.csv(sens$gene_level, "afterglow_half_life_sensitivity_gene_level.csv",
          row.names = FALSE)
```

Interpretation:

- High `pearson_vs_nominal` and small `median_abs_delta_log2_fc` indicate stable correction.
- Large perturbation sensitivity usually appears in genes with high amplification factors.
- Treat sensitivity as a robustness diagnostic, not as proof that the selected kinetic model is correct.

## Output Columns

`run_afterglow_de()` returns an `afterglow_de` object. The main table is `fit$results`.

Important result columns include:

- `gene_id`: gene identifier.
- `log2_fc`: model-specific corrected log2 fold-change from the DE method.
- `p_value`, `adjusted_p_value`: nominal and FDR-adjusted p-values.
- `method`: DE method actually used.
- `half_life_hours`: matched half-life.
- `model_amplification_factor`: selected model's error amplification factor.
- `pulse_amplification_factor`: pulse-model amplification factor.
- `synthesis_rate_amplification_factor`: synthesis-rate amplification factor.
- `n_non_positive_inferred`: stored in `fit$settings`, not per gene. This is the number of non-positive inferred matrix entries across genes and samples.

Amplification-factor field mapping:

| Field | Where it appears | Meaning |
| --- | --- | --- |
| `amplification_factor` | `fit$confidence` | Model-specific CF used by confidence scoring for the model currently being evaluated. |
| `model_amplification_factor` | `fit$results` after diagnostics are merged | Same selected-model CF, copied into DE results for filtering and reporting. |
| `pulse_amplification_factor` | diagnostics/confidence summaries when both model factors are retained | CF under the pulse model, `exp(lambda * time)`. |
| `synthesis_rate_amplification_factor` | diagnostics/confidence summaries when both model factors are retained | CF under the synthesis-rate model, `1 / (1 - exp(-lambda * time))`. |

Example candidate screen:

```r
fit <- run_afterglow_de(
  expr = count_matrix,
  group = group,
  half_life_ref = half_life_ref,
  time_hours = 2,
  control_group = "Control",
  treatment_group = "Stim",
  expression_scale = "linear",
  model = "pulse",
  de_method = "limma_voom",
  assay_type = "rna_seq"
)

fit <- add_confidence(fit, cv_measurement = 0.05)
res <- merge(fit$results, fit$confidence, by = "gene_id", all.x = TRUE)

candidates <- subset(
  res,
  adjusted_p_value < 0.05 &
    abs(log2_fc) >= 1 &
    confidence_tier %in% c("High", "Good") &
    model_amplification_factor <= 5
)

write.csv(candidates, "afterglow_candidates_high_good.csv", row.names = FALSE)
```

To identify genes gained by afterglow relative to a standard analysis:

```r
standard <- differential_expression(
  mat = log2(count_matrix + 1),
  group = group,
  control_group = "Control",
  treatment_group = "Stim",
  method = "limma",
  assay_type = "rna_seq",
  log_scale = TRUE
)

std_sig <- with(standard, gene_id[adjusted_p_value < 0.05 & abs(log2_fc) >= 1])
ag_sig <- with(res, gene_id[adjusted_p_value < 0.05 & abs(log2_fc) >= 1])
gain