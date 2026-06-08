test_that("half_life_sensitivity returns reproducible summaries", {
  expr <- matrix(
    c(
      10, 10, 12, 12,
      20, 20, 21, 21,
      15, 15, 18, 18
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("GENE1", "GENE2", "GENE3"), c("C1", "C2", "T1", "T2"))
  )
  group <- data.frame(
    sample = colnames(expr),
    group = c("C", "C", "T", "T"),
    stringsAsFactors = FALSE
  )
  ref <- data.frame(
    gene_id = rownames(expr),
    half_life_hours = c(1, 4, 8),
    n_records = c(3, 2, 1),
    stringsAsFactors = FALSE
  )

  sens <- half_life_sensitivity(
    expr = expr,
    group = group,
    half_life_ref = ref,
    time_hours = 1,
    control_group = "C",
    treatment_group = "T",
    expression_scale = "linear",
    model = "pulse",
    perturbation_sd = c(0, 0.1),
    n_iter = 2,
    seed = 123,
    verbose = FALSE
  )

  expect_s3_class(sens, "afterglow_sensitivity")
  expect_true(all(c("summary", "gene_level", "nominal", "settings") %in% names(sens)))
  expect_equal(sort(sens$summary$perturbation_sd), c(0, 0.1))
  expect_equal(sens$summary$n_genes, c(3, 3))
  expect_true(all(c("nominal_log2_fc", "perturbed_log2_fc", "delta_log2_fc") %in% names(sens$gene_level)))
  expect_equal(ref$half_life_hours, c(1, 4, 8))
})

test_that("half_life_sensitivity validates perturbation settings", {
  expr <- matrix(
    c(10, 10, 12, 12),
    nrow = 1,
    dimnames = list("GENE1", c("C1", "C2", "T1", "T2"))
  )
  group <- data.frame(sample = colnames(expr), group = c("C", "C", "T", "T"))
  ref <- data.frame(gene_id = "GENE1", half_life_hours = 2)

  expect_error(
    half_life_sensitivity(
      expr = expr,
      group = group,
      half_life_ref = ref,
      time_hours = 1,
      control_group = "C",
      treatment_group = "T",
      perturbation_sd = -0.1,
      verbose = FALSE
    ),
    "perturbation_sd"
  )
})
