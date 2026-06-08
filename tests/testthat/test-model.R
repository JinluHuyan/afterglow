test_that("pulse model reconstructs expected t0 expression", {
  expr <- matrix(
    c(10, 10, 11, 11),
    nrow = 1,
    dimnames = list("GENE1", c("C1", "C2", "T1", "T2"))
  )
  group <- data.frame(
    sample = c("C1", "C2", "T1", "T2"),
    group = c("C", "C", "T", "T")
  )
  ref <- data.frame(gene_id = "GENE1", half_life_hours = 1)

  out <- infer_turnover_expression(
    expr = expr,
    group = group,
    half_life_ref = ref,
    time_hours = 1,
    control_group = "C",
    treatment_group = "T",
    expression_scale = "linear",
    model = "pulse"
  )

  # half-life = 1 hour and time = 1 hour, so exp(lambda*t)=2.
  # Control baseline B=10, treatment t1=11:
  # inferred t0 = 10 + (11-10)*2 = 12.
  expect_equal(out$corrected_matrix[1, "T1"], 12)
  expect_equal(out$corrected_matrix[1, "T2"], 12)
})

test_that("full workflow returns result table", {
  expr <- matrix(
    c(10, 10.1, 11, 11.2, 20, 20.1, 20.2, 20.3),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("GENE1", "GENE2"), c("C1", "C2", "T1", "T2"))
  )
  group <- data.frame(
    sample = c("C1", "C2", "T1", "T2"),
    group = c("C", "C", "T", "T")
  )
  ref <- data.frame(gene_id = c("GENE1", "GENE2"), half_life_hours = c(1, 10))

  fit <- run_afterglow_de(
    expr = expr,
    group = group,
    half_life_ref = ref,
    time_hours = 1,
    control_group = "C",
    treatment_group = "T",
    expression_scale = "linear",
    model = "pulse",
    de_method = "welch"
  )

  expect_s3_class(fit, "afterglow_de")
  expect_true(all(c("gene_id", "log2_fc", "p_value", "adjusted_p_value") %in% colnames(fit$results)))
})

test_that("synthesis_rate model reconstructs expected equivalent steady-state expression", {
  expr <- matrix(
    c(10, 10, 12, 12),
    nrow = 1,
    dimnames = list("GENE1", c("C1", "C2", "T1", "T2"))
  )
  group <- data.frame(
    sample = c("C1", "C2", "T1", "T2"),
    group = c("C", "C", "T", "T")
  )
  ref <- data.frame(gene_id = "GENE1", half_life_hours = 2)

  out <- infer_turnover_expression(
    expr = expr,
    group = group,
    half_life_ref = ref,
    time_hours = 1,
    control_group = "C",
    treatment_group = "T",
    expression_scale = "linear",
    model = "synthesis_rate",
    verbose = FALSE
  )

  decay <- 2^(-1 / 2)
  expected <- (12 - 10 * decay) / (1 - decay)
  expect_equal(out$corrected_matrix[1, "T1"], expected, tolerance = 1e-8)
  expect_equal(out$diagnostics$model_amplification_factor[1], 1 / (1 - decay), tolerance = 1e-8)
})

test_that("model choice changes amplification diagnostics", {
  expr <- matrix(
    c(10, 10, 12, 12),
    nrow = 1,
    dimnames = list("GENE1", c("C1", "C2", "T1", "T2"))
  )
  group <- data.frame(
    sample = c("C1", "C2", "T1", "T2"),
    group = c("C", "C", "T", "T")
  )
  ref <- data.frame(gene_id = "GENE1", half_life_hours = 4)

  pulse <- infer_turnover_expression(expr, group, ref, 1, "C", "T",
                                     expression_scale = "linear", model = "pulse", verbose = FALSE)
  synthesis <- infer_turnover_expression(expr, group, ref, 1, "C", "T",
                                         expression_scale = "linear", model = "synthesis_rate", verbose = FALSE)

  expect_equal(pulse$diagnostics$model_amplification_factor[1], 2^(1 / 4), tolerance = 1e-8)
  expect_equal(synthesis$diagnostics$model_amplification_factor[1], 1 / (1 - 2^(-1 / 4)), tolerance = 1e-8)
  expect_gt(synthesis$diagnostics$model_amplification_factor[1], pulse$diagnostics$model_amplification_factor[1])
})

test_that("read_expression_csv handles blank top-left header cells", {
  tmp <- tempfile(fileext = ".csv")
  writeLines(c(
    " ,S1,S2",
    "GENE1,1,2",
    "GENE2,3,4"
  ), tmp)

  expr <- read_expression_csv(tmp, verbose = FALSE)

  expect_equal(rownames(expr), c("GENE1", "GENE2"))
  expect_equal(colnames(expr), c("S1", "S2"))
  expect_equal(expr["GENE1", "S1"], 1)
})
