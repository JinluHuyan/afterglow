test_that("compute_confidence is model-aware", {
  hl <- c(1, 4)

  pulse <- compute_confidence(
    half_life_hours = hl,
    time_hours = 1,
    observed_fc = c(1, 1),
    model = "pulse"
  )
  synthesis <- compute_confidence(
    half_life_hours = hl,
    time_hours = 1,
    observed_fc = c(1, 1),
    model = "synthesis_rate"
  )

  expect_equal(pulse$amplification_factor[1], 2)
  expect_equal(synthesis$amplification_factor[1], 2)
  expect_equal(pulse$amplification_factor[2], 2^(1 / 4), tolerance = 1e-8)
  expect_equal(
    synthesis$amplification_factor[2],
    1 / (1 - 2^(-1 / 4)),
    tolerance = 1e-8
  )
  expect_equal(as.character(pulse$confidence_tier[2]), "High")
  expect_equal(as.character(synthesis$confidence_tier[2]), "Moderate")
  expect_equal(unique(pulse$model), "pulse")
  expect_equal(unique(synthesis$model), "synthesis_rate")
})

test_that("add_confidence records raw and corrected fold changes separately", {
  expr <- matrix(
    c(10, 10, 12, 12),
    nrow = 1,
    dimnames = list("GENE1", c("C1", "C2", "T1", "T2"))
  )
  group <- data.frame(
    sample = c("C1", "C2", "T1", "T2"),
    group = c("C", "C", "T", "T")
  )
  ref <- data.frame(gene_id = "GENE1", half_life_hours = 2, n_records = 4)

  ag <- infer_turnover_expression(
    expr = expr,
    group = group,
    half_life_ref = ref,
    time_hours = 1,
    control_group = "C",
    treatment_group = "T",
    expression_scale = "linear",
    model = "pulse",
    verbose = FALSE
  )
  ag <- add_confidence(ag, cv_measurement = 0.05)

  expect_true(all(c(
    "raw_observed_log2_fc",
    "corrected_log2_fc",
    "confidence_input_log2_fc",
    "amplification_factor",
    "model"
  ) %in% colnames(ag$confidence)))
  expect_equal(ag$confidence$raw_observed_log2_fc, log2(12 / 10), tolerance = 1e-8)
  expect_equal(
    ag$confidence$corrected_log2_fc,
    log2((10 + (12 - 10) * sqrt(2)) / 10),
    tolerance = 1e-8
  )
  expect_equal(ag$confidence$confidence_input_log2_fc, ag$confidence$raw_observed_log2_fc)
})
