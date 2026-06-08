test_that("DESeq2 rejects log-scale input", {
  mat <- matrix(
    c(10, 10.2, 12, 12.3, 20, 20.1, 21, 21.2),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("GENE1", "GENE2"), c("C1", "C2", "T1", "T2"))
  )
  group <- data.frame(
    sample = colnames(mat),
    group = c("C", "C", "T", "T")
  )

  expect_error(
    differential_expression(
      mat = mat,
      group = group,
      control_group = "C",
      treatment_group = "T",
      method = "deseq2",
      log_scale = TRUE,
      verbose = FALSE
    ),
    "DESeq2 requires count-scale data"
  )
})

test_that("limma rejects RNA-seq count-scale input without log transformation", {
  mat <- matrix(
    c(10, 11, 12, 13, 20, 21, 22, 23),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("GENE1", "GENE2"), c("C1", "C2", "T1", "T2"))
  )
  group <- data.frame(
    sample = colnames(mat),
    group = c("C", "C", "T", "T")
  )

  expect_error(
    differential_expression(
      mat = mat,
      group = group,
      control_group = "C",
      treatment_group = "T",
      method = "limma",
      assay_type = "rna_seq",
      log_scale = FALSE,
      verbose = FALSE
    ),
    "limma with assay_type='rna_seq' expects log-scale"
  )
})

test_that("DESeq2 non-integer turnover counts require explicit rounding approximation", {
  mat <- matrix(
    c(10.2, 10.4, 15.6, 15.8, 20.1, 20.5, 21.3, 21.7),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("GENE1", "GENE2"), c("C1", "C2", "T1", "T2"))
  )
  group <- data.frame(
    sample = colnames(mat),
    group = c("C", "C", "T", "T")
  )

  skip_if_not_installed("DESeq2")
  expect_error(
    differential_expression(
      mat = mat,
      group = group,
      control_group = "C",
      treatment_group = "T",
      method = "deseq2",
      assay_type = "rna_seq",
      log_scale = FALSE,
      deseq_round_counts = FALSE,
      verbose = FALSE
    ),
    "DESeq2 requires integer counts"
  )
})

test_that("count-scale methods reject log-scale input", {
  mat <- matrix(
    c(10, 10.2, 12, 12.3, 20, 20.1, 21, 21.2),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("GENE1", "GENE2"), c("C1", "C2", "T1", "T2"))
  )
  group <- data.frame(
    sample = colnames(mat),
    group = c("C", "C", "T", "T")
  )

  expect_error(
    differential_expression(
      mat = mat,
      group = group,
      control_group = "C",
      treatment_group = "T",
      method = "limma_voom",
      log_scale = TRUE,
      verbose = FALSE
    ),
    "limma_voom requires count-scale data"
  )
  expect_error(
    differential_expression(
      mat = mat,
      group = group,
      control_group = "C",
      treatment_group = "T",
      method = "edgeR",
      log_scale = TRUE,
      verbose = FALSE
    ),
    "edgeR requires count-scale data"
  )
})

test_that("differential_expression validates group table samples", {
  mat <- matrix(
    c(10, 11, 12, 13, 20, 21, 22, 23),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("GENE1", "GENE2"), c("C1", "C2", "T1", "T2"))
  )
  group <- data.frame(
    sample = c("C1", "C2", "T1"),
    group = c("C", "C", "T")
  )

  expect_error(
    differential_expression(
      mat = mat,
      group = group,
      control_group = "C",
      treatment_group = "T",
      method = "welch",
      verbose = FALSE
    ),
    "Group table is missing samples"
  )
})

test_that("limma_voom and edgeR run when Bioconductor dependencies are installed", {
  mat <- matrix(
    c(
      100, 110, 230, 240,
      80, 82, 81, 84,
      50, 52, 110, 112,
      200, 210, 205, 212
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(paste0("GENE", 1:4), c("C1", "C2", "T1", "T2"))
  )
  group <- data.frame(sample = colnames(mat), group = c("C", "C", "T", "T"))

  skip_if_not_installed("limma")
  skip_if_not_installed("edgeR")

  voom <- differential_expression(
    mat = mat,
    group = group,
    control_group = "C",
    treatment_group = "T",
    method = "limma_voom",
    assay_type = "rna_seq",
    log_scale = FALSE,
    verbose = FALSE
  )
  edger <- differential_expression(
    mat = mat,
    group = group,
    control_group = "C",
    treatment_group = "T",
    method = "edgeR",
    assay_type = "rna_seq",
    log_scale = FALSE,
    verbose = FALSE
  )

  expect_true(all(c("gene_id", "log2_fc", "p_value", "adjusted_p_value", "method") %in% names(voom)))
  expect_true(all(c("gene_id", "log2_fc", "p_value", "adjusted_p_value", "method") %in% names(edger)))
  expect_equal(nrow(voom), nrow(mat))
  expect_equal(nrow(edger), nrow(mat))
})
