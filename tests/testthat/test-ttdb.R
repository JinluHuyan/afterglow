test_that("make_half_life_reference aggregates TTDB-like rows", {
  stability <- data.frame(
    gene_name_x = c("A", "A", "B"),
    species_name = c("Human", "Human", "Mouse"),
    cell_type = c("HepG2", "HepG2", "mESCs"),
    condition = c("WT", "WT", "WT"),
    half_life = c("1", "3", "10"),
    decay_rate = c("", "", ""),
    stringsAsFactors = FALSE
  )

  ref <- make_half_life_reference(
    stability_data = stability,
    species = "Human",
    tissue = "HepG2",
    condition = "WT",
    aggregate = "median"
  )

  expect_equal(nrow(ref), 1)
  expect_equal(ref$gene_id, "A")
  expect_equal(ref$half_life_hours, 2)
  expect_equal(ref$decay_rate_per_hour, log(2) / 2)
})

test_that("TTDB metadata options and study filters are available", {
  stability <- data.frame(
    gene_name_x = c("A", "A", "B", "C"),
    species_name = c("Human", "Human", "Human", "Mouse"),
    cell_type = c("HepG2", "K562", "HepG2", "mESCs"),
    condition = c("WT", "WT", "stress", "WT"),
    technique = c("Metabolic Labeling", "Metabolic Labeling", "Transcriptional Inhibition", "Metabolic Labeling"),
    study_id = c("s1", "s1", "s2", "s3"),
    sample_id = c("sample1", "sample2", "sample3", "sample4"),
    half_life = c(1, 2, 4, 8),
    stringsAsFactors = FALSE
  )

  opts <- ttdb_filter_options(stability, species = "Human", fields = c("cell_type", "condition"))
  expect_true(all(c("cell_type", "condition") %in% opts$field))
  expect_equal(opts$n_records[opts$field == "cell_type" & opts$value == "HepG2"], 2L)

  ref <- make_half_life_reference(
    stability_data = stability,
    species = "Human",
    cell_type = "HepG2",
    condition = "WT",
    study_id = "s1"
  )
  expect_equal(ref$gene_id, "A")
  expect_equal(ref$half_life_hours, 1)
})

test_that("external TTDB directory option is honored", {
  td <- tempdir()
  old <- set_ttdb_dir(td)
  on.exit({
    options(afterglow.ttdb_dir = old)
  }, add = TRUE)

  expect_equal(normalizePath(get_ttdb_dir(), winslash = "/", mustWork = FALSE),
               normalizePath(td, winslash = "/", mustWork = FALSE))
  expect_equal(length(available_ttdb_files()), 0)
})

test_that("read_ttdb_downloads explains missing external TTDB files", {
  td <- tempdir()
  expect_error(
    read_ttdb_downloads(td, verbose = FALSE),
    "Download TTDB"
  )
})
