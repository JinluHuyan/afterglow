test_that("make_half_life_reference aggregates TTDB-like rows", {
  stability <- data.frame(
    gene_name_x = c("A", "A", "B"),
    species_name = c("Human", "Human", "Mouse"),
    condition = c("WT", "WT", "WT"),
    half_life = c("1", "3", "10"),
    decay_rate = c("", "", ""),
    stringsAsFactors = FALSE
  )

  ref <- make_half_life_reference(
    stability_data = stability,
    species = "Human",
    condition = "WT",
    aggregate = "median"
  )

  expect_equal(nrow(ref), 1)
  expect_equal(ref$gene_id, "A")
  expect_equal(ref$half_life_hours, 2)
  expect_equal(ref$decay_rate_per_hour, log(2) / 2)
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
