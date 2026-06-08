if (!requireNamespace("testthat", quietly = TRUE)) {
  message("Package 'testthat' is not installed; skipping testthat tests.")
  q("no")
}

library(testthat)
library(afterglow)

test_check("afterglow")
