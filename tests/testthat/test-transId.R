library(genekitr)

context("transId")

test_that("transId", {
  expect_true("ENSMUSG00000059552" %in% transId(
    id = c("Cyp2c23", "Fhit", "Gal3st2b", "Trp53", "Tp53"),
    trans_to = "ensembl", org = "mouse", unique = TRUE
  )[5])
})
