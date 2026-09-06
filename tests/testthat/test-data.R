test_that("Vandeputte example data are aligned and ready for IRS", {
  data("vandeputte2017", package = "IRS")

  expect_identical(dim(vandeputte2017$counts), c(166L, 95L))
  expect_identical(
    colnames(vandeputte2017$counts),
    rownames(vandeputte2017$sample_data)
  )
  expect_identical(levels(vandeputte2017$sample_data$group), c("Control", "CD"))
  expect_equal(as.integer(table(vandeputte2017$sample_data$group)), c(66L, 29L))
  expect_true(all(vandeputte2017$sample_data$microbial_load > 0))
})

test_that("precomputed DAA benchmark contains every paired workflow", {
  data("irs_daa_benchmark", package = "IRS")

  expect_equal(nrow(irs_daa_benchmark), 16L)
  expect_setequal(
    unique(irs_daa_benchmark$dataset),
    c("simulation", "vandeputte2017")
  )
  expect_equal(
    unname(as.integer(table(irs_daa_benchmark$dataset))),
    c(8L, 8L)
  )
  expect_true(all(irs_daa_benchmark$irs_f1 >= irs_daa_benchmark$native_f1))
})
