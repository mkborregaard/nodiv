test_that("node_analysis",{
  set.seed(1337)
  res <- Node_analysis(coquettes, 50, "rdtable")
  expect_lt(abs(sum(SOS(res, 28), na.rm = T)-78.52672), 0.0001)
  expect_equal(sum(is.na(SOS(res, 31))), 49)
  
  expect_lt(abs(GND(res, 28) -  0.30157771), 0.0001)
})