test_that("mlpca_c() compute mlpca for case c error
          (independent errors, general heteroscedastic case)", {
  library(RMLPCA)
  data(data_clean)
  data(data_error_c)
  data(sds_c)
  data(data_cleaned_mlpca_c)
  # data that you will usually have on hands
  data_noisy <- data_clean + data_error_c

  # run mlpca_c with rank p = 5
  results <- RMLPCA::mlpca_c(
    X = data_noisy,
    Xsd = sds_c,
    p = 5
  )

  # estimated clean dataset
  data_cleaned_mlpca <- results$U %*% results$S %*% t(results$V)

  expect_equal(data_cleaned_mlpca, data_cleaned_mlpca_c)
})
