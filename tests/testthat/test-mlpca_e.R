test_that("mlpca_e() compute mlpca for case e error
          (correlated errors, with a different covariance matrix for each row,
          but no error correlation between the rows)", {
  library(RMLPCA)
  data(data_clean_e)
  data(data_error_e)
  # covariance matrix
  data(cov_e)
  data(data_cleaned_mlpca_e)
  # data that you will usually have on hands
  data_noisy <- data_clean_e + data_error_e

  # run mlpca_e with rank p = 1
  results <- RMLPCA::mlpca_e(
    X = data_noisy,
    Cov = cov_e,
    p = 1
  )

  # estimated clean dataset
  data_cleaned_mlpca <- results$U %*% results$S %*% t(results$V)

  expect_equal(data_cleaned_mlpca, data_cleaned_mlpca_e)
})
