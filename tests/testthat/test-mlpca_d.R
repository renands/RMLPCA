test_that("mlpca_d() compute mlpca for case d error
          (commom row covariance matrices)", {
  library(RMLPCA)
  data(data_clean)
  data(data_error_d)
  # covariance matrix
  data(cov_d)
  data(data_cleaned_mlpca_d)
  # data that you will usually have on hands
  data_noisy <- data_clean + data_error_d

  # run mlpca_d with rank p = 2
  results <- RMLPCA::mlpca_d(
    X = data_noisy,
    Cov = cov_d,
    p = 2
  )

  # estimated clean dataset
  data_cleaned_mlpca <- results$U %*% results$S %*% t(results$V)

  expect_equal(data_cleaned_mlpca, data_cleaned_mlpca_d)
})
