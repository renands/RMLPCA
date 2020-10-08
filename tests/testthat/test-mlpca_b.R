test_that("mlpca_b() compute mlpca for case b error
          (independent errors, homoscedastic within a column)", {
  library(RMLPCA)
  data(data_clean)
  data(data_error_b)
  data(sds_b)
  data(data_cleaned_mlpca_b)

  # data that you will usually have on hands
  data_noisy <- data_clean + data_error_b

  # run mlpca_b with rank p = 2
  results <- RMLPCA::mlpca_b(
    X = data_noisy,
    Xsd = sds_b,
    p = 2
  )

  # estimated clean dataset
  data_cleaned_mlpca <- results$U %*% results$S %*% t(results$V)

  expect_equal(data_cleaned_mlpca,data_cleaned_mlpca_b)


})
