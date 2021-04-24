# Select the folder where the inputs and expected outputs are located
rdsPath <- "../testdata"

# Test mlProjection
testthat::test_that("Test mlProjection",
{
  # Read files which contain the inputs
  prepareDataOutput <- readRDS(file.path(rdsPath, "expected_prepareData.rds"))

  # Perform mlProjection
  results <- mlProjection(
    expressionMatrix = prepareDataOutput$expressionMatrix,
    errorCovarianceMatrix = prepareDataOutput$errorCovarianceMatrix)


  # Get the expected results
  expectedResults <- readRDS(file.path(rdsPath, "expected_mlProjection.rds"))

  # Check that results are the same as the expected results
  names <- names(expectedResults)

  for (name in names)
  {
    testthat::expect_equal(results[[name]], expectedResults[[name]])
  }
})
