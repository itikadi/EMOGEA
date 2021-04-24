# Select the folder where the inputs and expected outputs are located
rdsPath <- "../testdata"

# Test multivariateCurveResolution
testthat::test_that("Test multivariateCurveResolution",
{
  # Read files which contain the inputs
  prepareDataOutput <- readRDS(file.path(rdsPath, "expected_prepareData.rds"))

  # Perform multivariateCurveResolution
  results <- multivariateCurveResolution(
    expressionMatrix = prepareDataOutput$expressionMatrix,
    residualMatrix = prepareDataOutput$residualMatrix)


  # Get the expected results
  expectedResults <- readRDS(file.path(rdsPath, "expected_multivariateCurveResolution.rds"))

  # Check that results are the same as the expected results
  names <- names(expectedResults)

  for (name in names)
  {
    testthat::expect_equal(results[[name]], expectedResults[[name]])
  }
})

# Test multivariateCurveResolution
testthat::test_that("Test multivariateCurveResolutionNoResiduals",
{
  # Read files which contain the inputs
  prepareDataOutput <- readRDS(file.path(rdsPath, "expected_prepareData.rds"))

  # Perform multivariateCurveResolution
  results <- multivariateCurveResolution(expressionMatrix = prepareDataOutput$expressionMatrix)


  # Get the expected results
  expectedResults <- readRDS(file.path(rdsPath, "expected_multivariateCurveResolutionNoResiduals.rds"))

  # Check that results are the same as the expected results
  names <- names(expectedResults)

  for (name in names)
  {
    testthat::expect_equal(results[[name]], expectedResults[[name]])
  }
})

