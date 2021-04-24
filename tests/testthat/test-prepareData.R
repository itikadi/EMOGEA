# Libraries
library(data.table)

# Select the folder where the inputs and expected outputs are located
rdsPath <- "../testdata"

# Test prepareData
testthat::test_that("Test prepareData",
{
  # Read files which contain the inputs
  expressionData <- fread(file.path(rdsPath, "input_expressionData.csv"))
  metaData <- fread(file.path(rdsPath, "input_metaData.csv"))
  sampleColumn <- "ID"
  conditionColumn <- "condition"


  # Perform prepareData
  results <- prepareData(
    expressionData = expressionData,
    metaData = metaData,
    sampleColumn = sampleColumn,
    conditionColumn = conditionColumn,
  	applyLogTransformation = FALSE)


  # Get the expected results
  expectedResults <- readRDS(file.path(rdsPath, "expected_prepareData.rds"))

  # Check that results are the same as the expected results
  names <- names(expectedResults)

  for (name in names)
  {
    testthat::expect_equal(results[[name]], expectedResults[[name]])
  }
})
