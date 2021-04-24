#' This function is meant to prepare the input to curve resolution function that calculates the profiles of a
#' gene expression matrix measured as a function of time. The input data is an expression matrix which must
#' have a metadata with replicate information. If there are no replicates, the total least squares implementation
#' of the algorithm cannot be performed.
#'
#' @author Tobias K. Karakach, Federico Taverna
#'
#' @import data.table
#'
#' @param expressionData A data.frame containing the expression data (genes as rows, samples as columns).
#' The first column of the data.frame must contain the gene identifiers/names.
#' Alternatively, the gene identifiers/names can be specified as row names.
#' Duplicated or missing gene identifiers/names will be dropped.
#'
#' @param metaData A data.frame containing the metadata (samples as rows, information as columns).
#' The sample names must match the ones preset in the expression data.
#'
#' @param sampleColumn The name of the metadata column that contains the sample names.
#'
#' @param conditionColumn The name of the metadata column that contains the condition names.
#'
#' @param genesAsRowNames Whether the genes are specified as row names or not.
#'
#' @param applyLogTransformation Whether to apply the log2 transformation to the expression data or not.
#' If the data is already log transformed, the value should be FALSE.
#'
#' @return A list which contains: the processed expression data (expressionMatrix),
#' the residual matrix (residualMatrix), the error covariance matrix (errorCovarianceMatrix).
#'
#' @export
prepareData <- function(expressionData, metaData, sampleColumn = "Label", conditionColumn = "Condition", genesAsRowNames = FALSE, applyLogTransformation = TRUE)
{
  # ----------------------------------------------------
  # Check/process the data and metadata.
  # This is just a preparatory step.
  # ----------------------------------------------------

  # Sanity check for the user input
  if (!any(c("data.frame", "data.table") %in% class(expressionData))) stop ("The expression data must be in the data.frame format.")
  if (!any(c("data.frame", "data.table") %in% class(metaData))) stop ("The metadata must be in the data.frame format.")
  if (!"character" %in% class(sampleColumn)) stop("The sample column must be a string.")
  if (!"character" %in% class(conditionColumn)) stop("The condition column must be a string.")

  # If the expression data and metadata are data.table, convert them as data.frame
  if ("data.table" %in% class(expressionData)) expressionData <- as.data.frame(expressionData)
  if ("data.table" %in% class(metaData)) metaData <- as.data.frame(metaData)

  # Check if the chosen columns are actually present in the metadata
  metadataColumns <- colnames(metaData)
  if (!sampleColumn %in% metadataColumns) stop("sampleColumn was not specified correctly (it is not present in the metadata).")
  if (!conditionColumn %in% metadataColumns) stop("conditionColumn was not specified correctly (it is not present in the metadata).")

  # Subset metadata and change column names
  metaData <- metaData[, c(sampleColumn, conditionColumn)]
  colnames(metaData) <- c("Label", "Condition")

  # Change expression data first column name
  if (genesAsRowNames == TRUE) expressionData <- cbind(Gene = rownames(expressionData), expressionData)
  colnames(expressionData)[1] <- "Gene"

  # Check if the gene columns has any NA, and in case remove them
  naIndex <- which(is.na(expressionData$Gene))

  if (length(naIndex) > 0)
  {
    warning("Some genes are NAs. Removing the NAs.")
    expressionData <- expressionData[-naIndex, ]
  }

  # Check if the genes are unique, if not remove the duplicates
  duplicateIndex <- which(duplicated(expressionData$Gene))

  if (length(duplicateIndex) > 0)
  {
    warning("Duplicate genes detected. Removing the duplicates.")
    expressionData <- expressionData[-duplicateIndex, ]
  }

  # Check if the number/name of samples in the expression data and metadata is the same
  dataSamples <- colnames(expressionData)[-1]
  metadataSamples <- metaData$Label
  if (length(setdiff(dataSamples, metadataSamples)) != 0) stop ("The number/name of samples in the data and metadata are different")

  # Reorder the samples in the metadata to match the samples in the data
  if (!identical(dataSamples, metadataSamples))
  {
    reorderIndex <- match(dataSamples, metadataSamples)
    metaData <- metaData[reorderIndex, ]
  }

  # Convert the condition column to factors
  metaData$Condition <- as.factor(metaData$Condition)

  # Determine replication levels for each condition
  metaData$replicate <- NA

  for (conditionLevel in levels(metaData$Condition))
  {
    levelIndex <- which(metaData$Condition %in% conditionLevel)
    metaData$replicate[levelIndex]<- 1:length(levelIndex)
  }

  # ----------------------------------------------------
  # Calculate the residual and error covariance matrices
  # ----------------------------------------------------

  # Transform the expression data in a matrix
  X <- expressionData
  rownames(X) <- X$Gene
  X <- X[, -1] # Remove the gene column
  X <- as.matrix(X)

  # Transform the 0s to NAs in the expression matrix
  X[X == 0] <- NA

  # Apply the log2 transformation to the expression matrix (if needed)
  if (applyLogTransformation == TRUE) X <- log2(X)

  # Calculate the measurement errors given the replicate information from the metadata
  n <- dim(X)[1]
  m <- dim(X)[2]
  Xres <- matrix(0, nrow = n, ncol = m)
  rownames(Xres) <- rownames(X)
  colnames(Xres) <- colnames(X)

  for (conditionLevel in levels(metaData$Condition))
  {
    levelIndex <- which(metaData$Condition %in% conditionLevel)
    centeredValues <- scale(X[, levelIndex], center = TRUE, scale = FALSE)
    Xres[, levelIndex] <- centeredValues
  }

  # This part differentiates true zeros from those arising from undetected samples
  # where 'NA's are all replaced by zeros and dispropotionately large measurement
  # errors are associated with these values.
  Xres[is.na(X)] <- 999
  X[is.na(X)] <- 0

  # Calculate the error-covariance matrix to use with mlproj
  Xwt <- 1 / Xres
  Xcov <- Xwt %*% t(Xwt)

  # Create the output list
  outputList <- list("expressionMatrix" = X, "residualMatrix" = Xres, "errorCovarianceMatrix" = Xcov)

  return(outputList)
}
