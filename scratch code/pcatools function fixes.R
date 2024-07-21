# This is from the PCAtools package and was rewritten to function with our data. Here is the code
library(BiocSingular)
library(matrixStats)
library(PCAtools)

pca <- function (mat, metadata = NULL, center = TRUE, scale = FALSE, 
                 rank = NULL, removeVar = NULL, transposed = FALSE, BSPARAM = ExactParam()) 
{
  if (is.data.frame(mat)) {
    mat <- as.matrix(mat)
  }
  if (!transposed) {
    mat <- t(mat)
  }
  if (!is.null(metadata)) {
    if (!identical(rownames(mat), rownames(metadata))) {
      stop("'colnames(mat)' is not identical to 'rownames(metadata)'")
    }
  }
  .center <- if (center) 
    NULL
  else 0
  vars <- colVars(mat, center = .center, useNames = FALSE) # Edits were made to this line.
  if (!is.null(removeVar)) {
    message("-- removing the lower ", removeVar * 100, "% of variables based on variance")
    varorder <- order(vars, decreasing = TRUE)
    keep <- head(varorder, max(1, ncol(mat) * (1 - removeVar)))
    mat <- mat[, keep, drop = FALSE]
    vars <- vars[keep]
  }
  if (is.null(rank)) {
    if (is(BSPARAM, "ExactParam")) {
      rank <- min(dim(mat))
    }
    else {
      stop("'rank' must be specified for approximate PCA methods")
    }
  }
  pcaobj <- runPCA(mat, center = center, scale = scale, rank = rank, 
                   BSPARAM = BSPARAM)
  if (scale) {
    total.var <- length(vars)
  }
  else {
    total.var <- sum(vars)
  }
  proportionvar <- (pcaobj$sdev^2)/total.var * 100
  pcaobj <- list(rotated = data.frame(pcaobj$x), loadings = data.frame(pcaobj$rotation), 
                 variance = proportionvar, sdev = pcaobj$sdev, metadata = metadata, 
                 xvars = colnames(mat), yvars = rownames(mat), components = colnames(pcaobj$x))
  rownames(pcaobj$rotated) <- pcaobj$yvars
  rownames(pcaobj$loadings) <- pcaobj$xvars
  names(pcaobj$variance) <- pcaobj$components
  class(pcaobj) <- "pca"
  return(pcaobj)
}
