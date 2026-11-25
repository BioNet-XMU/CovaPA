# Shared-permutation Global Test for an intercept-only binomial model.
# Adapted from Supporting_Data.zip, meth_global_test/reproduce_meth_global_test.R.
# The archived calculation was checked against Bioconductor globaltest 5.67.0.
# Goeman et al. (2004), https://doi.org/10.1093/bioinformatics/btg382

# x: participant-by-metabolite matrix on the desired, pooled-autoscaled scale.
# labels: two group labels, in the order of the rows of x.
# pathways: a named list of one-based column indices for the measured sets.
# Returns the statistic, Monte Carlo p, and BH q over all supplied pathways.
global_test <- function(x, labels, pathways, permutations = 999L, seed = 20261119L) {
  x <- as.matrix(x)
  if (!is.numeric(x) || any(!is.finite(x)) || length(labels) != nrow(x) ||
      anyNA(labels) || length(unique(labels)) != 2L) {
    stop("Expected a finite numeric matrix and exactly two group labels.")
  }
  if (!is.list(pathways) || !length(pathways) || is.null(names(pathways)) ||
      anyNA(names(pathways)) || any(!nzchar(names(pathways))) || anyDuplicated(names(pathways))) {
    stop("Provide a list of pathways with unique, non-empty names.")
  }
  valid_members <- vapply(pathways, function(indices) {
    is.numeric(indices) && length(indices) > 0L && all(is.finite(indices)) &&
      all(indices >= 1L & indices <= ncol(x) & indices == floor(indices)) &&
      !anyDuplicated(indices)
  }, logical(1))
  if (!all(valid_members)) {
    stop("Pathway members must be distinct, valid one-based column indices.")
  }
  if (length(permutations) != 1L || !is.numeric(permutations) || !is.finite(permutations) ||
      permutations < 1L || permutations != floor(permutations)) {
    stop("The number of permutations must be a positive integer.")
  }

  # Set the same RNG used in the archived METH test and restore the caller's RNG.
  old_seed <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_kind <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (is.null(old_seed)) {
      rm(".Random.seed", envir = .GlobalEnv)
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  })
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)

  y <- as.integer(factor(labels)) - 1L
  centered <- y - mean(y)
  shuffled <- replicate(permutations, centered[sample.int(length(y))])
  group_ss <- sum(centered^2)

  results <- lapply(names(pathways), function(pathway) {
    features <- x[, pathways[[pathway]], drop = FALSE]
    feature_ss <- sum(features^2)
    if (feature_ss == 0) {
      stop(sprintf("Pathway '%s' has no nonzero data on the analysis scale.", pathway))
    }
    score <- function(scores) {
      colSums(crossprod(features, scores)^2) * 100 / (feature_ss * group_ss)
    }
    observed <- score(centered)
    exceedances <- sum(score(shuffled) >= observed * (1 - 1e-14))
    data.frame(pathway = pathway, statistic = observed,
               p = (1 + exceedances) / (permutations + 1))
  })
  result <- do.call(rbind, results)
  result$q <- p.adjust(result$p, method = "BH")
  result
}
