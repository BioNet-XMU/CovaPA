# Welch screening followed by one-sided Fisher over-representation analysis.
# The equivalent hypergeometric tail avoids building a table for every pathway.
ora <- function(x, labels, compounds, pathways, background, alpha = 0.05,
                universe_size = length(unique(background))) {
  x <- as.matrix(x)
  groups <- unique(labels)
  if (!is.numeric(x) || any(!is.finite(x)) || length(labels) != nrow(x) ||
      length(compounds) != ncol(x) || anyNA(labels) || length(groups) != 2L ||
      any(table(labels) < 2L)) {
    stop("Expected a finite participant-by-feature matrix and two groups.")
  }
  p <- vapply(seq_len(ncol(x)), function(j) {
    a <- x[labels == groups[1L], j]
    b <- x[labels == groups[2L], j]
    s1 <- var(a) / length(a)
    s2 <- var(b) / length(b)
    delta <- mean(a) - mean(b)
    if (s1 + s2 == 0) return(if (delta == 0) 1 else 0)
    df <- (s1 + s2)^2 / (s1^2 / (length(a) - 1L) + s2^2 / (length(b) - 1L))
    2 * pt(-abs(delta) / sqrt(s1 + s2), df)
  }, numeric(1))
  # Apply BH on records first. Multiple records for one compound count once in enrichment.
  significant <- intersect(background, compounds[p.adjust(p, method = "BH") < alpha])
  if (universe_size < length(unique(background)) ||
      !all(unlist(pathways, use.names = FALSE) %in% background)) {
    stop("Every tested pathway member must belong to the enrichment background.")
  }
  members <- lengths(pathways)
  hits <- vapply(pathways, function(ids) length(intersect(ids, significant)), integer(1))
  n_hit <- length(significant)
  enriched <- phyper(hits - 1L, n_hit, universe_size - n_hit, members, lower.tail = FALSE)
  data.frame(pathway = names(pathways), hits = unname(hits), selected = n_hit,
             background = universe_size, p = unname(enriched),
             q = unname(p.adjust(enriched, method = "BH")))
}
