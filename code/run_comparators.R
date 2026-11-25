# Used by run_analysis.py. Run from the repository root.
source("code/global_test.R")
source("code/ora.R")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 8L)
read_table <- function(path) read.csv(path, check.names = FALSE,
                                      stringsAsFactors = FALSE, fileEncoding = "UTF-8")
gt <- read_table(args[1L])
features <- read_table(args[2L])
mapped <- read_table(args[3L])
members <- read_table(args[4L])
background <- scan(args[5L], what = character(), quiet = TRUE)
sets <- split(members$KEGG_compound, factor(members$KEGG_pathway,
                                          levels = unique(members$KEGG_pathway)))
gt_x <- as.matrix(gt[, -1L, drop = FALSE])
gt_x <- scale(gt_x, center = TRUE, scale = TRUE)
stopifnot(identical(gt$Group, features$Group),
          identical(names(features)[-1L], mapped$Feature))
gt_sets <- lapply(sets, function(ids) match(ids, colnames(gt_x)))
gt_out <- global_test(gt_x, gt$Group, gt_sets, as.integer(args[6L]), as.integer(args[7L]))
ora_out <- ora(features[, -1L, drop = FALSE], features$Group, mapped$Compound,
               sets, background, universe_size = as.integer(args[8L]))
out <- data.frame(pathway = gt_out$pathway, Global_Test_statistic = gt_out$statistic,
                  Global_Test_p = gt_out$p, Global_Test_q = gt_out$q,
                  ORA_hits = ora_out$hits, ORA_selected = ora_out$selected,
                  ORA_background = ora_out$background, ORA_p = ora_out$p, ORA_q = ora_out$q)
write.csv(out, stdout(), row.names = FALSE)
