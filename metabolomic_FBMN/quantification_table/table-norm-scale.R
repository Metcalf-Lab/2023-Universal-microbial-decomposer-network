library(MetabolAnalyze)
setwd("/Users/zacharyburcham/Dropbox/PMI_3_analyses/multi-omics_data/metabolomic/pieters_fbmn_run/quantification_table")

f = read.table("feature_table.txt", header=T, row.names = 1)
f2 = f/sum(f)
f3 = scaling(f2, type="pareto")
write.table(f3, file = "norm-scaled-feature-table.txt", sep = "\t", row.names = T, col.names = T)
