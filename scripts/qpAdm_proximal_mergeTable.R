setwd('/project/jnovembre/jhmarcus/ancient-sardinia/output/qpAdm_rev/proximal_ancinds')

tbl = read.table('proximal_1way_model.csv', sep = ',', header = T)
tbl$n_source = 1

for (i in 2:5){
  
  tbl2 = read.table(sprintf('proximal_%sway_model.csv', i), sep = ',', header = T)
  tbl2$n_source = i
  
  tbl = merge(x = tbl, y = tbl2, all = TRUE)
}

col_order <- names(tbl2)
tbl <- tbl[, col_order]
levels(tbl$status)[levels(tbl$status) == '-'] <- 'feasible'

write.table(tbl, 'merged_output.csv', sep = ',' , col.names = T, row.names = F, na = '-', quote = F)
