# Extract the only macrophage cells 
macro <- subset(scRNA,Cell_type == "Macrophage")

# Finding the differentially expressed genes in macrophages for male and female groups
table(macro@meta.data$gender)
macro.deg <- FindMarkers(macro, ident.1 = "Male", 
                       group.by = 'gender', 
                       ident.2 = "Female")
write.table(macro.deg,file = "macrophage_male_vs_female_deg.xls",sep = '\t',quote = F,col.names = NA)

summary(scRNA@assays$RNA@data["FPR2",])

# Separate the samples based on the expression of FPR2 more or less than 1
highCells=colnames(subset(x = macro, subset = FPR2 >= 1))
highORlow=ifelse(colnames(macro) %in% highCells,'high','low')
table(highORlow)
macro@meta.data$FPR2expr=highORlow
FPR2.deg <- FindMarkers(macro, ident.1 = "high", 
                                group.by = 'FPR2expr', 
                                ident.2 = "low")
write.table(FPR2.deg,file = "macrophage_FPR2pos_vs_FPR2neg_deg.xls",sep = '\t',quote = F,col.names = NA)