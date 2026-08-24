
############################################################################################################################################################################
#----------------------------------------monocle analysis--------------------------------#
############################################################################################################################################################################
filt_NBRC_altas_seu=readRDS('20251125_filt_NBRC_altas_seu.rds')
filt_NBRC_altas_seu@meta.data=readRDS('20251125_filt_NBRC_altas_seu_meta.rds' )

dir.create('NRBC_monocle_analysis')
dir.create('NRBC_monocle_analysis/res_data')
dir.create('NRBC_monocle_analysis/res_pic')

################################################## FL NRBC ################################################## 
FL_NRBC_seu=subset(filt_NBRC_altas_seu,tissue_stage %in% 'FL')

sample_cells=c( sample(rownames(FL_NRBC_seu@meta.data[FL_NRBC_seu$final_celltype %in% c( 'BFUE/CFUE') ,]),size = 1500),
                sample(rownames(FL_NRBC_seu@meta.data[FL_NRBC_seu$final_celltype %in% c( 'ProE') ,]),size = 1500),
                sample(rownames(FL_NRBC_seu@meta.data[FL_NRBC_seu$final_celltype %in% c( 'Bas') ,]),size = 1500),
                sample(rownames(FL_NRBC_seu@meta.data[FL_NRBC_seu$final_celltype %in% c( 'Poly') ,]),size = 1500),
                sample(rownames(FL_NRBC_seu@meta.data[FL_NRBC_seu$final_celltype %in% c( 'Orth') ,]),size = 1500)
)

mono_FL_NRBC_seu=subset(FL_NRBC_seu,cells=sample_cells)
table(mono_FL_NRBC_seu$final_celltype)
rm(FL_NRBC_seu);gc()


pd <- new("AnnotatedDataFrame", data =  mono_FL_NRBC_seu@meta.data[,c(1:10,27:28)])
fd <- new("AnnotatedDataFrame", data = data.frame(row.names = rownames(mono_FL_NRBC_seu),gene_short_name=rownames(mono_FL_NRBC_seu)))

cds <- newCellDataSet(GetAssayData(mono_FL_NRBC_seu,assay = 'RNA',layer = 'counts'),
                      phenoData =pd,
                      featureData = fd,
                      expressionFamily=negbinomial.size())


cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)



Idents(mono_FL_NRBC_seu)='final_celltype'
marker_gens=FindAllMarkers(mono_FL_NRBC_seu)
top500_marker_gens=marker_gens %>% group_by( cluster) %>% do(head(.,500))
length(unique(top500_marker_gens$gene)) # 
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10))
diff_test_res <- differentialGeneTest(cds[expressed_genes,],cores = 6,
                                      fullModelFormulaStr = "~final_celltype")
ordering_genes=diff_test_res[order(-log10(diff_test_res[,'pval']),decreasing = T),]

cds <- setOrderingFilter(cds, ordering_genes = unique(c(ordering_genes$gene_short_name[1:500],top500_marker_gens$gene)))
plot_ordering_genes(cds)


cds <- reduceDimension(cds, max_components = 2,method = 'DDRTree')
cds <- orderCells(cds)
cds@phenoData@data[,'stage']=filt_NBRC_altas_seu@meta.data[rownames(cds@phenoData@data),'age']

FL_cds=cds

p=plot_cell_trajectory(  FL_cds, color_by = "final_celltype",cell_size = 1)+scale_color_manual(values = cols) +ggtitle('FL NRBC')+
  plot_cell_trajectory(FL_cds, color_by = "Pseudotime",cell_size = 1)+scale_color_gradient(low ='navy',high ='firebrick3')+
  plot_cell_trajectory(FL_cds, color_by = "State",cell_size = 1)+scale_color_manual(values = cols) # 

p=plot_cell_trajectory(  FL_cds, color_by = "final_celltype",cell_size = 1)+scale_color_manual(values = cols) +ggtitle('FL NRBC')+
  plot_cell_trajectory(FL_cds, color_by = "Pseudotime",cell_size = 1)+scale_color_gradient(low ='navy',high ='firebrick3')
p
ggsave(p,filename='NRBC_monocle_analysis/res_pic/FL_nRBC_celltype_pseudotime_monocle.pdf',width = 10,height = 5)

saveRDS(FL_cds,file = 'NRBC_monocle_analysis/sub_FL_NRBC_cds.rds')

mono_FL_NRBC_seu@meta.data[,c('State','Pseudotime')]=FL_cds@phenoData@data[rownames(mono_FL_NRBC_seu@meta.data),c('State','Pseudotime')]
pse_pathway_df=mono_FL_NRBC_seu@meta.data[,c('final_celltype','Pseudotime',cho_pathways)]
pse_pathway_df=melt(pse_pathway_df,id.vars = c('final_celltype','Pseudotime'))
colnames(pse_pathway_df)=c("celltype" ,"Pseudotime","signature","scores" )

p2=ggplot(pse_pathway_df[pse_pathway_df$signature %in% c( "oxygen_transport_UCell","stress_response_to_metal_ion_UCell","regulation_of_autophagy_UCell","macroautophagy_UCell","vacuole_organization_UCell"),],aes(x=Pseudotime,y=scores,colour=signature))+theme_classic()+scale_color_manual(values = cols)+
  geom_smooth(method='loess',formula = y~x,se = F,size=2)+ggtitle('FL_NRBC')
p2

################################################## FBM NRBC ################################################## 
FBM_NRBC_seu=subset(filt_NBRC_altas_seu,tissue_stage %in% 'FBM')

sample_cells=c( sample(rownames(FBM_NRBC_seu@meta.data[FBM_NRBC_seu$final_celltype %in% c( 'BFUE/CFUE') ,]),size = 769), # 769 most
                sample(rownames(FBM_NRBC_seu@meta.data[FBM_NRBC_seu$final_celltype %in% c( 'ProE') ,]),size = 1500),
                sample(rownames(FBM_NRBC_seu@meta.data[FBM_NRBC_seu$final_celltype %in% c( 'Bas') ,]),size = 1500),
                sample(rownames(FBM_NRBC_seu@meta.data[FBM_NRBC_seu$final_celltype %in% c( 'Poly') ,]),size = 1500),
                sample(rownames(FBM_NRBC_seu@meta.data[FBM_NRBC_seu$final_celltype %in% c( 'Orth') ,]),size = 1500)
)

mono_FBM_NRBC_seu=subset(FBM_NRBC_seu,cells=sample_cells)
table(mono_FBM_NRBC_seu$final_celltype)
rm(FBM_NRBC_seu);gc()


pd <- new("AnnotatedDataFrame", data =  mono_FBM_NRBC_seu@meta.data[,c(1:10,27:28)])
fd <- new("AnnotatedDataFrame", data = data.frame(row.names = rownames(mono_FBM_NRBC_seu),gene_short_name=rownames(FBM_NRBC_seu)))

cds <- newCellDataSet(GetAssayData(mono_FBM_NRBC_seu,assay = 'RNA',layer = 'counts'),
                      phenoData =pd,
                      featureData = fd,
                      expressionFamily=negbinomial.size())


cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)



Idents(mono_FBM_NRBC_seu)='final_celltype'
marker_gens=FindAllMarkers(mono_FBM_NRBC_seu)
top500_marker_gens=marker_gens %>% group_by( cluster) %>% do(head(.,500))
length(unique(top500_marker_gens$gene)) # 
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10))
diff_test_res <- differentialGeneTest(cds[expressed_genes,],cores = 6,
                                      fullModelFormulaStr = "~final_celltype")
ordering_genes=diff_test_res[order(-log10(diff_test_res[,'pval']),decreasing = T),]

cds <- setOrderingFilter(cds, ordering_genes = unique(c(ordering_genes$gene_short_name[1:500],top500_marker_gens$gene)))
plot_ordering_genes(cds)


cds <- reduceDimension(cds, max_components = 2,method = 'DDRTree')
cds <- orderCells(cds)
cds@phenoData@data[,'stage']=filt_NBRC_altas_seu@meta.data[rownames(cds@phenoData@data),'age']


plot_cell_trajectory(cds, color_by = "final_celltype",cell_size = 1)+scale_color_manual(values = cols) +ggtitle('FBM NRBC')+
  plot_cell_trajectory(cds, color_by = "Pseudotime",cell_size = 1)+scale_color_gradient(low ='navy',high ='firebrick3')+
  plot_cell_trajectory(cds, color_by = "State",cell_size = 1)+scale_color_manual(values = cols) # 部分 YS  停止在bas 时期，部分Bas继续分化，


FBM_cds=cds
saveRDS(FBM_cds,file = 'NRBC_monocle_analysis/sub_FBM_NRBC_cds.rds')

p=plot_cell_trajectory(  FBM_cds, color_by = "final_celltype",cell_size = 1)+scale_color_manual(values = cols) +ggtitle('FBM NRBC')+
  plot_cell_trajectory(FBM_cds, color_by = "Pseudotime",cell_size = 1)+scale_color_gradient(low ='navy',high ='firebrick3')
p
ggsave(p,filename='NRBC_monocle_analysis/res_pic/FBM_nRBC_celltype_pseudotime_monocle.pdf',width = 10,height = 5)


mono_FBM_NRBC_seu@meta.data[,c('State','Pseudotime')]=FBM_cds@phenoData@data[rownames(mono_FBM_NRBC_seu@meta.data),c('State','Pseudotime')]
mono_FBM_NRBC_seu@meta.data[,cho_pathways]=filt_NBRC_altas_seu@meta.data[rownames(mono_FBM_NRBC_seu@meta.data),cho_pathways]
pse_pathway_df=mono_FBM_NRBC_seu@meta.data[,c('final_celltype','Pseudotime',cho_pathways)]
pse_pathway_df=melt(pse_pathway_df,id.vars = c('final_celltype','Pseudotime'))
colnames(pse_pathway_df)=c("celltype" ,"Pseudotime","signature","scores" )

p3=ggplot(pse_pathway_df[pse_pathway_df$signature %in% c( "oxygen_transport_UCell","stress_response_to_metal_ion_UCell","regulation_of_autophagy_UCell","macroautophagy_UCell","vacuole_organization_UCell"),],aes(x=Pseudotime,y=scores,colour=signature))+theme_classic()+scale_color_manual(values = cols)+
  geom_smooth(method='loess',formula = y~x,se = F,size=2)+ggtitle('FBM_NRBC')

p3

################################################## ABM NRBC ################################################## 
ABM_NRBC_seu=subset(filt_NBRC_altas_seu,tissue_stage %in% 'ABM')

sample_cells=c( sample(rownames(ABM_NRBC_seu@meta.data[ABM_NRBC_seu$final_celltype %in% c( 'BFUE/CFUE') ,]),size = 1500), # 769 most
                sample(rownames(ABM_NRBC_seu@meta.data[ABM_NRBC_seu$final_celltype %in% c( 'ProE') ,]),size = 1162),
                sample(rownames(ABM_NRBC_seu@meta.data[ABM_NRBC_seu$final_celltype %in% c( 'Bas') ,]),size = 1500),
                sample(rownames(ABM_NRBC_seu@meta.data[ABM_NRBC_seu$final_celltype %in% c( 'Poly') ,]),size = 1500),
                sample(rownames(ABM_NRBC_seu@meta.data[ABM_NRBC_seu$final_celltype %in% c( 'Orth') ,]),size = 1500)
)

mono_ABM_NRBC_seu=subset(ABM_NRBC_seu,cells=sample_cells)
table(mono_ABM_NRBC_seu$final_celltype)
rm(ABM_NRBC_seu);gc()


pd <- new("AnnotatedDataFrame", data =  mono_ABM_NRBC_seu@meta.data[,c(1:10,27:28)])
fd <- new("AnnotatedDataFrame", data = data.frame(row.names = rownames(mono_ABM_NRBC_seu),gene_short_name=rownames(mono_ABM_NRBC_seu)))

cds <- newCellDataSet(GetAssayData(mono_ABM_NRBC_seu,assay = 'RNA',layer = 'counts'),
                      phenoData =pd,
                      featureData = fd,
                      expressionFamily=negbinomial.size())


cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)



Idents(mono_ABM_NRBC_seu)='final_celltype'
marker_gens=FindAllMarkers(mono_ABM_NRBC_seu)
top500_marker_gens=marker_gens %>% group_by( cluster) %>% do(head(.,500))
length(unique(top500_marker_gens$gene)) # 
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10))
diff_test_res <- differentialGeneTest(cds[expressed_genes,],cores = 6,
                                      fullModelFormulaStr = "~final_celltype")
ordering_genes=diff_test_res[order(-log10(diff_test_res[,'pval']),decreasing = T),]

cds <- setOrderingFilter(cds, ordering_genes = unique(c(ordering_genes$gene_short_name[1:500],top500_marker_gens$gene)))
plot_ordering_genes(cds)


cds <- reduceDimension(cds, max_components = 2,method = 'DDRTree')
cds <- orderCells(cds)
cds@phenoData@data[,'id']=mono_ABM_NRBC_seu@meta.data[rownames(cds@phenoData@data),'id']
cds@phenoData@data[,'age']=mono_ABM_NRBC_seu@meta.data[rownames(cds@phenoData@data),'age']

plot_cell_trajectory(cds, color_by = "id",cell_size = 1)+scale_color_manual(values = cols)
plot_cell_trajectory(cds, color_by = "age",cell_size = 1)+scale_color_manual(values = cols)


plot_cell_trajectory(cds, color_by = "final_celltype",cell_size = 1)+scale_color_manual(values = cols) +ggtitle('ABM NRBC')+
  plot_cell_trajectory(cds, color_by = "Pseudotime",cell_size = 1)+scale_color_gradient(low ='navy',high ='firebrick3')+
  plot_cell_trajectory(cds, color_by = "State",cell_size = 1)+scale_color_manual(values = cols) # 部分 YS  停止在bas 时期，部分Bas继续分化，

table(cds@phenoData@data[cds@phenoData@data$State=='3','age']) # 主要来自41 岁
table(cds@phenoData@data$age)


ABM_cds=cds
saveRDS(ABM_cds,file = 'NRBC_monocle_analysis/sub_ABM_NRBC_cds.rds')

p=plot_cell_trajectory(  ABM_cds, color_by = "final_celltype",cell_size = 1)+scale_color_manual(values = cols) +ggtitle('ABM NRBC')+
  plot_cell_trajectory(ABM_cds, color_by = "Pseudotime",cell_size = 1)+scale_color_gradient(low ='navy',high ='firebrick3')
p
ggsave(p,filename='NRBC_monocle_analysis/res_pic/ABM_nRBC_celltype_pseudotime_monocle.pdf',width = 10,height = 5)


ABM_cds$Phase =filt_NBRC_altas_seu@meta.data[rownames(ABM_cds@phenoData@data),'Phase']
p=plot_cell_trajectory(ABM_cds, color_by = 'Phase',cell_size = 0.6)+scale_color_manual(values = cols[-2]) +
  plot_cell_trajectory(ABM_cds, markers =  c("CDK1"),use_color_gradient = T,cell_size = 0.6)
p
ggsave(p,filename='NRBC_monocle_analysis/res_pic/ABM_monocle_cellcycle_analysis.pdf',width =10 ,height =5 ,dpi = 300)


mono_ABM_NRBC_seu@meta.data[,c('State','Pseudotime')]=ABM_cds@phenoData@data[rownames(mono_ABM_NRBC_seu@meta.data),c('State','Pseudotime')]
mono_ABM_NRBC_seu@meta.data[,cho_pathways]=filt_NBRC_altas_seu@meta.data[rownames(mono_ABM_NRBC_seu@meta.data),cho_pathways]
pse_pathway_df=mono_ABM_NRBC_seu@meta.data[,c('final_celltype','Pseudotime',cho_pathways)]
pse_pathway_df=melt(pse_pathway_df,id.vars = c('final_celltype','Pseudotime'))
colnames(pse_pathway_df)=c("celltype" ,"Pseudotime","signature","scores" )

p4=ggplot(pse_pathway_df[pse_pathway_df$signature %in% c( "oxygen_transport_UCell","stress_response_to_metal_ion_UCell","regulation_of_autophagy_UCell","macroautophagy_UCell","vacuole_organization_UCell"),],aes(x=Pseudotime,y=scores,colour=signature))+theme_classic()+scale_color_manual(values = cols)+
  geom_smooth(method='loess',formula = y~x,se = F,size=2)+ggtitle('ABM_NRBC')

p4

library(patchwork)
p=(p1+p2+NoLegend()+p3+NoLegend()+p4+NoLegend()) + plot_layout(ncol = 1);p
ggsave(p,filename='NRBC_monocle_analysis/res_pic/niche_NRBC_pathway_with_pseudotime.pdf',width =5 ,height =18 ,dpi = 300)





