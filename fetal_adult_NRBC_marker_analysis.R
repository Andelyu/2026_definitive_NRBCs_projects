
########################################################################################################################################################################
#-------------------------------------------- fetal NRBC marker#-------------------------------------------- #
########################################################################################################################################################################
setwd('definitive_nRBC_marker_project')
filt_NBRC_altas_seu=readRDS('../20251125_filt_NBRC_altas_seu.rds')
filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu,tissue_stage !='YS')

CC_hDEGs_df=readRDS('res_data/CC_hDEGs_df.rds')
trans_CC_hDEGs_df=CC_hDEGs_df[CC_hDEGs_df$transmembrane=='yes',]

#DLK1: transmenbrane
top_fetal_unique_marker_genes=c('HBG1','HBG2','HBZ', "TUBB6","HSPA1A","HSPA1B",'IGF2BP1','IGF2BP3','DLK1','CISH','HIF3A')# 还可以考虑CHD7,TIMP3,CISH,后面好像有获得
# PDZK1IP1:transmenbrane
ABM_cho_topmarkers=c( 'HBD',"CA1","PDZK1IP1", "ANXA1","NECAB1","ANKRD28","TSC22D3","IFIT1B",'LGALS9',"HLA-B","HLA-DRA","HLA-DRB1") # top10 中挑选基因,

fetal_adult_NRBC_whole_marker=readRDS('Protein_NRBC_marker/res_data/fetal_adult_NRBC_whole_marker.rds')
sub_fetal_adult_all_Ery_tissue_markers=read.csv('Protein_NRBC_marker/DE_marker/fetal_adult_all_Ery_RNA_markers.csv')
sub_fetal_adult_all_Ery_tissue_markers=sub_fetal_adult_all_Ery_tissue_markers[,-1]

#fetal_adult_NRBC_whole_marker$avg_log2FC>1 & fetal_adult_NRBC_whole_marker$pct.1>0.05 & fetal_adult_NRBC_whole_marker$pct.2<0.1
fetal_adult_NRBC_whole_marker[fetal_adult_NRBC_whole_marker$gene=='PDZK1IP1',] # ok
#sub_fetal_adult_all_Ery_tissue_markers$avg_log2FC >1 & sub_fetal_adult_all_Ery_tissue_markers$pct.2 <0.1 & sub_fetal_adult_all_Ery_tissue_markers$pct.1 >0.1 
sub_fetal_adult_all_Ery_tissue_markers[sub_fetal_adult_all_Ery_tissue_markers$gene=='PDZK1IP1' & sub_fetal_adult_all_Ery_tissue_markers$avg_log2FC >0,] # 满足

# check
tmp_df=sub_fetal_adult_all_Ery_tissue_markers[sub_fetal_adult_all_Ery_tissue_markers$gene=='PDZK1IP1' & sub_fetal_adult_all_Ery_tissue_markers$celltype=='early_Ery'& sub_fetal_adult_all_Ery_tissue_markers$avg_log2FC >0,]
tmp_df$transmembrane='yes'
trans_CC_hDEGs_df=rbind(trans_CC_hDEGs_df,tmp_df)

p1=ggplot(trans_CC_hDEGs_df,aes(x=gene,y=avg_log2FC,fill=cluster))+geom_bar(stat = 'identity')+theme_classic()+theme(axis.text.x = element_text(angle = 90,hjust = 1))+NoLegend()+
  geom_hline(yintercept = 1,linetype = "dashed", color = "red", linewidth = 0.4)+NoLegend()

p2=DotPlot(filt_NBRC_altas_seu,group.by = 'source_celltype',features =trans_CC_hDEGs_df$gene,cols = c('white','firebrick3'))+RotatedAxis()

p=p1+p2+plot_layout(ncol = 1,heights = c(0.6,1.2));p # 后续添加上taget genes 中同类满足条件的基因
ggsave(p,filename='res_pic/main_figure5/specific_transmenbrane_profile_expression.pdf',width=8,height = 8)


#-----------------------------------------------------1.1 primary bulk RNAseq level validated ---------------------------------------------#
# ------------------------bulk RNAseq to validated the unique markers-----------------------#
library(ggplot2)
library(reshape2)
library(edgeR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(SparseArray)
library(ggplotify)
library(Seurat)
#------------------------------ABM nRBC  progenitor bulkRNAseq-----------------#
GSE301441_BM_NRBC_df1=read.table('../ref_data/bulk_RNAseq/GSE301441_BM_primary_nRBC/progenitor_processeddata.txt')
GSE301441_BM_NRBC_df=read.table('../ref_data/bulk_RNAseq/GSE301441_BM_primary_nRBC/terminal_processeddata.txt')

all_genes=unique(c(rownames(GSE301441_BM_NRBC_df1),rownames(GSE301441_BM_NRBC_df)))

GSE301441_BM_NRBC_df1=rbind(GSE301441_BM_NRBC_df1,data.frame(row.names = all_genes[!all_genes %in% rownames(GSE301441_BM_NRBC_df1)]))
GSE301441_BM_NRBC_df1[all_genes[!all_genes %in% rownames(GSE301441_BM_NRBC_df1)],1:dim(GSE301441_BM_NRBC_df1)[2]]=0
GSE301441_BM_NRBC_df=rbind(GSE301441_BM_NRBC_df,data.frame(row.names = all_genes[!all_genes %in% rownames(GSE301441_BM_NRBC_df)]))
GSE301441_BM_NRBC_df[all_genes[!all_genes %in% rownames(GSE301441_BM_NRBC_df)],1:dim(GSE301441_BM_NRBC_df)[2]]=0
GSE301441_BM_NRBC_df=cbind(GSE301441_BM_NRBC_df1,GSE301441_BM_NRBC_df);dim(GSE301441_BM_NRBC_df)

GSE301441_BM_NRBC_df=as.matrix(GSE301441_BM_NRBC_df)
GSE301441_BM_NRBC_df=log2(cpm(GSE301441_BM_NRBC_df)+1)
GSE301441_BM_NRBC_df[c('GAPDH','ACTB'),] # ACTB stable， filter out MM
GSE301441_BM_NRBC_df=GSE301441_BM_NRBC_df[,-grep('MM',colnames(GSE301441_BM_NRBC_df))]
GSE301441_BM_NRBC_df[c('GAPDH','ACTB'),] 


#------------------------------FL nRBC bulkRNAseq-----------------#
FL_primary_NRBC_df=read.table('../ref_data/bulk_RNAseq/Anxiuli_lab_bulkRNAseq/FL_nRBC_invivo.txt',sep="\t",header = T)
rownames(FL_primary_NRBC_df)=FL_primary_NRBC_df$X;FL_primary_NRBC_df=FL_primary_NRBC_df[,-1]
FL_primary_NRBC_df=as.matrix(FL_primary_NRBC_df)
FL_primary_NRBC_df=log2(cpm(FL_primary_NRBC_df)+1)
FL_primary_NRBC_df[c('GAPDH','ACTB'),] # ACTB stable, as reference


#----------------------merge-------------------------#
all_genes=unique(c(rownames(GSE301441_BM_NRBC_df),rownames(FL_primary_NRBC_df)))

FL_primary_NRBC_df=data.frame(FL_primary_NRBC_df)
GSE301441_BM_NRBC_df=data.frame(GSE301441_BM_NRBC_df)
FL_primary_NRBC_df[all_genes[!all_genes %in%  rownames(FL_primary_NRBC_df)],1:dim(FL_primary_NRBC_df)[2]]=0
GSE301441_BM_NRBC_df[all_genes[!all_genes %in%  rownames(GSE301441_BM_NRBC_df)],1:dim(GSE301441_BM_NRBC_df)[2]]=0
GSE301441_BM_NRBC_df['DLK1',]
FL_primary_NRBC_df['DLK1',]

FL_ABM_nRBC_df_df=cbind(FL_primary_NRBC_df,GSE301441_BM_NRBC_df[rownames(FL_primary_NRBC_df),])
FL_ABM_nRBC_df_df['DLK1',]

dir.create('res_data/bulk_RNAseq')
colnames(FL_ABM_nRBC_df_df)=gsub(pattern = 'healthy_',replacement = 'ABM_',colnames(FL_ABM_nRBC_df_df))
colnames(FL_ABM_nRBC_df_df)=gsub(pattern = 'FLP_',replacement = 'FL_',colnames(FL_ABM_nRBC_df_df))

saveRDS(FL_ABM_nRBC_df_df,file = 'res_data/bulk_RNAseq/FL_ABM_nRBC_df_df_nr_exp.rds')

table(trans_CC_hDEGs_df$cluster) # fetal/adult 5/9
an_row_df=data.frame(row.names =as.character(trans_CC_hDEGs_df$gene),type=c(rep('fetal','5'),rep('adult',9)) ) 
p=pheatmap(FL_ABM_nRBC_df_df[as.character(trans_CC_hDEGs_df$gene),],annotation_row = an_row_df,cluster_cols = F,cluster_rows = F,color = colorRampPalette(colors = c('white','firebrick3'))(100))
ggsave(as.ggplot(p),filename='res_pic/main_figure5/fetal_adult_trans_markers_bulk_FL_ABM_NRBC_heatmap.pdf',width = 8,height = 6)

trans_CC_hDEGs_df$bulkRNAseq='pass'
trans_CC_hDEGs_df[trans_CC_hDEGs_df$gene %in% c('SLC6A9','ADORA2B','LRP6','PLXND1','PDZK1IP1'),'bulkRNAseq']='fail'
trans_CC_hDEGs_df[trans_CC_hDEGs_df$gene %in% c('CD34','CXCR4'),'bulkRNAseq']='filt'

#---------UCB NRBC expression level---------------------#
new_UCB_NRBC_altas=readRDS('../NRBC_UCB_altas/dealt_symbol_UCB_NRBC.rds')
prebulkRNA_df=AggregateExpression(new_UCB_NRBC_altas,group.by = 'source_celltype',features = rownames(new_UCB_NRBC_altas))$RNA
prebulkRNA_df=log2(cpm(prebulkRNA_df)+1)
DLK1_UCB_df=rbind(prebulkRNA_df['DLK1',grep('fetal',colnames(prebulkRNA_df))],
      prebulkRNA_df['DLK1',grep('pre',colnames(prebulkRNA_df))],
      prebulkRNA_df['DLK1',grep('term',colnames(prebulkRNA_df))])
colnames(DLK1_UCB_df)=gsub(pattern = 'fetal-UCB ',replacement = '',colnames(DLK1_UCB_df))
rownames(DLK1_UCB_df)=c('fetal-NRBC','pre-NRBC','term-NRBC')
p=pheatmap(DLK1_UCB_df,cluster_rows = F,cluster_cols = F,color = colorRampPalette(colors = c('white','firebrick3'))(100))
ggsave(as.ggplot(p),filename='res_pic/main_figure5/DLK1_UCB_NRBC_heatmap.pdf',width = 6,height = 3)


# K562 
K562_RNAseq_df=read.csv('..//ref_data/bulk_RNAseq/erythroblast_cellline_data/GSE311284_HRG1_K562_counts.csv')
rownames(K562_RNAseq_df)=K562_RNAseq_df[,'X']
K562_RNAseq_df=K562_RNAseq_df[,-1]
library(edgeR)
K562_RNAseq_df=log2(cpm(K562_RNAseq_df)+1)
an_df=data.frame(row.names = c(top_fetal_unique_marker_genes,ABM_cho_topmarkers),type=c(rep('fetal',length(top_fetal_unique_marker_genes)),rep('adult',length(ABM_cho_topmarkers)) ))
p=pheatmap(K562_RNAseq_df[c(top_fetal_unique_marker_genes,ABM_cho_topmarkers)[c(top_fetal_unique_marker_genes,ABM_cho_topmarkers) %in% rownames(K562_RNAseq_df)],],
         color=colorRampPalette(colors = c('navy','white','firebrick3'))(100),annotation_row =an_df ,cluster_rows = F,cluster_cols = F,main = 'K562')
ggsave(as.ggplot(p),filename='res_pic/main_figure5/K562_tp10stage_marker_expression_heatmap.pdf',width = 6,height = 9)



#-----------------------------------------------------1.2 protein level validated ---------------------------------------------#
# primary tissue data------------------#
# Cell-specific proteome analyses of human bone marrow reveal molecular features of age-dependent functional decline
BM_nRBC_MS_df=read.csv('../ref_data/Protein_NRBC/NC_2018_BM_hemo_celltype_LF_MS_data.tsv',header = T,sep="\t",skip = 1)
BM_nRBC_MS_df=BM_nRBC_MS_df[,c(colnames(BM_nRBC_MS_df)[1:4],'ERP.number.of.donors','ERP.normalized.LF.sum')]
BM_nRBC_MS_df=BM_nRBC_MS_df[!is.na(BM_nRBC_MS_df$ERP.normalized.LF.sum),]
BM_Ery_MS_protein_LF_genes=unique(unlist(strsplit(BM_nRBC_MS_df$gene.name[BM_nRBC_MS_df$ERP.normalized.LF.sum >0.1],split = ';')));length(BM_Ery_MS_protein_LF_genes)

#Pedro L. Moura Erythroid Differentiation Enhances RNA Mis-Splicing in  SF3B1-Mutant Myelodysplastic Syndromes with Ring  Sideroblasts
adult_HSPC_Ery_MS_protein_df=read.table(file = '../Protein_NRBC_marker/recent_MS_Ery_protein_omics/dealt_adult_HSPC_Ery_MS_protein_df.tsv',sep="\t")
adult_HSPC_Ery_MS_protein_genes=unique(unlist(strsplit(adult_HSPC_Ery_MS_protein_df$Gene.names,split = ';')));length(adult_HSPC_Ery_MS_protein_genes)

cho_trans_genes=c('DLK1',)

gene_an_df=data.frame(row.names = as.character(trans_CC_hDEGs_df$gene[trans_CC_hDEGs_df$bulkRNAseq=='pass']), gene=as.character(trans_CC_hDEGs_df$gene[trans_CC_hDEGs_df$bulkRNAseq=='pass']))
gene_an_df$gene=factor(gene_an_df$gene,levels = gene_an_df$gene)
gene_an_df$MS='no'
gene_an_df$MS[rownames(gene_an_df) %in% c(BM_Ery_MS_protein_LF_genes,adult_HSPC_Ery_MS_protein_genes)]='detected'
gene_an_df$HPA_IHC='un'
gene_an_df$HPA_IHC[gene_an_df$gene %in% c('ROBO2','HLA-DQB1','HLA-DRB5','CD40LG')]='Positive'
gene_an_df$HPA_IHC[gene_an_df$gene %in% c('DLK1')]='Negtive'# may low
gene_an_df$HPA_IHC[rownames(gene_an_df) %in% c('CLEC2B','CD69')]='no_sure'# may low

p1=ggplot(gene_an_df,aes(x=gene,y='MS',shape=MS))+geom_point()+theme_classic()+RotatedAxis()+scale_shape_manual(values = c(3, 1,4))+ylab('')+theme(axis.text.x = element_blank())

p2=ggplot(gene_an_df,aes(x=gene,y='HPA_IHC',shape=HPA_IHC))+geom_point()+theme_classic()+RotatedAxis()+scale_shape_manual(values = c(4, 1,3))+ylab('')


p=p1/p2;p
ggsave(p,filename='res_pic/main_figure5/protein_fetal_adult_canididated_specific_markers.pdf',width = 6,height = 6)


#-------------------------------adult F cell-----------------------------------------#
cell_name_df=read.csv('ref_data/erythrocyte_cut_cell_name.csv',sep=',')
gene_name_df=read.csv('ref_data/erythrocyte_cut_genes_name.csv',sep=',')
library(Matrix)
count_mtx=readMM(file = 'ref_data/erythrocyte_cut.mtx')
count_mtx=t(count_mtx)

rownames(count_mtx)= gene_name_df$feature_name
colnames(count_mtx)=cell_name_df$X0

sample_info_df=read.csv('ref_data/erythrocyte_info.csv')
rownames(sample_info_df)=sample_info_df$X;sample_info_df=sample_info_df[,-1]

length(unique(gene_name_df$feature_name));length(gene_name_df$feature_name) # symbol 存在重复
count_mtx_uniq=count_mtx[!duplicated(rownames(count_mtx)),]

F_cell_seu=CreateSeuratObject(counts =count_mtx_uniq,project = 'blood_F_cells')
F_cell_seu=NormalizeData(F_cell_seu) %>%FindVariableFeatures() %>% ScaleData()
F_cell_seu@meta.data=sample_info_df

fetal_markers=c('HBG1','HBG2','HBZ', "TUBB6","HSPA1A","HSPA1B",'IGF2BP1','IGF2BP3','DLK1','CISH','HIF3A')
adult_markers=c( 'HBD',"CA1","PDZK1IP1", "ANXA1","NECAB1","ANKRD28","TSC22D3","IFIT1B",'LGALS9','CXCR4')

p=VlnPlot(F_cell_seu,group.by = 'sex',features = c('HBG1','HBG2','DLK1','IGF2BP1','IGF2BP3','HBD','CA1',"PDZK1IP1","TSC22D3","IFIT1B"),stack = F,ncol = 5,layer = 'data')
ggsave(p,filename='res_pic/main_figure5/adult_F_cell_marker_vlnplot.pdf',width = 12,height = 6)




############################################################################################################################################################################
#-------------------------------------validated the DLK1 marker-------------------------#
############################################################################################################################################################################
library(tibble)
library(edgeR)

fetal_markers=c('HBG1','HBG2','HBZ', "TUBB6","HSPA1A","HSPA1B",'IGF2BP1','IGF2BP3','DLK1','CISH','HIF3A')
adult_markers=c( 'HBD',"CA1","PDZK1IP1", "ANXA1","NECAB1","ANKRD28","TSC22D3","IFIT1B",'LGALS9','CXCR4')


# UCB  wpc NRBCs

UCB_NRBCs_df=read.csv('../zx_lab_NRBC/experiment_data/UCB_41_wpc_NRBCs_df.csv')
UCB_NRBCs_df=UCB_NRBCs_df[!duplicated(UCB_NRBCs_df$gene_symbol),]
rownames(UCB_NRBCs_df)=UCB_NRBCs_df$gene_symbol
UCB_NRBCs_df=UCB_NRBCs_df[,2:3]
UCB_NRBCs_df=log2(cpm(UCB_NRBCs_df)+1)
cho_df=UCB_NRBCs_df[rownames(UCB_NRBCs_df) %in% c(fetal_markers,adult_markers),]
cho_df=cho_df[c('HBG1','HBG2','HSPA1A','HSPA1B','IGF2BP1','CA1','HBD','TSC22D3','LGALS9','IFIT1B','CXCR4'),]
p=pheatmap(cho_df,col=colorRampPalette(colors = c('navy','white','firebrick3'))(100),cluster_rows = F,cluster_cols = F)
ggsave(as.ggplot(p),filename='res_pic/main_figure5/UCB_smartseq2_fNRBC_heatmap.pdf',height = 6,width = 3)


## CD235a+DRAQ5+

smartseq_NRBC_exp_df=read.csv('../zx_lab_NRBC/experiment_data/202606_CD235aDRAQ5_NRBC.csv')
smartseq_NRBC_exp_df=data.frame(row.names = smartseq_NRBC_exp_df$X,NRBC=smartseq_NRBC_exp_df$mPBMC1)
smartseq_NRBC_exp_df=log2(cpm(smartseq_NRBC_exp_df)+1)
markers=c(fetal_markers,adult_markers)[c(fetal_markers,adult_markers) %in% rownames(smartseq_NRBC_exp_df) ]
cho_smartseq_NRBC_exp_df=data.frame(row.names = markers,smartseq_NRBC_exp_df[markers,])
colnames(cho_smartseq_NRBC_exp_df)='NRBC'   

#------------------- CD235a+DRAQ5+DLK1+----------------------------#

DLK_pos_NRBC_samrtseq2_df=read.csv('../zx_lab_NRBC/experiment_data/202607_DLK1_pos_NRBC.csv')
rownames(DLK_pos_NRBC_samrtseq2_df)=DLK_pos_NRBC_samrtseq2_df$X
DLK_pos_NRBC_samrtseq2_df=DLK_pos_NRBC_samrtseq2_df[,-1]
DLK_pos_NRBC_samrtseq2_df=log2(edgeR::cpm(DLK_pos_NRBC_samrtseq2_df)+1)# 
DLK_pos_NRBC_samrtseq2_df=DLK_pos_NRBC_samrtseq2_df[c(fetal_markers,adult_markers)[c(fetal_markers,adult_markers) %in% rownames(DLK_pos_NRBC_samrtseq2_df)],]
DLK_pos_NRBC_samrtseq2_df=data.frame(DLK_pos_NRBC_samrtseq2_df)


DLK_pos_NRBC_samrtseq2_df[ rownames(cho_smartseq_NRBC_exp_df)[!rownames(cho_smartseq_NRBC_exp_df) %in% rownames(DLK_pos_NRBC_samrtseq2_df)],c('DLK1_rep1', 'DLK1_rep2')]=0
cho_smartseq_NRBC_exp_df[ rownames(DLK_pos_NRBC_samrtseq2_df)[!rownames(DLK_pos_NRBC_samrtseq2_df) %in% rownames(cho_smartseq_NRBC_exp_df)],c('NRBC')]=0

merged_df <- merge(cho_smartseq_NRBC_exp_df,DLK_pos_NRBC_samrtseq2_df, by = "row.names")
rownames(merged_df)=merged_df$Row.names;merged_df=merged_df[,-1]
merged_df=merged_df[rowSums(merged_df) >0,]
merged_df=merged_df[c('HBG1','HBG2','HSPA1B','IGF2BP3','CA1','HBD','ANXA1','TSC22D3','CXCR4'),]

p=pheatmap(merged_df,cluster_cols = F,col=colorRampPalette(colors = c('navy','white','firebrick3'))(100),cluster_rows = F)
ggsave(as.ggplot(p),filename='res_pic/main_figure5/mPBMC_smartseq2_fNRBC_heatmap.pdf',height = 6,width = 4)



