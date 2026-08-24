#---------------------------------------------------tissue RNA marker ----------------------------------------------------#

################################################################################################################################################################
#---------------------------------------------------prepare pkg and function ----------------------------------------------------#
################################################################################################################################################################
library(ggrepel)
library(Seurat)
library(SeuratData)
library(cowplot)
library(dplyr)
library(SingleCellExperiment)
library(RColorBrewer )
library(org.Hs.eg.db)
library(AnnotationDbi)
library(tibble)
library(pheatmap)
library(clusterProfiler)
library(ggplotify)
library(ggplot2)
library(reshape2)
library(corrplot)
library(APL)
library(patchwork)
library(yulab.utils)

#install.packages("patchwork")
library(RColorBrewer )
cols=c(brewer.pal(12,"Set3"),brewer.pal(6,"PiYG"),brewer.pal(6,"BrBG"),brewer.pal(8,"Set2"),
       brewer.pal(12,"Set3"),brewer.pal(8,"Pastel2"),brewer.pal(9,"Pastel1"),brewer.pal(8,"Accent"))
col=unique(cols)[-14]
getwd()
setwd('/home/gibh/2021_NRBC_chlyu')

dir.create('Protein_NRBC_marker')

# Custom functions ------------------------------------------------------
find_mDEGs_func=function(seu,group,sfile)
  early_Ery_tissue_markers=FindAllMarkers(subset(seu,Ery_stage=='early_Ery'))
  mid_Ery_tissue_markers=  FindAllMarkers(subset(seu,Ery_stage=='mid_Ery'))
  late_Ery_tissue_markers= FindAllMarkers(subset(seu,Ery_stage=='late_Ery'))
  
  early_Ery_tissue_markers$celltype='early_Ery'
  mid_Ery_tissue_markers$celltype='mid_Ery'
  late_Ery_tissue_markers$celltype='late_Ery'
  all_Ery_tissue_markers=rbind(early_Ery_tissue_markers,rbind(mid_Ery_tissue_markers,late_Ery_tissue_markers))
  rownames(all_Ery_tissue_markers)=NULL
  all_Ery_tissue_markers$group=group
   
  write.csv(all_Ery_tissue_markers,file = sfile,quote = F)
  rm(early_Ery_tissue_markers,mid_Ery_tissue_markers,late_Ery_tissue_markers);gc()
  
  pos_all_Ery_tissue_markers=all_Ery_tissue_markers[all_Ery_tissue_markers$avg_log2FC >0  & all_Ery_tissue_markers$p_val_adj <0.01,]
  if(dim(pos_all_Ery_tissue_markers)[2] <1 ){ stop('no DEGs are found!!!')}
  
  count_tmp_df=cbind( data.frame(table(pos_all_Ery_tissue_markers[,c('cluster','celltype')])[1,]), data.frame(table(pos_all_Ery_tissue_markers[,c('cluster','celltype')])[2,]))
  colnames(count_tmp_df)= strsplit(group,'_')[[1]]
  type1=strsplit(group,split = '_')[[1]][1]
  type2=strsplit(group,split = '_')[[1]][2]
  
  count_tmp_df['all',type1]=length(unique(pos_all_Ery_tissue_markers$gene[pos_all_Ery_tissue_markers$cluster==type1]))
  count_tmp_df['all',type2]=length(unique(pos_all_Ery_tissue_markers$gene[pos_all_Ery_tissue_markers$cluster==type2]))
  count_tmp_df['add_all',]=length(unique(pos_all_Ery_tissue_markers$gene))


  return(list(all_Ery_tissue_markers,count_tmp_df))
}


subcelltype_gseGO_func=function(RNA_markers,keyType='SYMBOL'){
  subcelltype_gseGO_glist=list()
  subcelltype_gseGO_list=list()
  tmp_df=data.frame()
  for (type in unique(RNA_markers$celltype)) {
    print(type)
    subcelltype_gseGO_glist[[type]]=RNA_markers[RNA_markers$celltype==type,'avg_log2FC']
    names(subcelltype_gseGO_glist[[type]])=RNA_markers[RNA_markers$celltype==type,'gene']
    subcelltype_gseGO_glist[[type]]=sort(subcelltype_gseGO_glist[[type]],decreasing = T)
    subcelltype_gseGO_list[[type]] =gseGO(geneList = subcelltype_gseGO_glist[[type]],ont ='all' ,OrgDb = org.Hs.eg.db,keyType = keyType)
    
    if(dim(subcelltype_gseGO_list[[type]]@result)[1] >0){
      subcelltype_gseGO_list[[type]]@result$celltype=type
      tmp_df=rbind(tmp_df,subcelltype_gseGO_list[[type]]@result)}
  }
  return(list(subcelltype_gseGO_list,tmp_df))
}
subcelltype_enrichGO_func=function(RNA_markers){
  subcelltype_enrichGO_glist=list()
  subcelltype_enrichGO_list=list()
  RNA_markers=RNA_markers[RNA_markers$avg_log2FC>0,]
  for (type in unique(RNA_markers$celltype)) {
    cluster=as.character(unique(RNA_markers$cluster))
    subcelltype_enrichGO_glist[[type]][[ cluster[1] ]]=unique(RNA_markers[RNA_markers$cluster==cluster[1] & RNA_markers$celltype==type ,'gene'])
    subcelltype_enrichGO_glist[[type]][[ cluster[2] ]]=unique(RNA_markers[RNA_markers$cluster==cluster[2] & RNA_markers$celltype==type ,'gene'])
    subcelltype_enrichGO_list[[type]][[ cluster[1] ]] =enrichGO(gene =subcelltype_enrichGO_glist[[type]][[ cluster[1] ]] ,OrgDb =org.Hs.eg.db ,keyType = 'SYMBOL',ont = 'all')
    subcelltype_enrichGO_list[[type]][[ cluster[2] ]] =enrichGO(gene =subcelltype_enrichGO_glist[[type]][[ cluster[2] ]] ,OrgDb = org.Hs.eg.db,keyType = 'SYMBOL',ont = 'all')
  }
  
  return(subcelltype_enrichGO_list)
}




#################################################################################################################################################################
#--------------------------------------------------prepare data ----------------------------------------------------#
#################################################################################################################################################################

#--------------------------prepare NRBC data-----------------------#
if(F){
  YS_altas_Ery_seu=readRDS('NRBC_YS_altas/YS_altas_Ery_seu.rds')
  FL_altas_Ery_seu=readRDS('NRBC_FL_altas/tmp_FL_altas_Ery_seu.rds')
  BM_NRBC_altas_seu=readRDS('NRBC_BM_altas/res_data/BM_NRBC_altas_seu.rds')
  
  (DimPlot(YS_altas_Ery_seu,group.by = 'final_celltype',reduction = 'ref.umap',cols = cols)+ggtitle('YS'))+
    (DimPlot(FL_altas_Ery_seu,group.by = 'subcelltype',reduction = 'umap2',cols = cols)+ggtitle('FL'))
  DimPlot(BM_NRBC_altas_seu,group.by = 'new_celltype',reduction = 'ref.umap',cols = cols,split.by = 'stage')
  
  # early Ery: BFUE/CFUE, ProE
  # mid Ery: Bas
  # late Ery: Poly, Orth
  
  # FL HBE+HBB+
  HBB_HBE1_Ery_seu=subset(FL_altas_Ery_seu,HBE1 >1 & HBB >1);HBB_HBE1_Ery_seu# 1924, FL 存在HBE1—>HBB转变的过程, 以后可以从这个数据中研究转变过程
  subset(FL_altas_Ery_seu,HBE1 >0 & HBB <1) # 2691
  table(HBB_HBE1_Ery_seu$id)# 主要来自CS18_F61样本，该样本也是YS——Ery的主要来源
  # 由高转低，在小鼠中发现FL 存在YS-EMP来源的Ery，但是人里面是缺乏的，我们之前的FL数据显示。这是跟小鼠的不同点，YS-EMP主要产生了mon-derived Max。
  table( subset(FL_altas_Ery_seu,HBE1 >0 & HBB <1)$id)/table(FL_altas_Ery_seu$id)
  
  head(BM_NRBC_altas_seu)
 
  
  NBRC_altas_seu=merge(YS_altas_Ery_seu,c(FL_altas_Ery_seu,BM_NRBC_altas_seu));
  NBRC_altas_seu[['prediction.score.celltype']]=NULL
  NBRC_altas_seu[['RNA']]=JoinLayers(NBRC_altas_seu[['RNA']])
  NBRC_altas_seu[['umap']]=merge(YS_altas_Ery_seu[['ref.umap']],c(FL_altas_Ery_seu[['umap2']],BM_NRBC_altas_seu[['ref.umap']]))
  
  NBRC_altas_seu$tissue_stage=NBRC_altas_seu$stage
  NBRC_altas_seu@meta.data[rownames(YS_altas_Ery_seu@meta.data),c('tissue_stage')]='YS'
  NBRC_altas_seu@meta.data[rownames(FL_altas_Ery_seu@meta.data),c('tissue_stage')]='FL'
  
  NBRC_altas_seu$final_celltype=as.character(NBRC_altas_seu$new_celltype)
  NBRC_altas_seu@meta.data[rownames(YS_altas_Ery_seu@meta.data),c('final_celltype')]=as.character(YS_altas_Ery_seu$final_celltype)
  NBRC_altas_seu@meta.data[rownames(FL_altas_Ery_seu@meta.data),c('final_celltype')]=as.character(FL_altas_Ery_seu$subcelltype)
  table(NBRC_altas_seu$final_celltype)
  
  # FILT OUT YS_NRBC
  filt_cells=rownames(NBRC_altas_seu@meta.data[NBRC_altas_seu$tissue_stage %in% c('FL','FBM') & NBRC_altas_seu$final_celltype %in% c("YS_Bas/Poly","YS_Orth"),])
  cho_cells=rownames(NBRC_altas_seu@meta.data)[!rownames(NBRC_altas_seu@meta.data) %in%  filt_cells ];length(cho_cells)
  filt_NBRC_altas_seu=subset(NBRC_altas_seu,cells=cho_cells)
  filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu, tissue_stage!='EBM')
  filt_NBRC_altas_seu$tissue_stage=factor(filt_NBRC_altas_seu$tissue_stage,levels = c('YS','FL','FBM','ABM'))
  filt_NBRC_altas_seu$final_celltype=factor(filt_NBRC_altas_seu$final_celltype,levels = c("BFUE/CFUE","ProE","Bas","Poly" ,"Orth"))
  DimPlot(filt_NBRC_altas_seu,group.by = 'final_celltype',split.by = 'tissue_stage',cols = cols,raster=FALSE)/
    FeaturePlot(filt_NBRC_altas_seu,features = 'HBE1',split.by = 'tissue_stage')
  rm(NBRC_altas_seu);gc()
  
  
  filt_NBRC_altas_seu$Ery_stage='late_Ery'
  filt_NBRC_altas_seu$Ery_stage[filt_NBRC_altas_seu$final_celltype %in% c('BFUE/CFUE','ProE')]='early_Ery'
  filt_NBRC_altas_seu$Ery_stage[filt_NBRC_altas_seu$final_celltype %in% c('Bas')]='mid_Ery'
  filt_NBRC_altas_seu$Ery_stage=factor(filt_NBRC_altas_seu$Ery_stage,levels = c('early_Ery','mid_Ery','late_Ery'))
  filt_NBRC_altas_seu$tissue_stage=factor(filt_NBRC_altas_seu$tissue_stage,levels = c('YS','FL','FBM','ABM'))
  
  table(filt_NBRC_altas_seu@meta.data[,c('Ery_stage','final_celltype')])
  
  rm(YS_altas_Ery_seu,FL_altas_Ery_seu,FBM_NRBC_altas_seu,BM_NRBC_altas_seu);gc()
  
  Idents(filt_NBRC_altas_seu)='tissue_stage'
  
  
  # FL/FBM FILT OUT HBE1+ high
  temp_seu= subset(filt_NBRC_altas_seu, tissue_stage %in% c('FL','FBM'))
  filt_cells2=WhichCells(temp_seu,expression =  HBE1 >2,slot = 'data')
  length(filt_cells2) # HBE1 >1 : 2245, HBE1 >2: 841 
  DimPlot(subset(filt_NBRC_altas_seu,cells=filt_cells2),group.by = 'Ery_stage',split.by = 'tissue_stage',cols = cols,raster=FALSE)+
    FeaturePlot(subset(filt_NBRC_altas_seu,cells=filt_cells2),features = 'HBE1',split.by = 'tissue_stage')
  cho_cells2=rownames(filt_NBRC_altas_seu@meta.data)[!rownames(filt_NBRC_altas_seu@meta.data) %in%  filt_cells2 ];length(cho_cells2)
  filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu,cells=cho_cells2)
  rm(temp_seu,cho_cells,cho_cells2,filt_cells,filt_cells2);gc()
  
  sort(unique(filt_NBRC_altas_seu$final_celltype))
  filt_NBRC_altas_seu$age[is.na(filt_NBRC_altas_seu$age)]=filt_NBRC_altas_seu$stage[is.na(filt_NBRC_altas_seu$age)]
  
  colnames(filt_NBRC_altas_seu@meta.data)=gsub(pattern = 'resource',replacement ='source' ,colnames(filt_NBRC_altas_seu@meta.data))
  filt_NBRC_altas_seu$source[is.na(filt_NBRC_altas_seu$source)]=filt_NBRC_altas_seu$orig.dataset[is.na(filt_NBRC_altas_seu$source)]
  filt_NBRC_altas_seu$age[filt_NBRC_altas_seu$tissue_stage=='YS']=filt_NBRC_altas_seu$stage[filt_NBRC_altas_seu$tissue_stage=='YS']
  
  filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu,source !='GSE253355')
  
  filt_NBRC_altas_seu_meta=filt_NBRC_altas_seu@meta.data
  filt_NBRC_altas_seu_meta=filt_NBRC_altas_seu_meta[,c('orig.ident', 'nCount_RNA', 'nFeature_RNA', 'component', 'age',  'sex', 'sort.ids', 'fetal.ids','orig.dataset',
                                                       'sequencing.type','id','tissue_stage', 'final_celltype', 'Ery_stage','source', 'source_celltype')]
  
  
  
  VlnPlot(filt_NBRC_altas_seu,features =as.character(cor_HIST_genes),stack = T )+NoLegend()
  VlnPlot(filt_NBRC_altas_seu,features =names(cor_HIST_genes)[!names(cor_HIST_genes) %in% c('HIST1H2BA','HIST1H4G')],stack = T )+NoLegend()
  
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  filt_NBRC_altas_seu <- CellCycleScoring(filt_NBRC_altas_seu, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
  filt_NBRC_altas_seu$source_celltype=paste(filt_NBRC_altas_seu$tissue_stage,filt_NBRC_altas_seu$final_celltype,sep = '_')
  filt_NBRC_altas_seu$source_celltype=factor(filt_NBRC_altas_seu$source_celltype,levels = c( "YS_ProE","YS_Bas","YS_Poly","YS_Orth", "FL_BFUE/CFUE" ,"FL_ProE", "FL_Bas","FL_Poly","FL_Orth",
                                                                                             "FBM_BFUE/CFUE" ,"FBM_ProE","FBM_Bas","FBM_Poly","FBM_Orth", "ABM_BFUE/CFUE","ABM_ProE","ABM_Bas","ABM_Poly","ABM_Orth"))
  filt_NBRC_altas_seu$final_celltype=factor(filt_NBRC_altas_seu$final_celltype,levels = c('BFUE/CFUE','ProE','Bas','Poly','Orth'))
  filt_NBRC_altas_seu$sourceid=paste(filt_NBRC_altas_seu$tissue_stage,filt_NBRC_altas_seu$id,sep = '_')
  
  filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu,subset=IGKC <=1 ) # FILT OUT IMMUNE B CELLS
  IGKC_nRBC_seu==subset(filt_NBRC_altas_seu,subset=IGKC <=1 )
  saveRDS(filt_NBRC_altas_seu,file = '20251125_filt_NBRC_altas_seu.rds')
  
}else{
  filt_NBRC_altas_seu=readRDS('20251125_filt_NBRC_altas_seu.rds')
}



#################################################################################################################################################################
#--------------------------------------------------figure 1: hematopoietic niche shape the differenation of nRBC---------------------------------------------------#
#################################################################################################################################################################
filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu,tissue_stage!='YS')
filt_NBRC_altas_seu=FindVariableFeatures(filt_NBRC_altas_seu)

setwd('definitive_nRBC_marker_project')
dir.create('res_pic/main_figure1')
dir.create('res_data')

p1=pheatmap(cor(data.frame(AverageExpression(filt_NBRC_altas_seu,group.by = 'source_celltype',layer = 'data',features =rownames(filt_NBRC_altas_seu) )$RNA)),
           color =colorRampPalette(colors = c('white','#C31E1F'))(100),main = 'based on all genes')
ggsave(as.ggplot(p1),filename = 'res_pic/main_figure1/tissue_basedonall_cor_heatmap.pdf',width = 6,height = 6,dpi = 300)

p2=pheatmap(cor(data.frame(AverageExpression(filt_NBRC_altas_seu,group.by = 'source_celltype',layer = 'data',features =VariableFeatures(filt_NBRC_altas_seu) )$RNA)),
           cluster_cols =F ,color =colorRampPalette(colors = c('white','#C31E1F'))(100),main = 'based on variable genes')
ggsave(as.ggplot(p2),filename = 'res_pic/main_figure1/tissue_basedonvariable_cor_heatmap.pdf',width = 6,height = 6,dpi = 300)


#   subcelltype in steady erythropoeisis  across niches
filt_NBRC_altas_seu_meta=filt_NBRC_altas_seu@meta.data
FL_CS18post_sample_df=filt_NBRC_altas_seu_meta[filt_NBRC_altas_seu_meta$tissue_stage=='FL' & !filt_NBRC_altas_seu_meta$age %in% c('CS14_4PCW', "CS15_5PCW" ,"CS17_6PCW" ),]
FL_celltype_mratio_df=data.frame( colMeans(prop.table(table(FL_CS18post_sample_df[,c('id','final_celltype')]),margin = 1)))
colnames(FL_celltype_mratio_df)='FL'

FBM_sample_df=filt_NBRC_altas_seu_meta[filt_NBRC_altas_seu_meta$tissue_stage=='FBM',]
FBM_celltype_mratio_df=data.frame( colMeans(prop.table(table(FBM_sample_df[,c('id','final_celltype')]),margin = 1)))
colnames(FBM_celltype_mratio_df)='FBM'

ABM_sample_df=filt_NBRC_altas_seu_meta[filt_NBRC_altas_seu_meta$tissue_stage=='ABM',]
table(ABM_sample_df$source)# CD34- CD235+: GSE150774:CD34+,GSE133181,GSE135194,GSE169426,  Mononuclear cells:GSE165645, Mononuclear:CD34+=4:1:GSE181989
ABM_celltype_mratio_CD34pos_df=data.frame( round(colMeans(table(ABM_sample_df[!ABM_sample_df$source %in% c('GSE150774','GSE181989','GSE165645'),c('id','final_celltype')])),digits = 0) )
colnames(ABM_celltype_mratio_CD34pos_df)='mcount'
ABM_celltype_mratio_GYPAposCD34neg_df=data.frame( round(colMeans(table(ABM_sample_df[ABM_sample_df$source=='GSE150774',c('id','final_celltype')])),digits = 0))
colnames(ABM_celltype_mratio_GYPAposCD34neg_df)='mcount'
ABM_celltype_mratio_GYPAposCD34neg_df['BFUE/CFUE','mcount']=0

ABM_NRBC_ratio_df=prop.table(ABM_celltype_mratio_CD34pos_df+ ABM_celltype_mratio_GYPAposCD34neg_df[rownames(ABM_celltype_mratio_CD34pos_df),])
colnames(ABM_NRBC_ratio_df)='ABM'

NRBC_ratio_df= cbind(cbind(FL_celltype_mratio_df,FBM_celltype_mratio_df),ABM_NRBC_ratio_df)
NRBC_ratio_df$celltytpe=rownames(NRBC_ratio_df)
NRBC_ratio_df= melt(NRBC_ratio_df)

NRBC_ratio_df$celltytpe=factor(NRBC_ratio_df$celltytpe,levels = c("BFUE/CFUE" ,"ProE","Bas","Poly","Orth"  ))
p=ggplot(NRBC_ratio_df ,aes(x=variable,y=value,fill=celltytpe))+geom_bar(position = "fill",stat = 'identity')+
  scale_fill_manual(values =  cols)+theme_classic()+theme(axis.text.x = element_text(hjust =1 ,angle = 45,face = 'bold'))
p
ggsave(p,file='res_pic/main_figure1/tissue_celltype_ratio_barplot.pdf',width =4 ,height = 6)


# NRBC reference marker
p=VlnPlot(filt_NBRC_altas_seu,group.by = 'source_celltype',features =c('KIT','KLF1','TFRC','GYPA','CCNB1','NCL'),stack = T,cols = col,split.by = 'tissue_stage')+NoLegend();p
ggsave(p,filename='res_pic/main_figure1/classical_nRBC_marker_expression_in_niches_subceltlype_vlnplot.pdf',width = 8,height = 8)



#################################################################################################################################################################
#------------------------------------------figure 2 ：the comparision of HSPC_derived_NRBC -------------------------------------------#
#################################################################################################################################################################

Idents(filt_NBRC_altas_seu)='tissue_stage'
HSPC_derived_NRBC_DE_res=FindAllMarkers(subset(filt_NBRC_altas_seu,NRBC_type=='definitive'))
saveRDS(HSPC_derived_NRBC_DE_res,file ='Protein_NRBC_marker/DE_marker/HSPC_derived_nRBC_wholelevel_RNA_markers.rds' )

# expression of top markers

top_gene_markers_HSPC_derived_NRBC=HSPC_derived_NRBC_DE_res[-grep('^AC[0-9]|CH507-|^LIN',HSPC_derived_NRBC_DE_res$gene),] %>% filter(avg_log2FC>1 & pct.2 < 0.2 & pct.1 > 0.2) %>% group_by(cluster)  %>%top_n(wt=avg_log2FC,10)# %>%  do(head(., n = 10))
top_gene_markers_HSPC_derived_NRBC=top_gene_markers_HSPC_derived_NRBC[order(top_gene_markers_HSPC_derived_NRBC$cluster,top_gene_markers_HSPC_derived_NRBC$avg_log2FC,decreasing = T),]
p=DotPlot(subset(filt_NBRC_altas_seu,NRBC_type=='definitive'),group.by = 'source_celltype',features =unique(top_gene_markers_HSPC_derived_NRBC$gene)[c(21:30,11:20,1:10)],cols = c('gray','firebrick3'),scale = F)+RotatedAxis() 
p

FL_top10_markers=unique(top_gene_markers_HSPC_derived_NRBC$gene)[21:30]

#Top 10 marker analysis revealed that FBM-derived NRBCs expressed significantly higher levels of replication-dependent histones, which contradicts known biological expectations.
#A likely explanation is that the 10x Genomics platform relies on poly(A)-tail capture of mRNA, whereas replication-dependent histone genes lack poly(A) tails and 
# instead use a 3′ terminal stem–loop structure bound by stem–loop binding protein (SLBP).

#Histone genes were excluded from the marker analysis.
HIS_genes=mapIds(x = org.Hs.eg.db,keys = rownames(
filt_NBRC_altas_seu)[grep('^HIS',rownames(filt_NBRC_altas_seu))],keytype = 'ALIAS',column = 'SYMBOL') 
HIS_genes=HIS_genes[as.character(HIS_genes)!=names(HIS_genes)] # 全部为复制依赖行组蛋白，仅仅在FBM中检测高，
HSPC_derived_NRBC_DE_res=HSPC_derived_NRBC_DE_res[!HSPC_derived_NRBC_DE_res$gene %in% names(HIS_genes),]

table(HSPC_derived_NRBC_DE_res[HSPC_derived_NRBC_DE_res$avg_log2FC >0,'cluster'])
#   FL   FBM   ABM 
# 10593   969  2186

# 查看top20
top_gene_markers_HSPC_derived_NRBC=HSPC_derived_NRBC_DE_res[-grep('^AC[0-9]|CH507-|^LIN',HSPC_derived_NRBC_DE_res$gene),] %>% filter(avg_log2FC>1 & pct.2 < 0.2 & pct.1 > 0.2) %>% group_by(cluster)  %>%top_n(wt=avg_log2FC,20)
p=DotPlot(filt_NBRC_altas_seu,group.by = 'source_celltype',features =top_gene_markers_HSPC_derived_NRBC$gene[top_gene_markers_HSPC_derived_NRBC$cluster=='FBM'],cols = c('gray','firebrick3'),scale = F)+RotatedAxis() #   colorRampPalette(colors = c('gray','firebrick3'))(100)
p
# The FBM absent of unique marker, FBM markers were still largely non-specific, with most genes involved in translation and RNA processing. filtout 
VlnPlot(UCB_NRBC_altas,group.by = 'celltype',features = top_gene_markers_HSPC_derived_NRBC$gene[top_gene_markers_HSPC_derived_NRBC$cluster=='FBM'],stack = T,split.by = 'type')
del_gene=c('EIF3C','EIF3CL','U2AF1','ATP6V0C','DDTL','BOLA2','BOLA2B','U2AF1L5','GET4','SMN1','SMN2','SLX1A','SLX1B', "SBF2-AS1",'PHOSPHO1')# 
FBM_unique_genes=top_gene_markers_HSPC_derived_NRBC[!top_gene_markers_HSPC_derived_NRBC$gene %in% del_gene & top_gene_markers_HSPC_derived_NRBC$cluster=='FBM', ]$gene

# FBM NRBCs similarity to FL NRBCs， absent unique highly expressive genes
if(F){
  fl_fbm_pos_marker=FindMarkers(filt_NBRC_altas_seu,ident.1 = 'FBM',ident.2 = 'FL')
  fbm_pos_marker=fl_fbm_pos_marker[fl_fbm_pos_marker$avg_log2FC >1 & fl_fbm_pos_marker$pct.1>0.1 & fl_fbm_pos_marker$pct.2 <0.3,]
  fbm_pos_marker=fbm_pos_marker[!rownames(fbm_pos_marker) %in% c(names(HIS_genes),del_gene),]
  fbm_pos_marker=fbm_pos_marker[-grep('^CH507|AC00',rownames(fbm_pos_marker)),]
  dim(fbm_pos_marker)# 69
  DotPlot(filt_NBRC_altas_seu,group.by = 'source_celltype',features =rownames(fbm_pos_marker),cols = c('gray','firebrick3'),scale = F)+RotatedAxis() #   colorRampPalette(colors = c('gray','firebrick3'))(100)
  # FBM NRBC缺乏特异
  
  
  FL_pos_marker=fl_fbm_pos_marker[abs(fl_fbm_pos_marker$avg_log2FC )>1 & fl_fbm_pos_marker$pct.2>0.2 & fl_fbm_pos_marker$pct.1 <0.3,]
  FL_pos_marker=FL_pos_marker[!rownames(FL_pos_marker) %in% c(names(HIS_genes),del_gene),]
  FL_pos_marker=FL_pos_marker[-grep('^CH507|AC00',rownames(FL_pos_marker)),]
  dim(FL_pos_marker) # 1563
  FL_pos_marker=FL_pos_marker[order(FL_pos_marker$avg_log2FC),]
  DotPlot(filt_NBRC_altas_seu,group.by = 'source_celltype',features =rownames(FL_pos_marker)[1:60],cols = c('gray','firebrick3'),scale = F)+RotatedAxis() #   colorRampPalette(colors = c('gray','firebrick3'))(100)
  # FL NRBC缺乏特异高表达基因, 差异表达
  
  # substage comparition
  group='FL_FBM'
  sfile='Protein_NRBC_marker/DE_marker/FL_FBM_all_Ery_RNA_markers.csv'
  res=find_DEGs_func(seu = subset(filt_NBRC_altas_seu,tissue_stage %in% c('FL','FBM')),group = group,sfile = sfile)
  FL_FBM_all_Ery_tissue_markers=res[[1]]
  FL_FBM_count_df=res[[2]]
  rm(res);gc()
  # FL 与FBM吧nRBC 更为相似可能是导致FL nRBC 相较于BM nRBC无显著差异基因的原因 
  
  # FBM中组蛋白表达太高，影响富集结果 
  FL_FBM_all_Ery_tissue_markers=FL_FBM_all_Ery_tissue_markers[!FL_FBM_all_Ery_tissue_markers$gene %in% names(HIS_genes),] 
  table(FL_FBM_all_Ery_tissue_markers[FL_FBM_all_Ery_tissue_markers$avg_log2FC >1 & FL_FBM_all_Ery_tissue_markers$pct.1 >0.1 & FL_FBM_all_Ery_tissue_markers$pct.2 <0.3,c('cluster','celltype')])
  
  
  FL_FBM_subcelltype_gseGO_list=subcelltype_gseGO_func(RNA_markers =FL_FBM_all_Ery_tissue_markers[FL_FBM_all_Ery_tissue_markers$cluster=='FL',],keyType = 'ALIAS' )
  FL_FBM_subcelltype_gseGO_list[[2]]$group='FL_FBM'
  degs_gseGO_res_df2=data.frame(FL_FBM_subcelltype_gseGO_list[[2]])
  write.csv(degs_gseGO_res_df2,file = 'Protein_NRBC_marker/res_data/main_figure2/HSPC_derived_nRBC_degs_gseGO_res_df2_ALIAS.csv')
  
  degs_gseGO_res_df2=degs_gseGO_res_df2[degs_gseGO_res_df2$ONTOLOGY=='BP',]   
  degs_gseGO_res_df2$res='up'
  degs_gseGO_res_df2$res[degs_gseGO_res_df2$NES <0]='down'
  
  top_degs_gseGO_bp_res_df=degs_gseGO_res_df2[degs_gseGO_res_df2$p.adjust <0.01,] %>% group_by(celltype,res) %>%top_n(wt = -log10(p.adjust),n=10)  %>% do(head(.,10))
  head(sort(table(top_degs_gseGO_bp_res_df$Description),decreasing = T),40)
  top_degs_gseGO_bp_res_df$Description=factor(top_degs_gseGO_bp_res_df$Description,levels = unique(top_degs_gseGO_bp_res_df$Description))
  top_degs_gseGO_bp_res_df$celltype=factor(top_degs_gseGO_bp_res_df$celltype,levels = c('early_Ery','mid_Ery','late_Ery'))
  
  
  temp_df=top_degs_gseGO_bp_res_df
  temp_df=temp_df[order(temp_df$celltype,temp_df$res,temp_df$Description),]
  temp_df$Description=factor(temp_df$Description,levels = unique(temp_df$Description))
  p2=ggplot(temp_df,aes(x=celltype,y=Description,color=NES,size=-log10( p.adjust)))+geom_point()+theme_bw()+scale_color_gradient2(low = 'blue',mid = 'white',high = 'firebrick3')+
    theme(axis.text.x = element_text(angle = 45,hjust = 1))+ggtitle(label = 'FL vs FBM nRBC gseGO of DEGs')+theme(text = element_text(face = 'bold'),panel.grid = element_blank())
  p2
  ggsave(p,width =8 ,height =12,filename='../Protein_NRBC_marker/res_pic/main_figure2/DEGS_gseGOBP_substate_FL_FBM_dotplot.pdf' )
  
  temp_df=top_degs_gseGO_bp_res_df
  temp_df=temp_df[order(temp_df$celltype,temp_df$res,temp_df$Description),]
  temp_df$Description=factor(temp_df$Description,levels = unique(temp_df$Description))
  filt_NBRC_altas_seu=AddModuleScore_UCell(filt_NBRC_altas_seu,features = list(angiogenesis=unlist(strsplit(temp_df$core_enrichment[temp_df$Description=='angiogenesis'],'/')),
                                                                               'blood_vessel_morphogenesis'=unlist(strsplit(temp_df$core_enrichment[temp_df$Description=='blood vessel morphogenesis'],'/')),
                                                                               'endothelial_cell_migration'=unlist(strsplit(temp_df$core_enrichment[temp_df$Description=='endothelial cell migration'],'/')),
                                                                               'coagulation'=unlist(strsplit(temp_df$core_enrichment[temp_df$Description=='coagulation'],'/')),
                                                                               ncores = 6))
  
  
  
  saveRDS(filt_NBRC_altas_seu@meta.data[,c('angiogenesis_UCell','blood_vessel_morphogenesis_UCell','endothelial_cell_migration_UCell','coagulation_UCell')],file = 'NRBC_angiogenesis_pathway.rds')
  
  p=VlnPlot(subset(filt_NBRC_altas_seu,tissue_stage %in% c('FL','FBM','ABM')),group.by = 'source_celltype',stack = T,cols = cols,
            features = c('angiogenesis_UCell','blood_vessel_morphogenesis_UCell','endothelial_cell_migration_UCell','coagulation_UCell', paste0(names(gene_sets),'_UCell')))
  
  ggsave(p,filename = 'Protein_NRBC_marker/res_pic/main_figure3/key_pathway_UCell_score_definitive_nRBC_vlnplot.pdf',width = 10,height = 10)
  

}

# FL & FBM NRBCs merged as fetal NRBCs, obtain top markers
candidated_ABM_specific_genes1  =HSPC_derived_NRBC_DE_res[HSPC_derived_NRBC_DE_res$cluster=='ABM' &HSPC_derived_NRBC_DE_res$avg_log2FC >2 & HSPC_derived_NRBC_DE_res$pct.1 >0.1 & HSPC_derived_NRBC_DE_res$pct.2<0.2 ,]
candidated_ABM_specific_genes1=candidated_ABM_specific_genes1[order(candidated_ABM_specific_genes1$avg_log2FC,decreasing = T),];length(candidated_ABM_specific_genes1$gene)
VlnPlot(filt_NBRC_altas_seu,group.by = 'source_celltype',features =candidated_ABM_specific_genes1$gene[1:20],stack = T,split.by = 'tissue_stage')
#"AC005943.2"    "RP11-411B6.6"     "RP11-111K18.1" , "LINC00570"， ， "RP11-354E11.2"非编码RNA特异高表达，暂时不考虑这些基因,KIAA0125 : USP45,泛素化酶,DNA损伤修复，这是USP45最明确的功能之一。
ABM_cho_topmarkers=c( 'HBD',"CA1","PDZK1IP1", "ANXA1","NECAB1","ANKRD28","TSC22D3","IFIT1B",'LGALS9',"HLA-B","HLA-DRA","HLA-DRB1") # top10 中挑选基因,
VlnPlot(filt_NBRC_altas_seu,group.by = 'source_celltype',features =ABM_cho_topmarkers,stack = T,split.by = 'tissue_stage')

whole_celltype_mexp=AverageExpression(filt_NBRC_altas_seu,group.by = 'tissue_stage')$RNA
fetal_unique_marker_genes=HSPC_derived_NRBC_DE_res[HSPC_derived_NRBC_DE_res$cluster=='ABM' & HSPC_derived_NRBC_DE_res$avg_log2FC <0,]
fetal_unique_marker_genes=fetal_unique_marker_genes[fetal_unique_marker_genes$avg_log2FC < -2 & fetal_unique_marker_genes$pct.1<0.3 & fetal_unique_marker_genes$pct.2 >0.2, ]
fetal_unique_marker_genes$score=-1*fetal_unique_marker_genes$avg_log2FC *fetal_unique_marker_genes$pct.2/(fetal_unique_marker_genes$pct.1+0.001)*whole_celltype_mexp[fetal_unique_marker_genes$gene,'FL']/(whole_celltype_mexp[fetal_unique_marker_genes$gene,'ABM']+0.001)*(fetal_unique_marker_genes$pct.2-fetal_unique_marker_genes$pct.1)
fetal_unique_marker_genes$score=log2(fetal_unique_marker_genes$score+1)
fetal_unique_marker_genes=fetal_unique_marker_genes[order(fetal_unique_marker_genes$score,-1*fetal_unique_marker_genes$avg_log2FC,decreasing = T),]
fetal_unique_marker_genes=fetal_unique_marker_genes[!fetal_unique_marker_genes$gene %in% top_gene_markers_HSPC_derived_NRBC$gene , ]
dim(fetal_unique_marker_genes)# 133
saveRDS(fetal_unique_marker_genes,file = 'Protein_NRBC_marker/res_data/main_figure2/fetal_unique_marker_genes.rds')

DotPlot(filt_NBRC_altas_seu,group.by = 'source_celltype',features =unique(fetal_unique_marker_genes$gene)[1:40],cols = c('gray','firebrick3'),scale = F)+RotatedAxis() #   colorRampPalette(colors = c('gray','firebrick3'))(100)
VlnPlot(filt_NBRC_altas_seu,group.by = 'source_celltype',features =unique(fetal_unique_marker_genes$gene)[1:20],stack = T)
temp_df=fetal_unique_marker_genes[1:20,]
temp_df=temp_df[order(temp_df$pct.2,decreasing = T),]
temp_df

VlnPlot(filt_NBRC_altas_seu,group.by = 'source_celltype',features =unique(temp_df$gene)[1:20],stack = T)

top_fetal_unique_marker_genes=c('HBG1','HBG2','HBZ', "TUBB6","HSPA1A","HSPA1B",'IGF2BP1','IGF2BP3','DLK1','CISH','HIF3A')# 还可以考虑CHD7,TIMP3,CISH,后面好像有获得

p=DotPlot(subset(filt_NBRC_altas_seu,NRBC_type=='definitive'),group.by = 'source_celltype',features =c(FL_top10_markers,FBM_unique_genes,
                                                                                                       top_fetal_unique_marker_genes,ABM_cho_topmarkers),cols = c('gray','firebrick3'),scale = F)+RotatedAxis() #   colorRampPalette(colors = c('gray','firebrick3'))(100)
p
ggsave(p,filename='Protein_NRBC_marker/res_pic/main_figure2/top_marker_HSPC_derived_nRBC.pdf',width = 16,height = 8)


HSPC_nRBC_enrichgo_res=compareCluster(geneClusters =list('FL'=unique(HSPC_derived_NRBC_DE_res[HSPC_derived_NRBC_DE_res$cluster=='FL' & HSPC_derived_NRBC_DE_res$avg_log2FC >0,'gene']),
                                                         'FBM'=unique(HSPC_derived_NRBC_DE_res[HSPC_derived_NRBC_DE_res$cluster=='FBM'& HSPC_derived_NRBC_DE_res$avg_log2FC >0,'gene']),
                                                         'ABM'=unique(HSPC_derived_NRBC_DE_res[HSPC_derived_NRBC_DE_res$cluster=='ABM'& HSPC_derived_NRBC_DE_res$avg_log2FC >0,'gene']) ),keyType = 'ALIAS',ont = "BP",
                                      fun = 'enrichGO',  OrgDb='org.Hs.eg.db' )
saveRDS(HSPC_nRBC_enrichgo_res,file = 'Protein_NRBC_marker/res_data/main_figure2/HSPC_nRBC_enrichgo_res.rds')

p=dotplot(HSPC_nRBC_enrichgo_res,showCategory=10);p
# ggplot
top10_enrichgo_res_df=HSPC_nRBC_enrichgo_res@compareClusterResult %>% group_by(Cluster) %>% do(head(.,15))
top10_enrichgo_res_df$ratio=as.numeric(data.frame(strsplit(top10_enrichgo_res_df$GeneRatio,split = '/'))[1,])/as.numeric(data.frame(strsplit(top10_enrichgo_res_df$GeneRatio,split = '/'))[2,])
top10_enrichgo_res_df$Description=factor(top10_enrichgo_res_df$Description,levels =unique(top10_enrichgo_res_df$Description) )
p=ggplot(top10_enrichgo_res_df,aes(x=Cluster,y=Description,color=-log10(p.adjust),size=ratio))+geom_point()+theme_bw()+scale_color_gradient(low = '#4387B5',high = 'firebrick3')
p 

ggsave(p,filename='Protein_NRBC_marker/res_pic/main_figure2/comapreenrichGO_HSPC_derived_nRBC.pdf',width = 8,height = 10)


#------------------------------------------------------ HSPC_derived fetal vs adult DE analysis---------------------------------------------------#
if(F){
  # fetal NRBC: YS、FL、FBM，三个阶段各抽取1:1:1的细胞，构成early、mid、late NRBC， 与整体不抽样分析，结果差异很小
    
  FL_early_Ery_filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu,Ery_stage=='early_Ery' & tissue_stage=='FL',downsample =1500)
  FL_mid_Ery_filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu,Ery_stage=='mid_Ery' & tissue_stage=='FL',downsample =1500)
  FL_late_Ery_filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu,Ery_stage=='late_Ery' & tissue_stage=='FL',downsample =1500)
  FL_filt_NBRC_altas_seu=merge(FL_early_Ery_filt_NBRC_altas_seu,c(FL_mid_Ery_filt_NBRC_altas_seu,FL_late_Ery_filt_NBRC_altas_seu))
  rm(FL_early_Ery_filt_NBRC_altas_seu,FL_mid_Ery_filt_NBRC_altas_seu,FL_late_Ery_filt_NBRC_altas_seu)
  
  FBM_early_Ery_filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu,Ery_stage=='early_Ery' & tissue_stage=='FBM',downsample =1500)
  FBM_mid_Ery_filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu,Ery_stage=='mid_Ery' & tissue_stage=='FBM',downsample =1500)
  FBM_late_Ery_filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu,Ery_stage=='late_Ery' & tissue_stage=='FBM',downsample =1500)
  FBM_filt_NBRC_altas_seu=merge(FBM_early_Ery_filt_NBRC_altas_seu,c(FBM_mid_Ery_filt_NBRC_altas_seu,FBM_late_Ery_filt_NBRC_altas_seu))
  rm(FBM_early_Ery_filt_NBRC_altas_seu,FBM_mid_Ery_filt_NBRC_altas_seu,FBM_late_Ery_filt_NBRC_altas_seu)
  
  ABM_early_Ery_filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu,Ery_stage=='early_Ery' & tissue_stage=='ABM',downsample =3000)
  ABM_mid_Ery_filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu,Ery_stage=='mid_Ery' & tissue_stage=='ABM',downsample =3000)
  ABM_late_Ery_filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu,Ery_stage=='late_Ery' & tissue_stage=='ABM',downsample =3000)
  ABM_filt_NBRC_altas_seu=merge(ABM_early_Ery_filt_NBRC_altas_seu,c(ABM_mid_Ery_filt_NBRC_altas_seu,ABM_late_Ery_filt_NBRC_altas_seu))
  rm(ABM_early_Ery_filt_NBRC_altas_seu,ABM_mid_Ery_filt_NBRC_altas_seu,ABM_late_Ery_filt_NBRC_altas_seu)
  
  #subset_filt_NBRC_altas_seu=merge(ABM_filt_NBRC_altas_seu,c(FBM_filt_NBRC_altas_seu,FL_filt_NBRC_altas_seu,YS_filt_NBRC_altas_seu))
  #rm(ABM_filt_NBRC_altas_seu,FBM_filt_NBRC_altas_seu,FL_filt_NBRC_altas_seu,YS_filt_NBRC_altas_seu)
  
  subset_filt_NBRC_altas_seu=merge(ABM_filt_NBRC_altas_seu,c(FBM_filt_NBRC_altas_seu,FL_filt_NBRC_altas_seu))
  rm(ABM_filt_NBRC_altas_seu,FBM_filt_NBRC_altas_seu,FL_filt_NBRC_altas_seu)
  
  gc()
  
  subset_filt_NBRC_altas_seu <- JoinLayers(subset_filt_NBRC_altas_seu)
  subset_filt_NBRC_altas_seu$type_stage='fetal'
  subset_filt_NBRC_altas_seu$type_stage[subset_filt_NBRC_altas_seu$tissue_stage=='ABM']='adult'
  Idents(subset_filt_NBRC_altas_seu)='type_stage'
  subset_filt_NBRC_altas_seu$source_celltype=factor(subset_filt_NBRC_altas_seu$source_celltype,levels =c( "FL_BFUE/CFUE","FL_ProE","FL_Bas","FL_Poly", "FL_Orth","FBM_BFUE/CFUE","FBM_ProE","FBM_Bas","FBM_Poly",  "FBM_Orth",
                                                                                                          "ABM_BFUE/CFUE", "ABM_ProE","ABM_Bas","ABM_Poly","ABM_Orth") )
  saveRDS(subset_filt_NBRC_altas_seu,file = 'Protein_NRBC_marker/res_data/temp_subset_filt_NBRC_altas_seu.rds')
  
}else{
  subset_filt_NBRC_altas_seu=readRDS('Protein_NRBC_marker/res_data/temp_subset_filt_NBRC_altas_seu.rds')
}


# 或者采用FL、FBM以及ABM中ABM vs FL/FBM 得到的degs list
fetal_adult_NRBC_whole_marker=FindAllMarkers(subset_filt_NBRC_altas_seu,group.by = 'type_stage')

fetal_marker_list=fetal_adult_NRBC_whole_marker[fetal_adult_NRBC_whole_marker$cluster=='fetal','avg_log2FC']
names(fetal_marker_list)=fetal_adult_NRBC_whole_marker[fetal_adult_NRBC_whole_marker$cluster=='fetal','gene']
fetal_marker_list=sort(fetal_marker_list,decreasing = T)
fetal_adult_gsego_res=gseGO(geneList = fetal_marker_list,OrgDb = org.Hs.eg.db,ont = 'BP',keyType = 'ALIAS')
dotplot(fetal_adult_gsego_res)
saveRDS(fetal_adult_gsego_res,file='Protein_NRBC_marker/res_data/main_figure2/wholelevel_fetal_adult_gsego_res.rds')

fetal_adult_gsego_res_df=fetal_adult_gsego_res@result;fetal_adult_gsego_res_df$res='up';fetal_adult_gsego_res_df$res[fetal_adult_gsego_res_df$NES <0]='down'
fetal_adult_gsego_res_df=fetal_adult_gsego_res_df[fetal_adult_gsego_res_df$p.adjust <0.05,] %>% group_by(res) %>% do(head(.,10))
fetal_adult_gsego_res_df=fetal_adult_gsego_res_df[order(fetal_adult_gsego_res_df$res,fetal_adult_gsego_res_df$NES),]
fetal_adult_gsego_res_df$Description=factor(fetal_adult_gsego_res_df$Description,levels = fetal_adult_gsego_res_df$Description)
p=ggplot(fetal_adult_gsego_res_df,aes(x=res,y=Description,size=-log(p.adjust),color=NES))+geom_point()+theme_bw()+scale_color_gradient(low = '#4387B5',high = 'firebrick3')+ggtitle(label = 'Fetal vs Adult')
p
ggsave(p,filename='Protein_NRBC_marker/res_pic/main_figure2/top10_fetal_adult_gseGO.pdf',width = 6,height = 6)



Idents(subset_filt_NBRC_altas_seu)='type_stage'
group='fetal_adult'
sfile='Protein_NRBC_marker/DE_marker/fetal_adult_all_Ery_RNA_markers.csv'
res=find_DEGs_func(seu = subset(subset_filt_NBRC_altas_seu,type_stage %in% c('fetal','adult')),group = group,sfile = sfile)
sub_fetal_adult_all_Ery_tissue_markers=res[[1]]
sub_fetal_addult_count_df=res[[2]]
rm(res);gc()
sub_fetal_adult_all_Ery_tissue_markers=sub_fetal_adult_all_Ery_tissue_markers[!sub_fetal_adult_all_Ery_tissue_markers$gene %in% names(HIS_genes),]

late_positive_sub_fetal_adult_all_Ery_tissue_markers=sub_fetal_adult_all_Ery_tissue_markers[sub_fetal_adult_all_Ery_tissue_markers$celltype=='late_Ery' & sub_fetal_adult_all_Ery_tissue_markers$avg_log2FC >0,]
late_positive_sub_fetal_adult_all_Ery_tissue_markers=late_positive_sub_fetal_adult_all_Ery_tissue_markers[order(late_positive_sub_fetal_adult_all_Ery_tissue_markers$avg_log2FC,decreasing = T),]
top_late_positive_sub_fetal_adult_all_Ery_tissue_markers=late_positive_sub_fetal_adult_all_Ery_tissue_markers[1:100,]
top_late_positive_sub_fetal_adult_all_Ery_tissue_markers=top_late_positive_sub_fetal_adult_all_Ery_tissue_markers[order(top_late_positive_sub_fetal_adult_all_Ery_tissue_markers$avg_log2FC*top_late_positive_sub_fetal_adult_all_Ery_tissue_markers$pct.1,decreasing = T),]
DotPlot(filt_NBRC_altas_seu,group.by = 'source_celltype',late_positive_sub_fetal_adult_all_Ery_tissue_markers$gene[1:50])+RotatedAxis()+scale_color_gradient(low = 'gray',high = 'firebrick3')
top_late_positive_sub_fetal_adult_all_Ery_tissue_markers=top_late_positive_sub_fetal_adult_all_Ery_tissue_markers[1:8,]
top_late_positive_sub_fetal_adult_all_Ery_tissue_markers=top_late_positive_sub_fetal_adult_all_Ery_tissue_markers[!top_late_positive_sub_fetal_adult_all_Ery_tissue_markers$gene %in% ABM_cho_topmarkers,]
DotPlot(filt_NBRC_altas_seu,group.by = 'source_celltype',top_late_positive_sub_fetal_adult_all_Ery_tissue_markers$gene)+RotatedAxis()+scale_color_gradient(low = 'gray',high = 'firebrick3')
# LGALS3 
ABM_cho_topmarkers=c(ABM_cho_topmarkers[1:9],'LGALS3',ABM_cho_topmarkers[10:12])
p=DotPlot(filt_NBRC_altas_seu,group.by = 'source_celltype',features =c(FL_top10_markers,FBM_unique_genes,top10_fetal_unique_marker_genes,ABM_cho_topmarkers),cols = c('gray','firebrick3'),scale = F)+RotatedAxis() #   colorRampPalette(colors = c('gray','firebrick3'))(100)
p
ggsave(p,filename='Protein_NRBC_marker/res_pic/main_figure2/top_marker_HSPC_derived_nRBC.pdf',width = 16,height = 8)


fetal_positive_all_Ery_tissue_markers=sub_fetal_adult_all_Ery_tissue_markers[sub_fetal_adult_all_Ery_tissue_markers$cluster=='fetal',]
fetal_positive_all_Ery_tissue_markers=fetal_positive_all_Ery_tissue_markers[fetal_positive_all_Ery_tissue_markers$avg_log2FC >0,]
fetal_positive_all_Ery_tissue_markers=fetal_positive_all_Ery_tissue_markers[fetal_positive_all_Ery_tissue_markers$pct.1>0.1,]
fetal_positive_all_Ery_tissue_markers=fetal_positive_all_Ery_tissue_markers[order(fetal_positive_all_Ery_tissue_markers$avg_log2FC,decreasing = T),]
top_fetal_positive_all_Ery_tissue_markers=fetal_positive_all_Ery_tissue_markers[1:50,]
top_fetal_positive_all_Ery_tissue_markers=top_fetal_positive_all_Ery_tissue_markers[order(top_fetal_positive_all_Ery_tissue_markers$pct.1,decreasing = T),]
DotPlot(filt_NBRC_altas_seu,features = unique(top_fetal_positive_all_Ery_tissue_markers$gene),cols = c('gray','firebrick3'),scale=F)+RotatedAxis()

top_fetal_unique_marker_genes=c(top_fetal_unique_marker_genes[1:9],'MEG3','GATA5','LIN28B','HMGA2',top_fetal_unique_marker_genes[10:11])

defintive_markers=list(FL_top10_markers=FL_top10_markers,FBM_unique_genes=FBM_unique_genes,top_fetal_unique_marker_genes=top_fetal_unique_marker_genes,ABM_cho_topmarkers=ABM_cho_topmarkers)
p=DotPlot(subset(filt_NBRC_altas_seu,NRBC_type=='definitive'),group.by = 'source_celltype',features =as.character(unlist(defintive_markers)),cols = c('gray','firebrick3'),scale = F)+RotatedAxis() #   colorRampPalette(colors = c('gray','firebrick3'))(100)
p
ggsave(p,filename='Protein_NRBC_marker/res_pic/main_figure2/top_marker_HSPC_derived_nRBC.pdf',width = 16,height = 8)

saveRDS(defintive_markers,file = 'Protein_NRBC_marker/res_data/main_figure2/defintive_markers.rds')


sub_fetal_adult_subcelltype_gseGO_list=subcelltype_gseGO_func(RNA_markers =sub_fetal_adult_all_Ery_tissue_markers[sub_fetal_adult_all_Ery_tissue_markers$cluster=='fetal',] )
sub_fetal_adult_subcelltype_gseGO_list[[2]]$group='fetal_adult'
saveRDS(sub_fetal_adult_subcelltype_gseGO_list,file = 'Protein_NRBC_marker/res_data/main_figure2/sub_fetal_adult_subcelltype_gseGO_list.rds')

top_sub_fetal_adult_subcelltype_gseGO_df=sub_fetal_adult_subcelltype_gseGO_list[[2]]
top_sub_fetal_adult_subcelltype_gseGO_df=top_sub_fetal_adult_subcelltype_gseGO_df[top_sub_fetal_adult_subcelltype_gseGO_df$ONTOLOGY=='BP', ]
top_sub_fetal_adult_subcelltype_gseGO_df$ratio=sapply(strsplit(top_sub_fetal_adult_subcelltype_gseGO_df$core_enrichment,split = '/'), length)/top_sub_fetal_adult_subcelltype_gseGO_df$setSize

top_sub_fetal_adult_subcelltype_gseGO_df$res='up'
top_sub_fetal_adult_subcelltype_gseGO_df$res[top_sub_fetal_adult_subcelltype_gseGO_df$NES <0]='down'
top_sub_fetal_adult_subcelltype_gseGO_df=top_sub_fetal_adult_subcelltype_gseGO_df %>% group_by(celltype,res) %>% do(head(.,10))
top_sub_fetal_adult_subcelltype_gseGO_df$celltype=factor(top_sub_fetal_adult_subcelltype_gseGO_df$celltype,levels = c('early_Ery','mid_Ery','late_Ery'))

top_sub_fetal_adult_subcelltype_gseGO_df=top_sub_fetal_adult_subcelltype_gseGO_df[order(top_sub_fetal_adult_subcelltype_gseGO_df$celltype,top_sub_fetal_adult_subcelltype_gseGO_df$NES),]
top_sub_fetal_adult_subcelltype_gseGO_df$Description=factor(top_sub_fetal_adult_subcelltype_gseGO_df$Description,levels = unique(top_sub_fetal_adult_subcelltype_gseGO_df$Description))
p5=ggplot(top_sub_fetal_adult_subcelltype_gseGO_df,aes(x=celltype,y=Description,size=ratio,color=NES))+geom_point()+scale_color_gradient2(low = 'navy',mid = 'white',high = 'firebrick3')+
  theme_bw()+theme(axis.text.x = element_text(angle = 45,hjust = 1,face = 'bold'))+ggtitle('fetal vs adult nRBC:gseGO of DEGs')
p5

ggsave(p5,width =6 ,height =8,filename='Protein_NRBC_marker/res_pic/main_figure2/DEGS_gseGOBP_fetal_adult_dotplot.pdf' )

###################################################################################################################################
#---------------------------the ratio of cell cycle phase--------------------------------# 
###################################################################################################################################
phase_stats_by_sample <- filt_NBRC_altas_seu@meta.data[filt_NBRC_altas_seu$tissue_stage!='YS',c('id','Phase','final_celltype','tissue_stage')] %>%
  group_by(tissue_stage, id, final_celltype, Phase) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(tissue_stage, id, final_celltype) %>%
  mutate(
    total = sum(count),
    proportion = count / total * 100,
    proportion_label = sprintf("%.1f%%", proportion)
  ) %>%
  ungroup()
# 选则样本中至少10个NRBC 存在
p_final <- ggplot(phase_stats_by_sample[phase_stats_by_sample$total >10 & phase_stats_by_sample$Phase=='G1',], aes(x = tissue_stage, y = proportion, fill = tissue_stage)) +
  geom_boxplot(width = 0.6, outlier.shape = 19, outlier.size = 1,alpha=0.7) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.7) +facet_wrap(~final_celltype, scales = "free_y", nrow = 1) +theme_classic() +
  scale_fill_manual(values = cols, name = "Tissue Stage") +labs(title = "Cell Cycle G1 Phase Distribution in subcelltype",x = "",y = "Proportion (%)") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 11),axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),legend.position = "bottom",
        strip.background = element_rect(fill = "lightgray"),strip.text = element_text(size = 11, face = "bold")) +ylim(c(0,115))+
  stat_compare_means(comparisons = list(c("FL", "FBM"),c("FL", "ABM"), c("FBM", "ABM")),method = "t.test",label = "p.signif" )

print(p_final)

# 计算样本数目
for( celltype in unique(phase_stats_by_sample$final_celltype)){
  print(celltype)
  print(table(phase_stats_by_sample[phase_stats_by_sample$total >10 & phase_stats_by_sample$Phase=='G1' &phase_stats_by_sample$final_celltype==celltype ,'tissue_stage']))
}

ggsave(p_final,filename = 'Protein_NRBC_marker/res_pic/main_figure2/definitive_G1_phase_wilcox_test.pdf',height = 6,width = 18)

p_final <- ggplot(phase_stats_by_sample[phase_stats_by_sample$total >10 & phase_stats_by_sample$Phase=='S',], aes(x = tissue_stage, y = proportion, fill = tissue_stage)) +
  geom_boxplot(width = 0.6, outlier.shape = 19, outlier.size = 1,alpha=0.7) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.7) +facet_wrap(~final_celltype, scales = "free_y", nrow = 1) +theme_classic() +
  scale_fill_manual(values = cols, name = "Tissue Stage") +labs(title = "Cell Cycle S Phase Distribution in subcelltype",x = "",y = "Proportion (%)") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 11),axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),legend.position = "bottom",
        strip.background = element_rect(fill = "lightgray"),strip.text = element_text(size = 11, face = "bold")) +ylim(c(0,100))+
  stat_compare_means(comparisons = list(c("FL", "FBM"),c("FL", "ABM"), c("FBM", "ABM")),method = "t.test",label = "p.signif" )

print(p_final)

ggsave(p_final,filename = 'Protein_NRBC_marker/res_pic/main_figure2/definitive_S_phase_wilcox_test.pdf',height = 6,width = 18)


p_final <- ggplot(phase_stats_by_sample[phase_stats_by_sample$total >10 & phase_stats_by_sample$Phase=='G2M',], aes(x = tissue_stage, y = proportion, fill = tissue_stage)) +
  geom_boxplot(width = 0.6, outlier.shape = 19, outlier.size = 1,alpha=0.7) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.7) +facet_wrap(~final_celltype, scales = "free_y", nrow = 1) +theme_classic() +
  scale_fill_manual(values = cols, name = "Tissue Stage") +labs(title = "Cell Cycle G2M Phase Distribution in subcelltype",x = "",y = "Proportion (%)") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 11),axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),legend.position = "bottom",
        strip.background = element_rect(fill = "lightgray"),strip.text = element_text(size = 11, face = "bold")) +ylim(c(0,115))+
  stat_compare_means(comparisons = list(c("FL", "FBM"),c("FL", "ABM"), c("FBM", "ABM")),method = "t.test",label = "p.signif" )

print(p_final)

ggsave(p_final,filename = 'Protein_NRBC_marker/res_pic/main_figure2/definitive_G2M_phase_wilcox_test.pdf',height = 6,width = 18)






