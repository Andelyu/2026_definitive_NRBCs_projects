

source('/home/gibh/2021_NRBC_chlyu/zx_lab_NRBC/scripts/scRNAseq_pipline/scRNAseq_analysis_model.R')


####################################################################################################################################
#----------------------------bulk NRBC RNA-Seq as the singleR referrence ---------------------------------#
####################################################################################################################################
if(F){
  #做相似性热图分析时，建议使用nr之后的log数据进行相似性分析
  library(pheatmap)
  
  #Li J, Hale J, Bhagia P, Xue F et al. Isolation and transcriptome analyses of human erythroid progenitors: BFU-E and CFU-E. Blood 2014 Dec 4;124(24):3636-45. PMID: 25339359
        
        
        
        
  An_nrbc_df1=read.table('/home/gibh/2021_NRBC_chlyu/ref_data/bulk_RNAseq/Anxiuli_lab_bulkRNAseq/GSE61566_genes.count_tracking.txt',header = T)
  rownames(An_nrbc_df1)=An_nrbc_df1$tracking_id
  An_nrbc_df1=An_nrbc_df1[c('CD34_count','BFU_count','CFU_count','Pro_count')]
  colnames(An_nrbc_df1)=gsub(pattern = '_count',replacement = '',x = colnames(An_nrbc_df1))
  An_nrbc_df1=round(An_nrbc_df1)
  
  #	An X, Schulz VP, Li J, Wu K et al. Global transcriptome analyses of human and murine terminal erythroid differentiation. Blood 2014 May 29;123(22):3466-77. PMID: 24637361
        
        
        
        
  An_nrbc_df2=read.table('/home/gibh/2021_NRBC_chlyu/ref_data/bulk_RNAseq/Anxiuli_lab_bulkRNAseq/GSE53983_All_hs_countData.txt',header = T)
  rownames(An_nrbc_df2)=An_nrbc_df2$Gene;An_nrbc_df2=An_nrbc_df2[,-1]
  colnames(An_nrbc_df2)=c('ProE_1','ProE_2','ProE_3','eBas_1','eBas_2','eBas_3','lBas_1','lBas_2','lBas_3','Poly_1','Poly_2','Poly_3','Orth_1','Orth_2','Orth_3')
  
  #Schulz VP, Yan H, Lezon-Geyda K, An X et al. A Unique Epigenomic Landscape Defines Human Erythropoiesis. Cell Rep 2019 Sep 10;28(11):2996-3009.e7.
  An_nrbc_df3=read.table('/home/gibh/2021_NRBC_chlyu/ref_data/bulk_RNAseq/Anxiuli_lab_bulkRNAseq/GSE128268_ProgenitorCounts.txt',header = T)
  
  shared_gene_ids=intersect(intersect(rownames(An_nrbc_df1),rownames(An_nrbc_df2)),rownames(An_nrbc_df3));length(shared_gene_ids)
  all_bulk_rnaseq_df=cbind(cbind(An_nrbc_df1[shared_gene_ids,],An_nrbc_df2[shared_gene_ids,]),An_nrbc_df3[shared_gene_ids,])
  pheatmap(cor(log2(cpm(all_bulk_rnaseq_df)+1)),main = 'GSE61566 &GSE53983 & GSE128268') 
  write.csv(all_bulk_rnaseq_df,file = 'ref_data/NRBC_ref_bullk_data.csv')
  
  
  shared_ids=intersect(rownames(An_nrbc_df2),rownames(An_nrbc_df3))
  test_df=cbind(An_nrbc_df2[shared_ids,],An_nrbc_df3[shared_ids,])
  pheatmap(cor(log2(cpm(test_df)+1)),main = 'GSE53983 & GSE128268')
  
  shared_ids=intersect(rownames(An_nrbc_df1),rownames(An_nrbc_df3))
  test_df2=cbind(An_nrbc_df1[shared_ids,],An_nrbc_df3[shared_ids,])
  pheatmap(cor(log2(cpm(test_df2)+1)),main = 'GSE61566 & GSE128268')
  
  # 考虑三批样本之间的批次效应的时候，注意，做批次效应分析之前数据一般必须进行normalization；
  library(limma);library(edgeR)
  colnames(all_bulk_rnaseq_df)[1:4]=c('HSPC','BFUE','CFUE','ProE')
  type=strsplit2(colnames(all_bulk_rnaseq_df),split = '_')[,1]
  type_num=1:8;names(type_num)=unique(type)
  design= data.frame(row.names = colnames(all_bulk_rnaseq_df),type=as.numeric(type_num[type]))
  rb_all_bulk_rnaseq_df=removeBatchEffect(x =log2(cpm(all_bulk_rnaseq_df)+1),batch =  as.factor(c(rep('barch1',4),rep('batch2',dim(An_nrbc_df2)[2]),rep('batch3',dim(An_nrbc_df3)[2]))),design =design )
  pheatmap(cor(rb_all_bulk_rnaseq_df),main = 'removeBatchEffect GSE61566 &GSE53983 & GSE128268')
  write.csv(rb_all_bulk_rnaseq_df,file = 'ref_data/NRBC_ref_bullk_removebatcheffect.csv')
  
  
  # 构建se对象
  library(SummarizedExperiment)
  colData=data.frame(row.names =colnames(test_df),celltype=strsplit2(colnames(test_df),split = '_')[,1])
  metadata='GSE53983 & GSE128268 RNA-seq'
  nrbc_ref_se <- SummarizedExperiment(assays=list(logcounts=log2(cpm(test_df)+1)),
                                      colData=colData,
                                      metadata=metadata)
  library(limma)
  library(edgeR)
  colData=data.frame(row.names =colnames(all_bulk_rnaseq_df),celltype=strsplit2(colnames(all_bulk_rnaseq_df),split = '_')[,1])
  metadata='GSE61566 & GSE53983 & GSE128268 RNA-seq'
  nrbc_ref_se2 <- SummarizedExperiment(assays=list(logcounts=log2(cpm(all_bulk_rnaseq_df)+1)),
                                       colData=colData,
                                       metadata=metadata)
  # SeuratData中的数据
  library(SeuratData)
  library(celldex)
  hema.se=NovershternHematopoieticData()
  lable_new=list( 'Hematopoietic stem cells_CD133+ CD34dim'='HSC1',
                  'Hematopoietic stem cells_CD38- CD34+'='HSC2',                 
                  'Colony Forming Unit-Granulocytes'='Gran.CFU',
                  'Colony Forming Unit-Megakaryocytic'='Meg.CFU',
                  'Colony Forming Unit-Monocytes'='Mon.CFU',
                  'Common myeloid progenitors'='CMPs',
                  
                  'Megakaryocyte/erythroid progenitors'='MEP',
                  'Erythroid_CD34+ CD71+ GlyA-'='Ery1',
                  'Erythroid_CD34- CD71+ GlyA-'='Ery2',
                  'Erythroid_CD34- CD71+ GlyA+'='Ery3',
                  'Erythroid_CD34- CD71lo GlyA+'='Ery4',
                  'Erythroid_CD34- CD71- GlyA+'='Ery5',
                  
                  'Early B cells'='B.early',
                  'Naive B cells'='B.naive',
                  'Pro B cells'='B.pro',
                  'Mature B cells'='mB',
                  'Mature B cells class able to switch'='mB.switch',
                  'Mature B cells class switched'='mB.switched',
                  
                  'CD4+ Central Memory'='T.CD4_cM',
                  'CD4+ Effector Memory'='T.CD4_eM',
                  'Naive CD4+ T cells'='T.CD4_na',
                  'CD8+ Central Memory'='T.CD8_cM',
                  'CD8+ Effector Memory'='T.CD8_eM',
                  'CD8+ Effector Memory RA'='T.CD8_eMRA',
                  'Naive CD8+ T cells'='T.CD8_na',
                  
                  'Granulocyte/monocyte progenitors'='GMPs',
                  'Granulocytes (Neutrophilic Metamyelocytes)'='Gran.Meta',
                  'Granulocytes (Neutrophils)'='Gran',
                  'Eosinophils'='Eos',
                  'Basophils'='Baso',
                  'Monocytes'='Mon',
                  
                  'Mature NK cells_CD56- CD16- CD3-'='NK3',
                  'Mature NK cells_CD56- CD16+ CD3-'='NK1',
                  'Mature NK cells_CD56+ CD16+ CD3-'='NK2',
                  'NK T cells'='NK4',
                  
                  'Megakaryocytes'='Meg',
                  'Myeloid Dendritic Cells'='DC.Myeloid',
                  'Plasmacytoid Dendritic Cells'='DC.Plasma'
                  
  )
  hema.se$celltype=as.character(lable_new[hema.se$label.fine])
  save(hema.se,file = '/home/gibh/2021_NRBC_chlyu/ref_data/hema_ref_bullk_RNAseq_se.Rdata')
  
  save(nrbc_ref_se,nrbc_ref_se2,file = '/home/gibh/2021_NRBC_chlyu/ref_data/NRBC_ref_bullk_RNAseq_se.Rdata')
  
  
  
  
}else(load('/home/gibh/2021_NRBC_chlyu/ref_data/NRBC_ref_bullk_RNAseq_se.Rdata',verbose = T))


####################################################################################################################################
#------------------------------- human_Fetal_scRNAseq altas data as the reference------------------------------ #
####################################################################################################################################

# 不同数据来源数据，基因注释存在差异，统一gene symbol 信息
#YS altas 数据中，多个样本来源 基因注释信息差异比较大，我们的数据，ABM
all_shared_ensembl_id_info=read.csv('ref_data/all_shared_ensembl_id_info.csv',header=T)
all_shared_ensembl_id_info$org_symbol=as.character(mapIds(x = org.Hs.eg.db,keys = all_shared_ensembl_id_info$X,keytype = 'ENSEMBL',column = 'SYMBOL'))
all_shared_ensembl_id_info$org_symbol[is.na(all_shared_ensembl_id_info$org_symbol)]=all_shared_ensembl_id_info$ref_symbol[is.na(all_shared_ensembl_id_info$org_symbol)]
table(all_shared_ensembl_id_info$ref_symbol==all_shared_ensembl_id_info$org_symbol) # F:2482,T:33248 
write.csv(all_shared_ensembl_id_info,file = 'ref_data/all_shared_ensembl_id_info.csv')


all_shared_ensembl_id_info=read.csv('ref_data/all_shared_ensembl_id_info.csv',header=T)

#----------------------Mapping the developing human immune system across organs--------------------------------------#
# 来自同一个大项目,其基因注释信息保持一致
if(F){
    Convert('ref_data/ref_scRNAseq_data/human_cell_atlas/human_immune_system_across_organs/PAN.A01.v01.raw_count.20210429.PFI.embedding.h5ad',dest="h5seurat") 
    human_Fetal_altas_seu=LoadH5Seurat('ref_data/ref_scRNAseq_data/human_cell_atlas/human_immune_system_across_organs/PAN.A01.v01.raw_count.20210429.PFI.embedding.h5seurat', meta.data = FALSE,reduction='umap')
    human_Fetal_altas_seu_meta=read.csv('ref_data/ref_scRNAseq_data/human_cell_atlas/human_immune_system_across_organs/PAN.A01.v01.entire_data_normalised_log.20210429.full_obs.annotated.clean.csv',header =T )
    
    rownames(human_Fetal_altas_seu_meta)=human_Fetal_altas_seu_meta$X;human_Fetal_altas_seu_meta=human_Fetal_altas_seu_meta[,-1]
    dim(human_Fetal_altas_seu_meta)#  908178
    table(colnames(human_Fetal_altas_seu) %in%  rownames(human_Fetal_altas_seu_meta))# FALSE   TRUE ：3695 908178 
    human_Fetal_altas_seu=subset(human_Fetal_altas_seu,cells =rownames(human_Fetal_altas_seu_meta) )
    human_Fetal_altas_seu@meta.data=human_Fetal_altas_seu_meta[colnames(human_Fetal_altas_seu),]
    
    DimPlot(human_Fetal_altas_seu,group.by = c('organ','uniform_label_lvl0'),raster=FALSE)
    table(human_Fetal_altas_seu$organ)
    
    ys_fl_fbm_seu=subset(human_Fetal_altas_seu, organ %in% c('YS','LI','BM'))
    table(ys_fl_fbm_seu$organ)
    table(ys_fl_fbm_seu$uniform_label )
    rm(human_Fetal_altas_seu)
    
    
    table(ys_fl_fbm_seu$uniform_label_lvl0)
    table(ys_fl_fbm_seu$anno_lvl_2_final_clean)
    
    DimPlot(ys_fl_fbm_seu,group.by = c('organ','uniform_label_lvl0'),raster=FALSE,cols = cols)
    DimPlot(ys_fl_fbm_seu,group.by = c('organ','anno_lvl_2_final_clean'),raster=FALSE)
    
    library(rjson)
    celltype_dict=fromJSON('ref_data/ref_scRNAseq_data/humnan_cell_atlas/human_immune_system_across_organs/broad_annotation_dict.json')
    celltype_dict_lelves=as.character(unlist(sapply(names(celltype_dict),function(x){rep(x,length(celltype_dict[[x]]))})))
    celltype_dict_df=data.frame(celltype_dict=as.character(unlist(celltype_dict)),celltype_dict_lelves=celltype_dict_lelves)
    ys_fl_fbm_seu$anno_lvl_1_final_clean=celltype_dict_df[match(x = ys_fl_fbm_seu$anno_lvl_2_final_clean,table = celltype_dict_df$celltype_dict),'celltype_dict_lelves']
    ys_fl_fbm_seu$anno_lvl_1_final_clean[ys_fl_fbm_seu$anno_lvl_2_final_clean %in% c("EARLY_MK","LATE_MK" )]='MEGAKARYOCYTE'
    
    celltype_dict_lelves_color=read.csv('ref_data/ref_scRNAseq_data/humnan_cell_atlas/human_immune_system_across_organs/broad_annotation_colors.csv')
    ys_fl_fbm_seu$anno_lvl_1_final_clean=factor(ys_fl_fbm_seu$anno_lvl_1_final_clean,levels = celltype_dict_lelves_color$celltype)
    DimPlot(ys_fl_fbm_seu,group.by = c('organ','anno_lvl_1_final_clean'),cols = celltype_dict_lelves_color$color,raster=FALSE)
    
    DimPlot(ys_fl_fbm_seu,group.by = 'uniform_label_lvl0',split.by = 'organ',cols = cols,raster=FALSE)
    
    ys_fl_fbm_seu=subset(ys_fl_fbm_seu,anno_lvl_1_final_clean %in% celltype_dict_lelves_color$celltype[-10])
    DimPlot(ys_fl_fbm_seu,group.by = 'anno_lvl_1_final_clean',cols = celltype_dict_lelves_color$color[-9],raster=FALSE,label = T)
    DimPlot(ys_fl_fbm_seu,group.by = 'anno_lvl_1_final_clean',split.by = 'organ',cols = celltype_dict_lelves_color$color[-9],raster=FALSE)
    
    
    ys_fl_fbm_seu$anno_final_celltype=as.character(ys_fl_fbm_seu$anno_lvl_1_final_clean)
    ys_fl_fbm_seu$anno_final_celltype[grep(pattern = 'DC',ys_fl_fbm_seu$anno_lvl_2_final_clean)]='DC'
    ys_fl_fbm_seu$anno_final_celltype[grep(pattern = 'MONOCYTE',ys_fl_fbm_seu$anno_lvl_2_final_clean)]='MONOCYTE'
    ys_fl_fbm_seu$anno_final_celltype[grep(pattern = 'MACROPHAGE',ys_fl_fbm_seu$anno_lvl_2_final_clean)]='MACROPHAGE'
    ys_fl_fbm_seu$anno_final_celltype[grep(pattern = 'YS_ERY',ys_fl_fbm_seu$anno_lvl_2_final_clean)]='YS_ERY'
    
    DimPlot(ys_fl_fbm_seu,group.by = 'anno_final_celltype',split.by = 'organ',raster=FALSE,cols = cols)
    DimPlot(subset(ys_fl_fbm_seu,organ == 'BM'),group.by = 'anno_final_celltype',split.by = 'organ',raster=FALSE,cols = cols,label = T )
    
    DimPlot(subset(ys_fl_fbm_seu,anno_final_celltype == 'MYELOID'),group.by = 'anno_lvl_2_final_clean',split.by = 'organ',raster=FALSE,cols = cols,label = T )
    table(subset(ys_fl_fbm_seu,anno_final_celltype == 'MYELOID')@meta.data[,'anno_lvl_2_final_clean'])
    #DEVELOPING_NEPHRON_I DEVELOPING_NEPHRON_II：23，1
    
    ys_fl_fbm_seu$anno_final_celltype[grep(pattern = 'MOP',ys_fl_fbm_seu$anno_lvl_2_final_clean)]='MOP'
    ys_fl_fbm_seu$anno_final_celltype[grep(pattern = 'NEUTROPHIL',ys_fl_fbm_seu$anno_lvl_2_final_clean)]='NEUTROPHIL'
    ys_fl_fbm_seu$anno_final_celltype[grep(pattern = 'OSTEOCLAST',ys_fl_fbm_seu$anno_lvl_2_final_clean)]='OSTEOCLAST'
    ys_fl_fbm_seu$anno_final_celltype[grep(pattern = 'MYELOCYTE',ys_fl_fbm_seu$anno_lvl_2_final_clean)]='MYELOCYTE'
    
    ys_fl_fbm_seu$anno_final_celltype[grep(pattern = 'YS_STROMA',ys_fl_fbm_seu$anno_lvl_2_final_clean)]='YS_STROMA'
    ys_fl_fbm_seu$anno_final_celltype[grep(pattern = 'NEPHRON',ys_fl_fbm_seu$anno_lvl_2_final_clean)]='NEPHRON' # 极少量，总24
    ys_fl_fbm_seu$anno_final_celltype[grep(pattern = 'FIBROBLAST_XVII',ys_fl_fbm_seu$anno_lvl_2_final_clean)]='MYE_FIBROBLAST' # 极少量，总24
    ys_fl_fbm_seu$anno_final_celltype[grep(pattern = 'MACROPHAGE_ERY',ys_fl_fbm_seu$anno_lvl_2_final_clean)]='MACROPHAGE_ERY' # 极少量，总24
    ys_fl_fbm_seu$anno_final_celltype[grep(pattern = 'LANGERHANS_CELLS',ys_fl_fbm_seu$anno_lvl_2_final_clean)]='MACROPHAGE' # 极少量，总24
    
    table(subset(ys_fl_fbm_seu,anno_final_celltype == 'MYELOID')@meta.data[,'anno_lvl_2_final_clean'])
    
    table(ys_fl_fbm_seu@meta.data[ys_fl_fbm_seu$anno_final_celltype=='TISSUE STROMA','anno_lvl_2_final_clean'])
    ys_fl_fbm_seu$anno_final_celltype[grep(pattern = 'HEPATOCYTE',ys_fl_fbm_seu$anno_lvl_2_final_clean)]='HEPATOCYTE' 
    ys_fl_fbm_seu$anno_final_celltype[grep(pattern = 'NEURON',ys_fl_fbm_seu$anno_lvl_2_final_clean)]='NEURON' 
    ys_fl_fbm_seu$anno_final_celltype[grep(pattern = 'OSTEOBLAST',ys_fl_fbm_seu$anno_lvl_2_final_clean)]='OSTEOBLAST' 
    ys_fl_fbm_seu$anno_final_celltype[grep(pattern = 'EPITHELIUM',ys_fl_fbm_seu$anno_lvl_2_final_clean)]='EPITHELIUM' 
    ys_fl_fbm_seu$anno_final_celltype[grep(pattern = 'CHONDROCYTE',ys_fl_fbm_seu$anno_lvl_2_final_clean)]='CHONDROCYTE' 
    ys_fl_fbm_seu$anno_final_celltype[grep(pattern = 'CHONDROCYTE',ys_fl_fbm_seu$anno_lvl_2_final_clean)]='CHONDROCYTE' 
    ys_fl_fbm_seu$anno_final_celltype[grep(pattern = 'MESOTHELIUM',ys_fl_fbm_seu$anno_lvl_2_final_clean)]='MESOTHELIUM' 
    ys_fl_fbm_seu$anno_final_celltype[grep(pattern = 'GLIAL',ys_fl_fbm_seu$anno_lvl_2_final_clean)]='GLIAL' 
    ys_fl_fbm_seu$anno_final_celltype[grep(pattern = 'KERATINOCYTE',ys_fl_fbm_seu$anno_lvl_2_final_clean)]='KERATINOCYTE' 
    
    # CHONDROCYTE： 软骨细胞，OSTEOCLAST：破骨细胞，OSTEOBLAST成骨细胞， KERATINOCYTE ：角质细胞
    
    DimPlot(ys_fl_fbm_seu,group.by = 'anno_final_celltype',raster=FALSE,cols = cols,label=T,label.size=3)
    DimPlot(ys_fl_fbm_seu,group.by = 'anno_final_celltype',split.by = 'organ',raster=FALSE,cols = cols,label.size=3)
    table(ys_fl_fbm_seu@meta.data[,c('organ','anno_final_celltype')])
    
    
    ys_fl_fbm_seu$anno_final_celltype2=ys_fl_fbm_seu$anno_final_celltype
    ys_fl_fbm_seu@meta.data[ys_fl_fbm_seu$anno_final_celltype2 %in% c('YS_ERY','ERYTHROID CELLS'),'anno_final_celltype2']=ys_fl_fbm_seu@meta.data[ys_fl_fbm_seu$anno_final_celltype2 %in% c('YS_ERY','ERYTHROID CELLS'),'anno_lvl_2_final_clean']
    
    saveRDS(ys_fl_fbm_seu,file = 'ref_data/ref_scRNAseq_data/ys_fl_fbm_seu.rds')
    
    # visualize the data 
    ys_fl_fbm_seu=NormalizeData(ys_fl_fbm_seu)
    FeaturePlot(ys_fl_fbm_seu,features = c('HBE1','HBZ','OAT','CD36'),split.by = 'organ',raster=FALSE)
    
    YS_Ery_vs_Ery_Markers=FindMarkers(ys_fl_fbm_seu,ident.1 = 'YS_ERY',ident.2 ='ERYTHROID CELLS' ,group.by = 'anno_final_celltype')
    FeaturePlot(ys_fl_fbm_seu,features = c('HBE1','HBZ','OAT','CD36'),split.by = 'organ',raster=FALSE)
    
    ys_fl_fbm_seu=NormalizeData(ys_fl_fbm_seu ) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData() %>% RunPCA()
    ys_fl_fbm_seu=RunHarmony(ys_fl_fbm_seu, group.by.vars=c('method'))
    ys_fl_fbm_seu=RunUMAP(ys_fl_fbm_seu,reduction.name = 'sub_umap',dims = 1:20,reduction = 'harmony')
    DimPlot(ys_fl_fbm_seu,group.by = c('anno_lvl_2_final_clean','organ','method'),cols = cols,pt.size = 0.6,alpha = 0.6,reduction = 'sub_umap')
    
  
    rm(ys_fl_fbm_seu,YS_Ery_seu)
    
 
    
}else{
  ys_fl_fbm_seu=readRDS('ref_data/ref_scRNAseq_data/ys_fl_fbm_seu.rds')
} 



