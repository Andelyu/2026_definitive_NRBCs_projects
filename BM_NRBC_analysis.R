set.seed(123)

library(SeuratDisk)
library(anndata)
library(Seurat)
library(RColorBrewer)
cols=c(brewer.pal(12,"Set3"),brewer.pal(6,"PiYG"),brewer.pal(6,"BrBG"),brewer.pal(8,"Set2"),
       brewer.pal(12,"Set3"),brewer.pal(8,"Pastel2"),brewer.pal(9,"Pastel1"),brewer.pal(8,"Accent"))
col=unique(cols)[-14]
setwd('/home/gibh/2021_NRBC_chlyu')
dir.create('NRBC_BM_altas/res_pic')
dir.create('NRBC_BM_altas/res_data')

source('/home/gibh/2021_NRBC_chlyu/zx_lab_NRBC/scripts/scRNAseq_pipline/scRNAseq_analysis_model.R')
load('/home/gibh/2021_NRBC_chlyu/ref_data/NRBC_ref_bullk_RNAseq_se.Rdata',verbose = T)  

########################################################################################################################################################### 
#-------------------------------------prepare  data------------------------------------
########################################################################################################################################################### 

all_shared_ensembl_id_info=read.csv('ref_data/all_shared_ensembl_id_info.csv',header=T);all_shared_ensembl_id_info=all_shared_ensembl_id_info[,-1]
rownames(all_shared_ensembl_id_info)=all_shared_ensembl_id_info$X
all_shared_ensembl_id_info=all_shared_ensembl_id_info[,-1]
table(duplicated(all_shared_ensembl_id_info$ref_symbol1)) #497 symbol replicated 

################################################################################################################# 
#--------------------------------------------------FBM  ALTAS DATA-------------------------------------------#
################################################################################################################# 
run_time=1
if(run_time==0){
  ys_fl_fbm_seu=readRDS('ref_data/ref_scRNAseq_data/ys_fl_fbm_seu.rds')
  FBM_altas_seu=subset(ys_fl_fbm_seu,organ=='BM');rm(ys_fl_fbm_seu)
  head(FBM_altas_seu)
  FBM_altas_seu=subset(FBM_altas_seu,predicted_doublets=='False')
  FBM_altas_seu$resource='E-MTAB-11343'
  table(rownames(FBM_altas_seu) %in% all_shared_ensembl_id_info$ref_symbol1) #  different genome annotation to fl_lab: FALSE/TRUE:10870/ 22668, 注意，seurat rownames中不使用_,而是采用-代替_
  table(rownames(FBM_altas_seu) %in% all_shared_ensembl_id_info$ref_symbol) # the same genome annotation to zx_lab: FALSE/TRUE:24/33514

  FBM_altas_seu_new=CreateSeuratObject(GetAssayData(FBM_altas_seu,assay = 'RNA',layer = 'counts'),meta.data = FBM_altas_seu@meta.data,min.cells = 10,min.features = 200)
  FBM_altas_seu_new[['umap']]=FBM_altas_seu[['umap']]
  rm(FBM_altas_seu);gc()
  FBM_altas_seu=FBM_altas_seu_new;rm(FBM_altas_seu_new);gc()
  
  table(rownames(FBM_altas_seu) %in% all_shared_ensembl_id_info$ref_symbol) # the same genome annotation to zx_lab: FALSE/TRUE:17/33514
  

  # 先整理不在本身注释基因组信息中，因为添加了后缀
  symbol_addpostfix=rownames(FBM_altas_seu)[!rownames(FBM_altas_seu) %in% all_shared_ensembl_id_info$ref_symbol]
  symbol_addpostfix=data.frame(strsplit2(symbol_addpostfix,split = '-'))[,1] 
  symbol_addpostfix_inf=all_shared_ensembl_id_info[all_shared_ensembl_id_info$ref_symbol %in% symbol_addpostfix,]
  symbol_addpostfix_inf=symbol_addpostfix_inf[order(symbol_addpostfix_inf$ref_symbol),]
  symbol_addpostfix_inf$number=1:dim(symbol_addpostfix_inf)[1]
  du_symbol_addpostfix=symbol_addpostfix_inf[match(symbol_addpostfix,symbol_addpostfix_inf$ref_symbol)+1,'ref_symbol1']
  names(du_symbol_addpostfix)=paste(symbol_addpostfix,'1',sep = '-')
  
  # 再处理在范围内的symbol
  uniq_symbols=rownames(FBM_altas_seu)[rownames(FBM_altas_seu) %in% all_shared_ensembl_id_info$ref_symbol]
  new_symbol_list=all_shared_ensembl_id_info[match(uniq_symbols,all_shared_ensembl_id_info$ref_symbol),'ref_symbol1']
  names(new_symbol_list)=uniq_symbols
  new_symbol_list=c(du_symbol_addpostfix,new_symbol_list)
  
  # 得到总list
  du_symbol=as.character(new_symbol_list)[duplicated(as.character(new_symbol_list))];length(du_symbol)
  tmp_assay=GetAssayData(FBM_altas_seu,assay = 'RNA',layer = 'counts')
  rownames(tmp_assay)=as.character(new_symbol_list[rownames(tmp_assay)])
  # 将symbol进行映射处理
  du_tmp_assay=tmp_assay[rownames(tmp_assay) %in% du_symbol,]
  tmp_assay=tmp_assay[!rownames(tmp_assay) %in% du_symbol,]
  
   # 再处理映射后重复symbol 
  du_tmp_assay=aggregate(du_tmp_assay,by=list(rownames(du_tmp_assay)), FUN=sum) 
  du_tmp_assay=column_to_rownames(du_tmp_assay,'Group.1')
  du_tmp_assay=as(as.matrix(du_tmp_assay),'dgCMatrix')
  #合并得到最后symbol处理矩阵
  tmp_assay=rbind(tmp_assay,du_tmp_assay)
  FBM_altas_seu=CreateSeuratObject(counts = tmp_assay,min.cells = 10,min.features = 200,meta.data =FBM_altas_seu@meta.data )
  rm(tmp_assay,du_tmp_assay)
  
  table(FBM_altas_seu$donor)

  FBM_altas_seu=NormalizeData(FBM_altas_seu) %>% FindVariableFeatures(nfeatures =3000 ) %>% ScaleData() %>% RunPCA() %>% RunHarmony(group.by.vars=c('method','donor','Sort_id')) 
  
  FBM_altas_seu=RunUMAP(FBM_altas_seu,reduction = 'harmony',dims = 1:30,reduction.name = 'umap2')
  DimPlot(FBM_altas_seu,group.by = 'subcelltype',reduction = 'umap2',cols = cols,label = T)
  
  
  FBM_altas_seu=RunUMAP(FBM_altas_seu,reduction = 'harmony',dims = 1:20,reduction.name = 'umap1',return.model = T)
  DimPlot(FBM_altas_seu,reduction = 'umap1',group.by = 'subcelltype',cols = col,label = T)  
  
  sort(table(FBM_altas_seu$subcelltype))
  FBM_altas_seu=subset(FBM_altas_seu,subcelltype %in% names(table(FBM_altas_seu$subcelltype))[table(FBM_altas_seu$subcelltype) >10])
  table(FBM_altas_seu$anno_lvl_2_final_clean[FBM_altas_seu$subcelltype=='PROGENITORS'])
  
  FBM_altas_seu$subcelltype[FBM_altas_seu$subcelltype=='PROGENITORS']= FBM_altas_seu$anno_lvl_2_final_clean[FBM_altas_seu$subcelltype=='PROGENITORS']
  FBM_altas_seu$subcelltype=gsub(pattern ='YS_Bas' ,replacement = 'YS_Bas/Poly',FBM_altas_seu$subcelltype)
  subcelltype_levels=c( 'HSC_MPP','CYCLING_MPP','CMP','GMP','LMPP_MLP','MEMP','CYCLING_MEMP','MEP' ,
                        "BFUE/CFUE","ProE","Bas","Poly","Orth","YS_Bas/Poly", "YS_Orth", "MEGAKARYOCYTE","EO/BASO/MAST",'NEUTROPHIL',
                       "MOP","MONOCYTE","DC","MACROPHAGE","MACROPHAGE_ERY","B CELLS", "ILC","NK/T CELLS","MYELOCYTE",
                       "ENDOTHELIUM","FIBROBLASTS","SKELETAL MUSCLE",
                       'CHONDROCYTE','OSTEOCLAST','OSTEOBLAST'
  );table(unique(FBM_altas_seu$subcelltype) %in%  subcelltype_levels)
  FBM_altas_seu$subcelltype=factor(FBM_altas_seu$subcelltype,levels =subcelltype_levels )
  
  DimPlot(FBM_altas_seu,group.by = 'subcelltype',reduction = 'umap2',cols = cols,label = F)+ggtitle('FBM ALTAS')
  table(FBM_altas_seu$donor) 
  
  
  
  saveRDS(FBM_altas_seu,file = 'NRBC_BM_altas/FBM_altas_seu.rds')# 93653 
  
  run_time=run_time+1
}else{
  FBM_altas_seu=readRDS('NRBC_BM_altas/FBM_altas_seu.rds')
}
table(rownames(FBM_altas_seu) %in% all_shared_ensembl_id_info$ref_symbol1)# FALSE/TRUE : 6/22875 
table(rownames(FBM_altas_seu) %in% all_shared_ensembl_id_info$ref_symbol)# FALSE/TRUE :4869/18012 



################################################################################################################# 
#--------------------------------------------------EBM NRBC data--------------------------------------------#
################################################################################################################# 
# 6 fetal BM erythopoisis,5 PCW: 3 buds , 8 PWC :3 longbone,GSE143753
# ref 8 fetal BM erythopoisis,5 PCW, 8 PWC,GSE143753

# create the EBM seuratObject 
if(F){
  # library(velocyto.R) #read.loom.matrices ，因为存储中只有一个matrix信息，数据为rawdata信息
  library(SeuratDisk)
  #Chromium Single cell 3’ Library and Gel Bead Kit V2 (10× Genomics
  # Convert('ref_data/ref_scRNAseq_data/GSE143753/GSM4274188_CS13_limbbud_rawdata.loom',dest = "h5seurat")
  
  fetal_bm_list=list()
  fetal_bm_limbbud_1=Connect(filename ='ref_data/ref_scRNAseq_data/GSE143753/GSM4274188_CS13_limbbud_rawdata.loom',mode = 'r' )
  fetal_bm_list[['limbbud_1']]=as.Seurat(fetal_bm_limbbud_1);Idents( fetal_bm_list[['limbbud_1']])='limbbud_1'
  fetal_bm_limbbud_2=Connect(filename ='ref_data/ref_scRNAseq_data/GSE143753/GSM4274189_CS15_limbbud_rawdata.loom',mode = 'r' )
  fetal_bm_list[['limbbud_2']]=as.Seurat(fetal_bm_limbbud_2);Idents( fetal_bm_list[['limbbud_2']])='limbbud_2'
  fetal_bm_limbbud_3=Connect(filename ='ref_data/ref_scRNAseq_data/GSE143753/GSM4274190_CS15_2__limbbud_rawdata.loom',mode = 'r' )
  fetal_bm_list[['limbbud_3']]=as.Seurat(fetal_bm_limbbud_3);Idents( fetal_bm_list[['limbbud_3']])='limbbud_3'
  
  fetal_bm_longbone_1=Connect(filename ='ref_data/ref_scRNAseq_data/GSE143753/GSM4274191_CS20_longbone_rawdata.loom',mode = 'r' )
  fetal_bm_list[['longbone_1']]=as.Seurat(fetal_bm_longbone_1);Idents( fetal_bm_list[['longbone_1']])='longbone_1'
  fetal_bm_longbone_2=Connect(filename ='ref_data/ref_scRNAseq_data/GSE143753/GSM4274192_CS22_longbone_rawdata.loom',mode = 'r' )
  fetal_bm_list[['longbone_2']]=fetal_bm_limbbud_3_seu=as.Seurat(fetal_bm_longbone_2);Idents( fetal_bm_list[['longbone_2']])='longbone_2'
  fetal_bm_longbone_3=Connect(filename ='ref_data/ref_scRNAseq_data/GSE143753/GSM4274193_CS22_2__longbone_rawdata.loom',mode = 'r' )
  fetal_bm_list[['longbone_3']]=as.Seurat(fetal_bm_longbone_3);Idents( fetal_bm_list[['longbone_3']])='longbone_3'
  
  rm(list = c('fetal_bm_limbbud_1','fetal_bm_limbbud_2','fetal_bm_limbbud_3','fetal_bm_longbone_1','fetal_bm_longbone_2','fetal_bm_longbone_3'))
  
  fetal_bm_list_raw=sapply(names(fetal_bm_list), function(x){ sRNA_qc_cal_func(proname = x,seurate_sRNA_raw = fetal_bm_list[[x]])})
  fetal_bm_list_raw[[2]]
  fetal_bm_list_raw[[5]]
  fetal_bm_list_raw[[8]]
  fetal_bm_list_raw[[11]]
  fetal_bm_list_raw[[14]]
  fetal_bm_list_raw[[17]]
  
  fetal_bm_list_qc=sapply(c(1,4,7,10,13,16), function(x){ qc_normalization_varfeature_func(fetal_bm_list_raw[[x]],qc_pare = c(200,6000,60000,10,60))})
  fetal_bm_list_qc= fetal_bm_list_qc[c(1,3,5,7,9,11)];
  
  rm(fetal_bm_list_raw,fetal_bm_list)
  fetal_bm_seu=merge(fetal_bm_list_qc[[1]],c(fetal_bm_list_qc[2:6]));rm(fetal_bm_list_qc)
  colnames(fetal_bm_seu@meta.data)=gsub(pattern='orig.ident',replacement='age',colnames(fetal_bm_seu@meta.data))
  colnames(fetal_bm_seu@meta.data)=gsub(pattern='old.ident',replacement='donor',colnames(fetal_bm_seu@meta.data))
  fetal_bm_seu$resource='GSE143753'
  
  fetal_bm_seu$stage_type='limbbud'
  fetal_bm_seu$stage_type[grep(x = fetal_bm_seu$donor,pattern = 'longbone')]='longbone'
  
  
  # --------------------------------------统一基因symbol信息 -------------------------------------------#
  symbol_addpostfix=rownames(fetal_bm_seu)[!rownames(fetal_bm_seu) %in% all_shared_ensembl_id_info$ref_symbol]
  symbol_addpostfix=data.frame(strsplit2(symbol_addpostfix,split = '-'))[,1] 
  symbol_addpostfix_inf=all_shared_ensembl_id_info[all_shared_ensembl_id_info$ref_symbol %in% symbol_addpostfix,]
  symbol_addpostfix_inf=symbol_addpostfix_inf[order(symbol_addpostfix_inf$ref_symbol),]
  symbol_addpostfix_inf$number=1:dim(symbol_addpostfix_inf)[1]
  du_symbol_addpostfix=symbol_addpostfix_inf[match(symbol_addpostfix,symbol_addpostfix_inf$ref_symbol)+1,'ref_symbol1']
  names(du_symbol_addpostfix)=paste(symbol_addpostfix,'1',sep = '-')
  
  uniq_symbols=rownames(fetal_bm_seu)[rownames(fetal_bm_seu) %in% all_shared_ensembl_id_info$ref_symbol]
  new_symbol_list=all_shared_ensembl_id_info[match(uniq_symbols,all_shared_ensembl_id_info$ref_symbol),'ref_symbol1']
  names(new_symbol_list)=uniq_symbols
  new_symbol_list=c(du_symbol_addpostfix,new_symbol_list)
  
  
  du_symbol=as.character(new_symbol_list)[duplicated(as.character(new_symbol_list))];length(du_symbol)
  tmp_assay=GetAssayData(fetal_bm_seu,assay = 'RNA',layer = 'counts')
  rownames(tmp_assay)=as.character(new_symbol_list[rownames(tmp_assay)])
  
  du_tmp_assay=tmp_assay[rownames(tmp_assay) %in% du_symbol,]
  tmp_assay=tmp_assay[!rownames(tmp_assay) %in% du_symbol,]
  
  du_tmp_assay=aggregate(du_tmp_assay,by=list(rownames(du_tmp_assay)), FUN=sum) 
  du_tmp_assay=column_to_rownames(du_tmp_assay,'Group.1')
  du_tmp_assay=as(as.matrix(du_tmp_assay),'dgCMatrix')
  
  tmp_assay=rbind(tmp_assay,du_tmp_assay)
  fetal_bm_seu=CreateSeuratObject(counts = tmp_assay,min.cells = 10,min.features = 200,meta.data =fetal_bm_seu@meta.data )
  rm(tmp_assay,du_tmp_assay);gc()
  
  #-------------------annotate the celltype of EBM altas and fetal_bm_seu mapping to FBM-------------------#
  
  
  fetal_bm_seu=NormalizeData(fetal_bm_seu) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% RunHarmony(group.by.vars='donor') %>% RunUMAP(dims=1:30,reduction='harmony')
  ElbowPlot(fetal_bm_seu,50)
  fetal_bm_seu=RunTSNE(fetal_bm_seu,reduction = 'harmony',dims = 1:30)
  (DimPlot(fetal_bm_seu,group.by = 'donor',cols = cols,reduction = 'tsne')+DimPlot(fetal_bm_seu,group.by = 'donor',cols = cols,reduction = 'umap'))/
    FeaturePlot(fetal_bm_seu,features = c('GYPA','HBE1'))
  
  
  anchor1=FindTransferAnchors(reference =FBM_altas_seu ,query =fetal_bm_seu,reference.reduction  = 'pca',dims = 1:20 )
  fetal_bm_seu=MapQuery(anchorset =anchor1,query = fetal_bm_seu,reference = FBM_altas_seu,reference.reduction = 'pca',reduction.model = 'umap1',refdata = list(celltype='subcelltype')) # 注意，command 保存的是umap1 信息，如果是umap2模块信息，会混淆
  DimPlot(fetal_bm_seu,reduction = 'ref.umap',group.by = 'predicted.celltype',cols = cols)+
    DimPlot(fetal_bm_seu,reduction = 'umap',group.by = 'predicted.celltype',cols = cols)
  
  rm(anchor1)
  
  colnames(fetal_bm_seu@meta.data)=gsub(pattern ='percent.mt' ,replacement = 'mito',colnames(fetal_bm_seu@meta.data))
  fetal_bm_seu$resource='GSE143753'
  fetal_bm_seu$method='3GEX'
  #fetal_bm_seu@meta.data=fetal_bm_seu@meta.data[,-grep(pattern = 'subcelltype.score',colnames(fetal_bm_seu@meta.data))]
  
  DimPlot(fetal_bm_seu,reduction = 'ref.umap',group.by = 'predicted.celltype',cols = cols)+
    FeaturePlot(fetal_bm_seu,features = c('HBE1','GYPA'),reduction = 'ref.umap')
  
  fetal_bm_seu=FindNeighbors(fetal_bm_seu,reduction = 'harmony',dims = 1:30)
  fetal_bm_seu=FindClusters(fetal_bm_seu,resolution = c(0.2,0.5))
  DimPlot(fetal_bm_seu,reduction = 'umap',group.by = c('RNA_snn_res.0.2','RNA_snn_res.0.5'),cols = cols,label = T)
  
  Idents(fetal_bm_seu)='RNA_snn_res.0.2'
  fetal_bm_seu_cluster_markers=FindAllMarkers(fetal_bm_seu,only.pos = T)
  gene_list=list()
  for(i in 0:11){
    tmp=fetal_bm_seu_cluster_markers[fetal_bm_seu_cluster_markers$cluster==i,];
    gene_list[[i+1]]=tmp[,'avg_log2FC'] ;names(gene_list[[i+1]])=tmp$gene
    gene_list[[i+1]]=sort(gene_list[[i+1]],decreasing = T)
  }
  
  library(clusterProfiler)
  geoGO_list=sapply(gene_list, function(x){gseGO(x,ont = 'BP',OrgDb = org.Hs.eg.db,keyType = 'SYMBOL')})
  names(geoGO_list)=as.character(0:11)
  saveRDS(geoGO_list,file = 'NRBC_BM_altas/res_data/BM_cluster_res0.2_geoGO_list.rds')
  
  top5_geoGO_df=data.frame()
  for (i in names(geoGO_list)) {
    tmp=geoGO_list[[i]][1:5,c('ID','Description','setSize','NES','p.adjust','core_enrichment')]
    tmp$cluster=i
    tmp$ratio=as.numeric(sapply(tmp$core_enrichment,function(x){length(strsplit(x,split = '/')[[1]])}))/tmp$setSize
    tmp=tmp[,c('Description','setSize','NES','p.adjust','cluster','ratio')]
    top5_geoGO_df=rbind(top5_geoGO_df,tmp)
  }
  
  top5_geoGO_df$Description=factor(top5_geoGO_df$Description,levels = unique(top5_geoGO_df$Description))
  top5_geoGO_df$cluster=factor(top5_geoGO_df$cluster,levels = as.character(0:12))
  ggplot(top5_geoGO_df,aes(x=cluster,y=Description ,size=ratio,color=NES))+geom_point()+theme_classic()+scale_color_gradient(low = 'gray',high = 'firebrick3')
  
  top5_fetal_bm_seu_cluster_markers=fetal_bm_seu_cluster_markers%>% group_by(cluster) %>% top_n(n = 5,wt = -log10(p_val_adj)) %>% top_n(n = 5,wt = avg_log2FC )
  DotPlot(fetal_bm_seu,group.by = 'RNA_snn_res.0.2',features = unique(top5_fetal_bm_seu_cluster_markers$gene))+RotatedAxis()
  
  EBM_gene_signature=c('MSX1', 'HOXA9', 'PRRX1', 'HOXC6', 'TWIST2', 'PDGFRA', 'RUNX2', 'OSR2', 
                       'NOV', 'SFRP2', 'SOX9', 'ACAN', 'SIX1', 'MYOG', 'CDH5', 'CD68', 'GYPA', 'FGF8', 'EPCAM', 'SOX10')
  
  
  DotPlot(fetal_bm_seu,group.by = 'RNA_snn_res.0.2',features =EBM_gene_signature,cols = c('gray','firebrick3') )+RotatedAxis()
  DotPlot(fetal_bm_seu,group.by = 'RNA_snn_res.0.5',features =EBM_gene_signature,cols = c('gray','firebrick3') )+RotatedAxis()
  table(fetal_bm_seu$predicted.celltype)
  DimPlot(subset(fetal_bm_seu,predicted.celltype %in% c("SKELETAL MUSCLE","FIBROBLASTS")),reduction = 'umap',group.by = 'predicted.celltype',cols = cols)
  
  
  #11:Schwan cell,  SOX10+
  # 10:10_0/1:Mac(CD68+)， 10_2:Meg
  #9: epithelial cells（EPCAM+）
  #8:Endo(CDH5+)
  #7:Ery(GYPA+)
  #6/1: perichondrial mesenchymal stromal cells(PMSC, OSR2+NOV+)
  #5:5_0/1:Myoprogenitor(11/9:Myoprogenitor,SIX1+，5_2:： Myocyte,SIX1+MYOG+)
  #3: 6/8,3_1: osteoprogenitors (RUNX2+),3_0: 分两个亚群:SFRP2+SOX9+ACAN+， chondroblasts,3_2： SOX9+ACAN+chondrocytes
  #0：osteo-chondrogenic progenitors (OCPs，PRRX1+SOX9lowPDGFRAhigh)
  #4/2: BMSC1/BMSC2（ PRRX1+PDGFRA+）
  
  Idents(fetal_bm_seu)='RNA_snn_res.0.2'
  fetal_bm_seu=FindSubCluster(fetal_bm_seu,cluster = '10',resolution = 0.2,graph.name = 'RNA_snn')
  DimPlot(fetal_bm_seu,reduction = 'umap',group.by = 'sub.cluster',cols = cols,label = T)
  Idents(fetal_bm_seu)='sub.cluster'
  fetal_bm_seu=FindSubCluster(fetal_bm_seu,cluster = '3',resolution = 0.2,graph.name = 'RNA_snn',subcluster.name ='sub.cluster2' )
  DimPlot(fetal_bm_seu,reduction = 'umap',group.by = 'sub.cluster2',cols = cols,label = T)
  Idents(fetal_bm_seu)='sub.cluster2'
  fetal_bm_seu=FindSubCluster(fetal_bm_seu,cluster = '5',resolution = 0.1,graph.name = 'RNA_snn',subcluster.name ='sub.cluster3' )
  DimPlot(fetal_bm_seu,reduction = 'umap',group.by = 'sub.cluster3',cols = cols,label = T)
  
  DimPlot( fetal_bm_seu,  reduction = 'umap',  group.by = c('RNA_snn_res.0.2', 'RNA_snn_res.0.5', 'sub.cluster3'),cols = cols,label = T)
  
  pheatmap(cor(data.frame(AverageExpression(fetal_bm_seu,group.by = 'RNA_snn_res.0.2')$RNA)),color = c('white','firebrick3'),display_numbers = T)
  pheatmap(cor(data.frame(AverageExpression(fetal_bm_seu,group.by = 'sub.cluster3')$RNA)),color = c('white','firebrick3'),display_numbers = T)
  
  Idents(fetal_bm_seu)='sub.cluster3'
  
  
  cluster_celltype_list=list('11'='Schwan cell','10_0'='MACROPHAGE','10_1'='MACROPHAGE','10_2'='MEGAKARYOCYTE',
                             '9'='EPITHELIUM','8'='ENDOTHELIUM','7'='Ery','6'='PMSC1','1'='PMSC2','5_0'='Myoprogenitor','5_1'='Myoprogenitor','5_2'='MYOCYTE',
                             '3_1'='Osteoprogenitors','3_0'='CHONDROBLAST','3_2'='CHONDROCYTE','0'='OCPs','4'='BMSC1','2'='BMSC2') # osteo-chondrogenic progenitors: OCPs
  
  fetal_bm_seu$celltype=as.character(cluster_celltype_list[fetal_bm_seu$sub.cluster3])               
  
  # 由于integated by harmony 与ref.umap结果有明显差异，再次使用cca进行整合
  if(F){
    fetal_bm_seu=JoinLayers(fetal_bm_seu)
    fetal_bm_seu[['RNA']]=split(fetal_bm_seu[['RNA']],fetal_bm_seu$donor)
    fetal_bm_seu <- IntegrateLayers(object = fetal_bm_seu, method = CCAIntegration,orig.reduction = "pca", new.reduction = "integrated.cca",verbose = FALSE)
    fetal_bm_seu=RunUMAP(fetal_bm_seu,reduction = 'integrated.cca',dims = 1:30,reduction.name = 'cca_umap')
    DimPlot(fetal_bm_seu,reduction = 'cca_umap',group.by =c( 'donor','predicted.celltype','celltype'),cols = cols)
    
  }
  
  fetal_bm_seu=FindSubCluster(fetal_bm_seu,cluster = '7',graph.name = 'RNA_snn',subcluster.name = 'Ery_subcluster')
  # EBM中的Ery主要是来自limbbud， 基本都是Bas，处于G1期，5WPC，13,15WPC 
  DimPlot(subset(fetal_bm_seu,RNA_snn_res.0.2=='7'),reduction = 'umap',group.by = c('Ery_subcluster','predicted.celltype'),cols = cols,label = T) 
  
  fetal_bm_Ery_seu=subset(fetal_bm_seu,RNA_snn_res.0.2=='7')
  fetal_bm_Ery_seu=RunUMAP(fetal_bm_Ery_seu,reduction = 'harmony',reduction.name = 'sub_umap',dims = 1:10)
  
  DimPlot(fetal_bm_Ery_seu,reduction = 'sub_umap',group.by = c('Ery_subcluster','stage_type','Phase','predicted.celltype'),cols = cols,label = T) 
  DimPlot(fetal_bm_Ery_seu,reduction = 'ref.umap',group.by = c('Ery_subcluster','stage_type','Phase','predicted.celltype'),cols = cols,label = T) 
  
  
  celltype_features=c('KIT',"TFRC","GYPA",'HBE1',"MALAT1",'NCL','CD63','MKI67','CD74')
  FeaturePlot(fetal_bm_Ery_seu,features = celltype_features,reduction = 'sub_umap')
  FeaturePlot(fetal_bm_Ery_seu,features = celltype_features,reduction = 'ref.umap')
  
  VlnPlot(fetal_bm_Ery_seu,features = celltype_features,group.by = c('Ery_subcluster'),stack = T)
  
  fetal_bm_Ery_seu$celltype='YS_Bas'
  fetal_bm_seu$celltype[fetal_bm_seu$celltype=='Ery']='YS_Bas'
  
  p=DimPlot(fetal_bm_seu,group.by = c('celltype','sub.cluster3'),reduction ='umap',cols = col,label = T )
  p
  ggsave(p,file='NRBC_BM_altas/res_pic/EBM_altas_celltype_subcluster_umap.pdf',width = 12,height = 6)
  
  rm(fetal_bm_seu_counts);gc()
  
  saveRDS(fetal_bm_seu,file = 'NRBC_BM_altas/tmp_embryo_bm_seu.rds')
}else{
  fetal_bm_seu=readRDS('NRBC_BM_altas/tmp_embryo_bm_seu.rds')
}

table(rownames(fetal_bm_seu) %in% all_shared_ensembl_id_info$ref_symbol1)# FALSE/TRUE : 6/20998
table(rownames(fetal_bm_seu) %in% all_shared_ensembl_id_info$ref_symbol)# FALSE/TRUE :3838/17166


 
################################################################################################################# 
#--------------------------------------------------ABM ALTAS DATA--------------------------------------------#
################################################################################################################# 

if(F){

  # 查看disco_bone_marrow_v01.rds数据是否满足要求,该数据集，是由不同数据整合而成，
  #缺乏MSC等基质细胞，部分数据来源是BM 单核细胞；舍弃该数据作为ABM altas数据
  
  run_time=2
  if(run_time==0){
  BM_seu=readRDS('ref_data/ref_scRNAseq_data/human_cell_atlas/disco_bone_marrow_v01.rds')
  DimPlot(BM_seu,group.by = 'ct',label = T)
  BM_seu$celltype1=BM_seu$ct
  unique(BM_seu$celltype1)
  BM_seu$celltype1[grep(pattern = 'monocyte',BM_seu$celltype1)]='monocyte'
  BM_seu$celltype1[grep(pattern = 'NK',BM_seu$celltype1)]='NK'
  BM_seu$celltype1[grep(pattern = 'B cell',BM_seu$celltype1)]='B cell'
  BM_seu$celltype1[grep(pattern = 'Immature B',BM_seu$celltype1)]='B cell'
  BM_seu$celltype1[BM_seu$celltype1 %in% c('Naive B','Memory B')]='B cell'
  BM_seu$celltype1[BM_seu$celltype1 %in% c('pDC','cDC')]='DC'
  
  
  BM_seu$celltype1[grep(pattern = 'CD8 T',BM_seu$celltype1)]='CD8 T'
  BM_seu$celltype1[grep(pattern = 'CD4 T',BM_seu$celltype1)]='CD4 T'
  #BM_seu$celltype1[grep(pattern = 'Treg',BM_seu$celltype1)]='T cell'
  
  BM_seu$celltype1[grep(pattern = "Megakaryocyte/erythroid progenitor",BM_seu$celltype1)]='MEP'
  BM_seu$celltype1[grep(pattern = "Common lymphoid progenitor",BM_seu$celltype1)]='CLP'
  BM_seu$celltype1[grep(pattern = "Proliferation HSC",BM_seu$celltype1)]='HSC'
  BM_seu$celltype1[grep(pattern = "Granulocyte progenitor",BM_seu$celltype1)]='Granulocyte'
  BM_seu$celltype1[grep(pattern =  "Megakaryocyte progenitor",BM_seu$celltype1)]='Megakaryocyte'
  
  unique(BM_seu$celltype1)
  DimPlot(BM_seu,group.by = 'celltype1',label = T,raster = F)
  
  # BM_seu存在data矩阵为空，导致细胞id对应不上，重新创见对象
  BM_seu_new=CreateSeuratObject(counts = BM_seu@assays$RNA@counts,assay = 'RNA',meta.data =BM_seu@meta.data,project = 'BM_aging')
  BM_seu_new[['pca']]=CreateDimReducObject(BM_seu@reductions$pca@cell.embeddings,loadings = BM_seu@reductions$pca@feature.loadings,stdev = BM_seu@reductions$pca@stdev)
  BM_seu_new[['umap']]=CreateDimReducObject(BM_seu@reductions$umap@cell.embeddings)
  rm(BM_seu)
  DimPlot(BM_seu_new,group.by = 'celltype1',label = T)
  
  saveRDS(BM_seu_new,file = 'NRBC_BM_altas/disco_BM_seu_new.rds')
  run_time=run_time+1
}else if(run_time==1 )
  {
  
  BM_seu_new= readRDS('NRBC_BM_altas/disco_BM_seu_new.rds')
  table(BM_seu_new$projectId)
  BM_seu_new=NormalizeData(BM_seu_new)
  table(subset(BM_seu_new,celltype1 %in% c("MEP","Proerythroblast" ,"Erythroblast" ,"RBC" ))$projectId )
  # E-MTAB-9139   GSE120221   GSE133181   GSE135194   GSE159624   GSE159929   GSE165645   GSE169426   GSE175604   GSE179346   GSE181989   GSE188222   GSE193138 
  # 13             6040        1546        3513          11          43         323          980           9           7         646         223         429
  # 排除E-MTAB-9139、GSE159624、GSE159929、 GSE175604   GSE179346
  
  VlnPlot(subset(BM_seu_new,celltype1 %in% c("MEP","Proerythroblast" ,"Erythroblast" ,"RBC" ) ),features = c('IGKC','IGHG1','IGHA1','TFRC','CA1'),group.by = 'projectId',stack=T) # BM_seu_new
  # 排除"GSE120221"   "GSE193138"  "GSE179346"    "GSE188222" ，高表达IGH/IGK基因
  
  sub_BM_seu_new=subset(BM_seu_new, projectId %in% c('GSE133181','GSE135194','GSE165645','GSE169426','GSE181989')) # GSE120221 样本抽取6个样本，数量应该差不多,
  
  #-------------------------check the bm altas project info------------------#
  {
    # E-MTAB-9139,都是60岁以上的老人，免疫细胞,不可用
    
    #GSE120221： 10X Genomics Single Cell 3’ Solution, version 2 ,mononuclear cells，可用,但是红细胞数据比较少
    # paper: Human bone marrow assessment by single-cell RNA sequencing, mass cytometry, and flow cytometry
    # GSE159624, NK cells，不可用
    # Single-cell profiling reveals the trajectories of natural killer cell differentiation in bone marrow and 
    # a stress signature induced by acute myeloid leukemia. Cell Mol Immunol 2021 May;18(5):1290-1304. PMID: 33239726
    
    #GSE159929, hg38, version 3.0.1; 10x Genomics,5'-end, 以T细胞为主，NRBC少
    # Single-cell transcriptome profiling of an adult human cell atlas of 15 major organs. Genome Biol 2020 Dec 7;21(1):294. PMID: 33287869
    # GSE175604 ： t cells，不可用
    # Luo Y, Xu C, Wang B, Niu Q et al. Single-cell transcriptomic analysis reveals disparate effector differentiation pathways
    #in human T(reg) compartment. Nat Commun 2021 Jun 23;12(1):3913. PMID: 34162888
    
    #GSE179346： CD45+ bone marrow cells,不可
    #GSE188222，bone marrow plasma cells，不可
    
    #GSE193138，remove erythrocytes and were then labeled with antibodies. 可用
    #Through flow sorting, the final cells used for single-cell RNA sequencing (scRNA-seq) contained 
    #CD45+CD66b- population (90%), CD66b+ population (5%), and CD45- population (5%).
    #paper ；Decoding lymphomyeloid divergence and immune hyporesponsiveness in G-CSF-primed human bone marrow by single-cell RNA-seq. 
    # Cell Discov 2022 Jun 22;8(1):59. PMID: 35732626
    
    #GSE133181:  CD34+ cells,32717, 5 samples, no age info，可选
    #Single-cell analysis of bone marrow-derived CD34+ cells from children with sickle cell disease and thalassemia.
    #Blood 2019 Dec 5;134(23):2111-2115. PMID: 31697810
    DimPlot(subset(BM_seu_new,projectId=='GSE133181'),group.by = 'celltype1',cols = cols)
    
    #GSE135194:Lineage(CD3CD14CD19)-CD34+,4 samples,14757
    #the hg19 genome with STAR ,paper: Single-cell RNA-seq reveals a distinct transcriptome signature of hematopoiesis in GATA2 deficiency
    #Genomics Chromium Single Cell ‘3 Solution Kit V2 from 10XGenomics.
    #, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135194
    DimPlot(subset(BM_seu_new,projectId=='GSE133181'),group.by = 'celltype1',cols = cols)
    
    
    # GSE165645, Mononuclear cells , the 10x Genomics Chromium Single Cell 5’ Library and Gel Bead kit v1 ,GRCh38，可用, 但是NRBC少
    # paper: Le Coz C, Nguyen DN, Su C, Nolan BE et al. Constrained chromatin accessibility in
    # PU.1-mutated agammaglobulinemia patients. J Exp Med 2021 Jul 5;218(7). PMID: 33951726
    DimPlot(subset(BM_seu_new,projectId=='GSE165645'),group.by = 'celltype1',cols = cols)
    
    # GSE169426， CD34+ cells 
    # paper：Hematopoietic differentiation is characterized by a transient peak of entropy at a single-cell level. BMC Biol 2022 Mar 9;20(1):60. PMID: 35260165
    DimPlot(subset(BM_seu_new,projectId=='GSE169426'),group.by = 'celltype1',cols = cols)
    
    # GSE181989： MNCs(Mononuclear cells ) mixed with CD34+ cells at 4:1 ratio and were analyzed by 10×Genomics.候选
    #： ，Single-Cell RNA-Seq of Bone Marrow Cells in Aplastic Anemia. Front Genet 2021;12:745483. PMID: 35046994
    DimPlot(subset(BM_seu_new,projectId=='GSE181989'),group.by = 'celltype1',cols = cols)
    
  }
  
   
   
  table(sub_BM_seu_new$projectId)
  sub_BM_seu_new$age[is.na(sub_BM_seu_new$age)]='unknown'
  dim(sub_BM_seu_new) # 64846
  
  DimPlot(sub_BM_seu_new,group.by = c('celltype1','projectId'),raster = F,cols = cols) 
  DimPlot(sub_BM_seu_new,group.by = c('ct','seurat_clusters'),label = T,raster = F,cols = cols)
  
  
  run_time=run_time+1
  rm(BM_seu_new);gc()
  
  #Proerythroblast: BFUE/CFUE,ENG:CD105, 红系细胞方向分化标志，ITGA2B：CD41，Megakaryocyte,标志，MEP就是MEP，不考虑其中的BFUE，MEP：CD36 低表达
  VlnPlot(subset(sub_BM_seu_new,celltype1 %in% c("MEP","Proerythroblast" ,"Erythroblast" ,"RBC" )),group.by = 'celltype1',features = c('CD34','KIT','CD36','TFRC','ENG','ITGA2B','GYPA'),stack = T) # ITGA2B:Megakaryocyte
  table(subset(sub_BM_seu_new,celltype1 %in% c("Proerythroblast" ,"Erythroblast" ,"RBC" ))$projectId )
  #GSE133181 GSE135194 GSE165645 GSE169426 GSE181989 
  #   806      2338       280     488        348 
  table(sub_BM_seu_new@meta.data[,c('projectId','sample')])
  
  sub_BM_seu_new=NormalizeData(sub_BM_seu_new) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData() %>% RunPCA() 
  anchor=FindTransferAnchors(reference =FBM_altas_seu,query = sub_BM_seu_new,reference.reduction = 'pca',dims = 1:30 )
  sub_BM_seu_new=MapQuery(anchorset =anchor,reference =FBM_altas_seu,query = sub_BM_seu_new, reference.reduction = "pca", reduction.model = 'umap1',refdata = list(celltype='subcelltype'))
  
  DimPlot(sub_BM_seu_new,group.by = c('ct','predicted.celltype'),cols = cols,reduction = 'umap')
  
  DimPlot(subset(subset(sub_BM_seu_new, celltype1 %in% c( 'MEP',"Proerythroblast" , "Erythroblast","RBC" )),predicted.celltype %in% c('BFUE/CFUE','ProE','Bas','Poly','Orth')),group.by = c('celltype1','predicted.celltype'),cols = cols,reduction = 'umap')
  
  sub_BM_seu_new$celltype1[sub_BM_seu_new$predicted.celltype %in% c('BFUE/CFUE','ProE','Bas','Poly','Orth')]=sub_BM_seu_new$predicted.celltype[sub_BM_seu_new$predicted.celltype %in% c('BFUE/CFUE','ProE','Bas','Poly','Orth')]
  table(sub_BM_seu_new$predicted.celltype[sub_BM_seu_new$celltype1=='Proerythroblast'])
  sub_BM_seu_new$celltype1[sub_BM_seu_new$celltype1=='Proerythroblast']='BFUE/CFUE'
  DimPlot(sub_BM_seu_new,group.by = c('celltype1','predicted.celltype'),cols = cols,reduction = 'umap')
  
  sub_BM_seu_new$subcelltype=sub_BM_seu_new$celltype1
  
  saveRDS(sub_BM_seu_new,'NRBC_BM_altas/tmp_sub_ABM_seu_new.rds')
  
}else{
    sub_BM_seu_new=readRDS('NRBC_BM_altas/tmp_sub_ABM_seu_new.rds')
    sub_BM_seu_new$resource=sub_BM_seu_new$projectId
  }
  
  # check gene symbol ananotation
  list_seu=list(GSE133181=subset(sub_BM_seu_new,projectId=='GSE133181'),GSE135194=subset(sub_BM_seu_new,projectId=='GSE135194'),
                GSE165645=subset(sub_BM_seu_new,projectId=='GSE165645'),GSE169426=subset(sub_BM_seu_new,projectId=='GSE169426'),
                GSE181989=subset(sub_BM_seu_new,projectId=='GSE181989')
                )
  list_seu=sapply(list_seu,function(x){subset(x, features = rownames(x)[Matrix::rowSums(GetAssayData(x, layer = "counts")) > 0])})
  
  
  #GSE135194: Genome_build: hg19, 原文中
  # GSE165645 Citeseq data
  # 文章已经对数据 gene symbol 进行了统一处理处，得的的symbol 同ourlab gene symbol 信息一致,nice 
  for (cho_ncol in colnames(all_shared_ensembl_id_info)) { print(cho_ncol); print(table(rownames(list_seu[['GSE133181']]) %in% all_shared_ensembl_id_info[,cho_ncol]))} #ref_symbol, zxlab_symbol
  for (cho_ncol in colnames(all_shared_ensembl_id_info)) { print(cho_ncol); print(table(rownames(list_seu[['GSE135194']]) %in% all_shared_ensembl_id_info[,cho_ncol]))} #ref_symbol,zxlab_symbol
  for (cho_ncol in colnames(all_shared_ensembl_id_info)) { print(cho_ncol); print(table(rownames(list_seu[['GSE165645']]) %in% all_shared_ensembl_id_info[,cho_ncol]))} #ref_symbol,zxlab_symbol
  for (cho_ncol in colnames(all_shared_ensembl_id_info)) { print(cho_ncol); print(table(rownames(list_seu[['GSE169426']]) %in% all_shared_ensembl_id_info[,cho_ncol]))} #ref_symbol,zxlab_symbol
  for (cho_ncol in colnames(all_shared_ensembl_id_info)) { print(cho_ncol); print(table(rownames(list_seu[['GSE181989']]) %in% all_shared_ensembl_id_info[,cho_ncol]))} #ref_symbol,zxlab_symbol
  
  for (cho_ncol in colnames(all_shared_ensembl_id_info)) { print(cho_ncol); print(table(rownames(sub_BM_seu_new) %in% all_shared_ensembl_id_info[,cho_ncol]))} #ref_symbol, zxlab_symbol
  
  # 原文件中gene 注释信息, 不统一
  if(F){
  GSE133181_gene_info=read.table('ref_data/gene_symbol_annotation/GSE133181_genes.tsv',sep="\t");rownames(GSE133181_gene_info)=GSE133181_gene_info$V1
  GSE165645_gene_info=read.table('ref_data/gene_symbol_annotation/GSE165645_features.tsv',sep="\t");rownames(GSE165645_gene_info)=GSE165645_gene_info$V1
  GSE169426_gene_info=read.table('ref_data/gene_symbol_annotation/GSE169426_genes.tsv',sep="\t");rownames(GSE169426_gene_info)=GSE169426_gene_info$V1
  GSE181989_gene_info=read.table('ref_data/gene_symbol_annotation/GSE181989_features.tsv',sep="\t");rownames(GSE181989_gene_info)=GSE181989_gene_info$V1
  
  GSE133181_gene_info$ref_symbol1=all_shared_ensembl_id_info[rownames(GSE133181_gene_info),'ref_symbol1']
  GSE165645_gene_info$ref_symbol1=all_shared_ensembl_id_info[rownames(GSE165645_gene_info),'ref_symbol1']
  GSE169426_gene_info$ref_symbol1=all_shared_ensembl_id_info[rownames(GSE169426_gene_info),'ref_symbol1']
  GSE181989_gene_info$ref_symbol1=all_shared_ensembl_id_info[rownames(GSE181989_gene_info),'ref_symbol1']
  
}
  
  
  # Putative regulators for the continum of erythroid differentiation revealed by single-cell transcriptome
  if(F){
  # create the ABM SeuratObject data
  # ratio：ABM altas cell
  Cell_GSE253355_ABM_seu=readRDS('ref_data/ref_scRNAseq_data/GSE253355_ABM/GSE253355_Normal_Bone_Marrow_Atlas_Seurat_SB_v2.rds')
  table(rownames(FBM_altas_seu) %in% rownames(Cell_GSE253355_ABM_seu)) # FALSE  TRUE ,4106 29432
  
  Cell_GSE253355_ABM_seu=NormalizeData(Cell_GSE253355_ABM_seu) %>% FindVariableFeatures(nfeatures =3000) %>% ScaleData() %>% RunPCA()
  anchor=FindTransferAnchors(reference =FBM_altas_seu ,query =Cell_GSE253355_ABM_seu,reference.reduction  = 'pca' )
  Cell_GSE253355_ABM_seu=MapQuery(anchorset =anchor,query = Cell_GSE253355_ABM_seu,reference = FBM_altas_seu,reference.reduction = 'pca',reduction.model = 'umap1',refdata = list(celltype='subcelltype'))
  ## reference.reduction = 'pca'或者harmony对结果几乎没有影响
  DimPlot(Cell_GSE253355_ABM_seu,reduction = 'ref.umap',group.by = c('predicted.celltype','cluster_anno_l2'),cols = cols) 
  # mapping效果比较差，可能是symbol没有矫正有关
  DimPlot(Cell_GSE253355_ABM_seu,reduction = 'UMAP_dim30',group.by = c('predicted.celltype','cluster_anno_l2'),cols = cols)
  
  Cell_GSE253355_ABM_NRBC_seu=subset(Cell_GSE253355_ABM_seu,cluster_anno_l2 %in% c('MEP','Megakaryocyte','Late Erythroid','Erythroblast','RBC'))
  DimPlot(Cell_GSE253355_ABM_NRBC_seu,reduction = 'UMAP_dim30',group.by = c('predicted.celltype','cluster_anno_l2'),cols = cols)
  Cell_GSE253355_ABM_NRBC_seu=RunPCA(Cell_GSE253355_ABM_NRBC_seu)
  Cell_GSE253355_ABM_NRBC_seu=RunUMAP(Cell_GSE253355_ABM_NRBC_seu,dims = 1:30,reduction.name = 'umap1')
  DimPlot(Cell_GSE253355_ABM_NRBC_seu,reduction = 'umap1',group.by = c('predicted.celltype','cluster_anno_l2'),cols = cols)
  Cell_GSE253355_ABM_NRBC_seu=FindClusters(Cell_GSE253355_ABM_NRBC_seu,resolution = 0.4)
  Idents(Cell_GSE253355_ABM_NRBC_seu)='seurat_clusters'
  Cell_GSE253355_ABM_NRBC_seu=FindSubCluster(Cell_GSE253355_ABM_NRBC_seu,resolution = 0.3,cluster = '4',graph.name ='RNA_snn' )
  
  DimPlot(Cell_GSE253355_ABM_NRBC_seu,reduction = 'umap1',group.by = c('sub.cluster','cluster_anno_l2'),cols = cols,label = T)/  
    FeaturePlot(Cell_GSE253355_ABM_NRBC_seu,features = c('CD34','PF4','HBB'),reduction = 'umap1',cols = c('gray','firebrick3'),ncol = 3)/
    VlnPlot(Cell_GSE253355_ABM_NRBC_seu,features =c('CD34','GYPA','PF4','HBB'),stack = T,group.by = 'sub.cluster' )
  
  pheatmap(cor(data.frame(AverageExpression(Cell_GSE253355_ABM_NRBC_seu,group.by = 'sub.cluster')$RNA)))
  
  # defined the Ery and Mek
  # 特点：MEP：CD34+HBB-PF-;BFUE/CFUE: CD34+HBBlowGYPA-,ProE:CD34-GPYAlowHBB, Bas:CD34-GYPAmidHBBmid, Poly:CD34-GYPAmidHBBmid,PolyNCLlow, Orth:GYPAhighHBBhighNCL-
  Cell_GSE253355_ABM__NRBC_celltype_list=list('3'='MEP','1'='BFUE/CFUE','4_1'='BFUE/CFUE','4_2'='BFUE/CFUE','6'='BFUE/CFUE','4_0'='Meg','4_3'='Meg','5'='ProE','0'='Bas','2'='Bas','8'='Poly','7'='Orth')
  Cell_GSE253355_ABM_NRBC_seu$subcelltype=as.character(Cell_GSE253355_ABM_NRBC_celltype_list[Cell_GSE253355_ABM__NRBC_seu$sub.cluster])
  DimPlot(Cell_GSE253355_ABM_NRBC_seu,reduction = 'umap1',group.by = c('subcelltype'),cols = cols,label = T)
  
  
  Cell_GSE253355_ABM_seu$subcelltype=as.character(Cell_GSE253355_ABM_seu$cluster_anno_l2)
  Cell_GSE253355_ABM_seu@meta.data[rownames(Cell_GSE253355_ABM_NRBC_seu@meta.data),'subcelltype']=as.character(Cell_GSE253355_ABM_NRBC_seu$subcelltype)
  DimPlot(Cell_GSE253355_ABM_seu,reduction = 'UMAP_dim30',group.by = c('subcelltype'),cols = cols,label = T,label.size = )+NoLegend()+
    DimPlot(Cell_GSE253355_ABM_seu,reduction = 'ref.umap',group.by = c('subcelltype'),cols = cols,label = T,label.size = 2)
  
  # 血管平滑肌细胞(vascular smooth muscle cell,VSMC)是构成血管中膜的主要细胞成分
  # Endo:'PECAM1',  ACTA2:SVAM, 缺乏红细胞晚期细胞类型,补充部分红细胞
  FeaturePlot(Cell_GSE253355_ABM_seu,features = c('PECAM1','ACTA2'))
  saveRDS(Cell_GSE253355_ABM_seu,file = 'NRBC_BM_altas/tmp_Cell_GSE253355_ABM_seu.rds')
    
}else{
      Cell_GSE253355_ABM_seu=readRDS( 'NRBC_BM_altas/tmp_Cell_GSE253355_ABM_seu.rds')
  }
    
  #------------------------ of human BM and UCB cells,GSE150774---------------------------------#
  # A 10-mL BM sample was collected from six young healthy male donors, health young man
  if(F){
    # CD34- CD235+
    GSE150774_BM1_seu=sRNA_qc_cal_func(inpath = '/home/gibh/2021_NRBC_chlyu/ref_data/ref_scRNAseq_data/GSE150774/HP2_filtered_gene_bc_matrices',proname = 'ABM1')
    GSE150774_BM2_seu=sRNA_qc_cal_func(inpath = '/home/gibh/2021_NRBC_chlyu/ref_data/ref_scRNAseq_data/GSE150774/HP3_filtered_gene_bc_matrices',proname = 'ABM2')
    GSE150774_BM3_seu=sRNA_qc_cal_func(inpath = '/home/gibh/2021_NRBC_chlyu/ref_data/ref_scRNAseq_data/GSE150774/HP4_filtered_gene_bc_matrices',proname = 'ABM3')
    GSE150774_BM1_seu[[2]]
    GSE150774_BM2_seu[[2]]
    GSE150774_BM3_seu[[2]]
    
    GSE150774_BM_seu=merge(GSE150774_BM1_seu[[1]],c(GSE150774_BM2_seu[[1]],GSE150774_BM3_seu[[1]]))
    GSE150774_BM_seu[['RNA']]=JoinLayers(GSE150774_BM_seu[['RNA']])
    rm(GSE150774_BM1_seu,GSE150774_BM2_seu,GSE150774_BM3_seu)
    
    qc_pare = c(200,2000,30000,10,60)
    GSE150774_BM_seu=subset(GSE150774_BM_seu,subset = nFeature_RNA > qc_pare[1] & nFeature_RNA < qc_pare[2] & nCount_RNA <qc_pare[3] & percent.mt < qc_pare[4] & percent.rb < qc_pare[5])
    colnames(GSE150774_BM_seu@meta.data)=gsub(pattern = 'orig.ident',replacement = 'donor',colnames(GSE150774_BM_seu@meta.data))
    colnames(GSE150774_BM_seu@meta.data)=gsub(pattern = 'percent.mt',replacement = 'mito',colnames(GSE150774_BM_seu@meta.data))
    dim(GSE150774_BM_seu) # 12352
    
    GSE150774_BM_seu=NormalizeData(GSE150774_BM_seu) 
    GSE150774_BM_seu=FindVariableFeatures(GSE150774_BM_seu,nfeatures = 2000) %>% ScaleData() %>% RunPCA() 
    GSE150774_BM_seu= RunHarmony(GSE150774_BM_seu,group.by.vars='donor') %>% RunUMAP(dims=1:30,reduction='harmony')
    DimPlot(GSE150774_BM_seu,group.by = 'donor')+FeaturePlot(GSE150774_BM_seu,features = 'GYPA')
    
    ElbowPlot(GSE150774_BM_seu)
    GSE150774_BM_seu=RunUMAP(GSE150774_BM_seu,dims=1:8,reduction='harmony')
    FeaturePlot(GSE150774_BM_seu,features =c('CD34','GYPA','HBB'))
    
    GSE150774_BM_seu=RunUMAP(GSE150774_BM_seu,dims=1:5,reduction='harmony',reduction.name = 'umap1')
    FeaturePlot(GSE150774_BM_seu,features =c('CD34','GYPA','HBB'),reduction = 'umap1')
    
    
    GSE150774_BM_seu=FindNeighbors(GSE150774_BM_seu,reduction = 'harmony',dims = 1:20)
    GSE150774_BM_seu=FindClusters(GSE150774_BM_seu,resolution = 0.4)
    DimPlot(GSE150774_BM_seu,cols = col)
    
    load('/home/gibh/2021_NRBC_chlyu/ref_data/NRBC_ref_bullk_RNAseq_se.Rdata')
    GSE150774_BM_seu=singleR_analysis_func(refdata =nrbc_ref_se2,test =as.matrix(GSE150774_BM_seu@assays$RNA$data)[VariableFeatures(ABM_NRBC_seu),],outdata = ABM_NRBC_seu,an_type1 ='celltype', an_type2 = 'predict.lable1')
    GSE150774_BM_seu=GSE150774_BM_seu[[1]]
    
    DimPlot(GSE150774_BM_seu,group.by = c('donor','predict.lable1'),cols = cols)
    
    FBM_altas_seu=readRDS('NRBC_BM_altas/FBM_altas_seu.rds')
    table(rownames(FBM_altas_seu) %in% rownames(GSE150774_BM_seu))# FALSE  TRUE , 21039 12499
    
    all_shared_ensembl_id_info=read.csv('ref_data/all_shared_ensembl_id_info.csv',header=T)
    all_shared_ensembl_id_info=all_shared_ensembl_id_info[,-1:-2]
    table(rownames(GSE150774_BM_seu) %in% all_shared_ensembl_id_info$abm_symbol)
    reedited_symbol=rownames(GSE150774_BM_seu)[!rownames(GSE150774_BM_seu) %in% all_shared_ensembl_id_info$abm_symbol];reedited_symbol
    
    anchor=FindTransferAnchors(reference =FBM_altas_seu,query = GSE150774_BM_seu,reference.reduction = 'pca',dims = 1:30 )
    GSE150774_BM_seu=MapQuery(anchorset =anchor,reference =FBM_altas_seu,query = GSE150774_BM_seu, reference.reduction = "pca", reduction.model = 'umap1',refdata = list(celltype='subcelltype'))
    DimPlot(GSE150774_BM_seu,group.by = c('predict.lable1','predicted.celltype',''),cols = cols,reduction = 'ref.umap')/
      DimPlot(GSE150774_BM_seu,group.by = c('predict.lable1','predicted.celltype',''),cols = cols,reduction = 'umap')
    
    GSE150774_BM_seu=FindNeighbors(GSE150774_BM_seu,reduction = 'harmony',dims = 1:20)
    GSE150774_BM_seu=FindClusters(GSE150774_BM_seu,resolution =c(0.2))
    DimPlot(GSE150774_BM_seu,group.by = c('predict.lable1','predicted.celltype','RNA_snn_res.0.2','RNA_snn_res.0.4'),cols = cols,reduction = 'umap',label = T)
    GSE150774_BM_seu=FindSubCluster(GSE150774_BM_seu,cluster = '1',resolution = 0.2,graph.name = 'RNA_snn')
    DimPlot(GSE150774_BM_seu,group.by = c('predict.lable1','predicted.celltype'),cols = cols,reduction = 'umap')/
      DimPlot(GSE150774_BM_seu,group.by = c('RNA_snn_res.0.2','sub.cluster'),cols = cols,reduction = 'umap',label = T)+
      VlnPlot(GSE150774_BM_seu,features = c('CD34','TFRC','NCL','GYPA','CD3D','CD79A','CD14','CD163'),group.by = 'sub.cluster',stack = T)
    # HSPC:CD34,NRBC:TFRC,NCL,HBB,GYPA, T cell :CD3D, B CELL: CD79A, Mono: CD14, Mac: CD163,
    
    table(GSE150774_BM_seu@meta.data[,c('seurat_clusters','predicted.celltype')])
    VlnPlot(subset(GSE150774_BM_seu,seurat_clusters=='8'),features = c('TFRC','NCL','GYPA','CD79A','CD14','CD163','CD68'),group.by = 'predicted.celltype',stack = T)
    
    pheatmap(cor(data.frame(AverageExpression(GSE150774_BM_seu,group.by = 'sub.cluster')$RNA)))
    
    cluster_list=list('0'='Orth','2'='Orth','3'='Poly','1_1'='Poly','1_2'='Poly','1_0'='Bas','1_3'='Bas','1_4'='ProE',
                      '7'='B CELLs','4'='NK/T CELLS','5'='NK/T CELLS','6'='MONOCYTES','8'='Mac_Ery') # cluster 8 或可以剔除
    
    
    GSE150774_BM_seu$subcelltype=as.character(cluster_list[GSE150774_BM_seu$sub.cluster])
    
    
    # 发现存在部分非红细胞，免疫细胞，monocyte, NKT cells，B cell
    saveRDS(GSE150774_BM_seu,file = 'NRBC_BM_altas/tmp_GSE150774_BM_seu.rds')
}else{
    GSE150774_BM_seu=readRDS('NRBC_BM_altas/tmp_GSE150774_BM_seu.rds')
  }
    
  table(rownames(FBM_altas_seu) %in% rownames(Cell_GSE253355_ABM_seu)) # FALSE  TRUE ,4106 29432
  table(rownames(GSE150774_BM_seu) %in% all_shared_ensembl_id_info$ref_symbol1) # FALSE  TRUE , 5072/17809
  table(rownames(GSE150774_BM_seu) %in% all_shared_ensembl_id_info$ref_symbol) # FALSE  TRUE ,1804/12806
  table(rownames(Cell_GSE253355_ABM_seu) %in% all_shared_ensembl_id_info$ref_symbol1) # FALSE  TRUE ,8754/20680 
  table(rownames(Cell_GSE253355_ABM_seu) %in% all_shared_ensembl_id_info$ref_symbol) # FALSE  TRUE 20/29434
  
  
  # ------------------------------统一meta.data 以及symbol信息-----------------------------------#
  cho_colnames=c('orig.ident','nCount_RNA' ,'nFeature_RNA','percent.mt','cluster_anno_l2','cluster_anno_l1', 'Sex','Age','subcelltype')
  Cell_GSE253355_ABM_seu@meta.data=Cell_GSE253355_ABM_seu@meta.data[,cho_colnames]
  Cell_GSE253355_ABM_seu$platform='10xv3'
  Cell_GSE253355_ABM_seu$resource='GSE253355'
  Cell_GSE253355_ABM_seu$sample=Cell_GSE253355_ABM_seu$orig.ident
  colnames(Cell_GSE253355_ABM_seu@meta.data)=gsub(pattern = 'percent.mt',replacement = 'mito',colnames(Cell_GSE253355_ABM_seu@meta.data))
  colnames(Cell_GSE253355_ABM_seu@meta.data)=gsub(pattern = 'cluster_anno_l2',replacement = 'celltype',colnames(Cell_GSE253355_ABM_seu@meta.data))
  
  GSE150774_BM_seu$resource='GSE150774'
  GSE150774_BM_seu$sample=GSE150774_BM_seu$donor
  GSE150774_BM_seu$platform='10x3v2'
  
  table(rownames(Cell_GSE253355_ABM_seu) %in% rownames(GSE150774_BM_seu)) #FALSE  TRUE,  16972 12480 
  
  
  GSE150774_gene_info=read.table('ref_data/ref_scRNAseq_data/GSE150774/HP2_filtered_gene_bc_matrices/genes.tsv')
  rownames(GSE150774_gene_info)=GSE150774_gene_info$V1
  GSE253355_gene_info=read.table(gzfile('ref_data/ref_scRNAseq_data/GSE253355_ABM/GSM8019212_H14_MACS_features.tsv.gz'),sep="\t")
  rownames(GSE253355_gene_info)=GSE253355_gene_info$V1;GSE253355_gene_info=GSE253355_gene_info[,-1]
  table(GSE150774_gene_info$V2 %in% GSE253355_gene_info$V2)# FALSE/TRUE:12200/20538，注释差异明显
  
  GSE150774_gene_info$ref_symbol1=all_shared_ensembl_id_info[rownames(GSE150774_gene_info),'ref_symbol1']
  GSE253355_gene_info$ref_symbol1=all_shared_ensembl_id_info[rownames(GSE253355_gene_info),'ref_symbol1']
  
 
  table(rowSums(GSE150774_BM_seu@assays$RNA@layers$counts)>0)# FALSE:49
  GSE150774_BM_seu=CreateSeuratObject(GetAssayData(GSE150774_BM_seu,assay = 'RNA',layer = 'counts'),min.cells = 10,min.features = 200,meta.data = GSE150774_BM_seu@meta.data)
  
  table(duplicated(GSE150774_gene_info[,c('V2','ref_symbol1')]))# TRUE: 26 
  table(duplicated(GSE253355_gene_info[,c('V2','ref_symbol1')]))# TRUE: 13 
  GSE150774_gene_info=unique(GSE150774_gene_info[,c('V2','ref_symbol1')])
  GSE253355_gene_info=unique(GSE253355_gene_info[,c('V2','ref_symbol1')])
  
  new_symbol_list=list()
  
  
  new_symbol_list[['GSE150774']]=GSE150774_gene_info$ref_symbol1
  names(new_symbol_list[['GSE150774']])=GSE150774_gene_info$V2
  new_symbol_list[['GSE253355']]=GSE253355_gene_info$ref_symbol1
  names(new_symbol_list[['GSE253355']])=GSE253355_gene_info$V2
  new_symbol_list[['collected_ABM']]=all_shared_ensembl_id_info$ref_symbol1
  names(new_symbol_list[['collected_ABM']])=all_shared_ensembl_id_info$ref_symbol1
  
  object_list=list('GSE150774'=GSE150774_BM_seu,'GSE253355'=Cell_GSE253355_ABM_seu,'collected_ABM'=sub_BM_seu_new)
  
  
  
  # 由于ananotation symbol有重复，所以symbol后面带有后缀，先去除重复
  # 先处理原本矩阵中的symbol重复
  library(tidyverse)
  for( refID in names(object_list)){
      table(rownames(object_list[[refID]]) %in% names(new_symbol_list[[refID]]))
      duplicated_symbol=rownames(object_list[[refID]])[ ! rownames(object_list[[refID]]) %in% names(new_symbol_list[[refID]])]
      togene_symbol=as.character(t(data.frame(strsplit(duplicated_symbol,split = '.',fixed = T)))[,1])
      
      tmp_assay=GetAssayData(object_list[[refID]],assay = 'RNA',layer = 'counts')
      du_tmp_assay=tmp_assay[rownames(tmp_assay) %in% c(togene_symbol,duplicated_symbol), ]
      tmp_assay=tmp_assay[!rownames(tmp_assay) %in% c(togene_symbol,duplicated_symbol),]
      
      du_tmp_assay=aggregate(du_tmp_assay,by=list(as.character(t(data.frame(strsplit(rownames(du_tmp_assay),split = '.',fixed = T)))[,1])),FUN=sum)
      du_tmp_assay=column_to_rownames(du_tmp_assay,'Group.1')
      du_tmp_assay=as(as.matrix(du_tmp_assay),'dgCMatrix')
      
      tmp_assay=rbind(tmp_assay,du_tmp_assay)
      object_list[[refID]]= CreateSeuratObject(counts = tmp_assay,meta.data =object_list[[refID]]@meta.data )
  }
  
  # 再对处理后的矩阵 gene symbol重命名，而重命名的gene symbol 也会存在重复，同样需要再处理
  
  for(refId in names(object_list)){
    tmp_assay=GetAssayData(object_list[[refId]],assay = 'RNA',layer = 'counts')
    rownames(tmp_assay)=as.character(new_symbol_list[[refId]][ rownames(tmp_assay)])
    
    du_symbol= rownames(tmp_assay)[duplicated( rownames(tmp_assay))];print(length(du_symbol))
    du_tmp_assay=tmp_assay[rownames(tmp_assay) %in% du_symbol,]
    tmp_assay=tmp_assay[!rownames(tmp_assay) %in% du_symbol,]
    
    du_tmp_assay=aggregate(du_tmp_assay,by=list(rownames(du_tmp_assay)), FUN=sum) 
    du_tmp_assay=column_to_rownames(du_tmp_assay,'Group.1')
    du_tmp_assay=as(as.matrix(du_tmp_assay),'dgCMatrix')
    
    object_list[[refId]]=CreateSeuratObject(counts = tmp_assay,min.cells = 10,min.features = 200,meta.data =object_list[[refId]]@meta.data )
  
}
  
  
  
  object_list[[1]]$resource='GSE150774';object_list[[1]]$sample=object_list[[1]]$donor;
  object_list[[2]]$resource='GSE253355';object_list[[2]]$sample=object_list[[2]]$orig.ident;
  object_list[['collected_ABM']]$resource=object_list[['collected_ABM']]$projectId;
  
  saveRDS(object_list[['GSE150774']],file='NRBC_BM_altas/tmp_GSE150774_BM_seu.rds')
  saveRDS(object_list[['GSE253355']],file='NRBC_BM_altas/tmp_Cell_GSE253355_ABM_seu.rds')# GSE253355 红细胞不在 ABM_altas_seu，filt_NBRC_altas_seu ？
  
  table(rownames(object_list[['GSE150774']]) %in% all_shared_ensembl_id_info$ref_symbol1)# FALSE  TRUE: 0/12357 
  table(rownames(object_list[['GSE253355']]) %in% all_shared_ensembl_id_info$ref_symbol1)#FALSE  TRUE:10 27872 
  table(rownames(object_list[['collected_ABM']]) %in% all_shared_ensembl_id_info$ref_symbol1)#FALSE  TRUE:10 27872 
  
  sub_BM_seu_new=object_list[['collected_ABM']]
  saveRDS(sub_BM_seu_new,file='NRBC_BM_altas/tmp_sub_ABM_seu_new.rds')
  
     
  rm(GSE150774_BM_seu,Cell_GSE253355_ABM_seu,sub_BM_seu_new)
  rm(du_tmp_assay,tmp_assay,new_symbol_list);gc()
  
  ABM_altas_seu=merge(object_list[[1]],c(object_list[2],object_list[[3]]))
  rm(object_list);gc()
  ABM_altas_seu$subcelltype[is.na(ABM_altas_seu$subcelltype)]=ABM_altas_seu$celltype1[is.na(ABM_altas_seu$subcelltype)]
  

   
  ABM_altas_seu[['RNA']]=JoinLayers(ABM_altas_seu[['RNA']])
  ABM_altas_seu=NormalizeData(ABM_altas_seu) %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData() %>% RunPCA() %>% RunUMAP(dims=1:30,reduction.name ='raw_umap' )
  DimPlot(ABM_altas_seu,reduction = 'raw_umap',group.by = 'resource',raster = F)

  anchor=FindTransferAnchors(reference =FBM_altas_seu ,query =ABM_altas_seu,reference.reduction  = 'pca',dims = 1:20 )
  ABM_altas_seu=MapQuery(anchorset =anchor,query = ABM_altas_seu,reference = FBM_altas_seu,reference.reduction = 'pca',reduction.model = 'umap1',refdata = list(celltype='subcelltype')) # 注意，command 保存的是umap1 信息，如果是umap2模块信息，会混淆
  rm(anchor);gc()
  p=DimPlot(ABM_altas_seu,reduction = 'ref.umap',group.by = c('predicted.celltype','subcelltype'),raster=FALSE,cols = cols)
  p
  ggsave(p,file='NRBC_BM_altas/res_pic/ABM_altas_predicted_subcelltype_refumap.pdf',width = 16,height = 6)
  
  saveRDS(ABM_altas_seu,file='NRBC_BM_altas/tmp_ABM_altas_seu_new.rds')

}else{
  ABM_altas_seu=readRDS('NRBC_BM_altas/tmp_ABM_altas_seu_new.rds')
}
table(rownames(ABM_altas_seu) %in% all_shared_ensembl_id_info$ref_symbol1)# FALSE/TRUE : 10/27090
table(rownames(ABM_altas_seu) %in% all_shared_ensembl_id_info$ref_symbol)# FALSE/TRUE : 7419/19681
subset(ABM_altas_seu,subset=GYPA >1 & IGKC >1) # 732 细胞满足这个条件，导致被删除
table(subset(ABM_altas_seu,subset=GYPA >1 & IGKC >1)$subcelltype)

VlnPlot(subset(ABM_altas_seu,predicted.celltype %in% c('BFUE/CFUE','ProE','Bas','Poly','Orth' ) ),features =c('TFRC','GYPA','IGHA1','IGHG1','IGKC') ,group.by = 'predicted.celltype',stack = T)
ABM_altas_seu[['prediction.score.celltype']]=NULL


##############################################################################################################################################
#--------------------------------------BM altas altas----------------------------------#
##############################################################################################################################################


if(F){
  ABM_altas_seu=subset(ABM_altas_seu, subcelltype  %in% unique(ABM_altas_seu$subcelltype)[!unique(ABM_altas_seu$subcelltype) %in% c("RBC","Erythroblast")])
  BM_altas_seu=merge(fetal_bm_seu,c(FBM_altas_seu,ABM_altas_seu))
  table(rownames(BM_altas_seu) %in%  all_shared_ensembl_id_info$ref_symbol1)# T/F:28130/10
  rownames(BM_altas_seu)[!rownames(BM_altas_seu) %in% all_shared_ensembl_id_info$ref_symbol1] # --代替下划线__
  BM_altas_seu[['RNA']]=JoinLayers(BM_altas_seu[['RNA']])
  BM_altas_seu[['prediction.score.celltype']]=NULL
  
  
  # FBM PCA有50 dimension,umap1 model 只是使用了20 个PCs
  FBM_pca_dim20=CreateDimReducObject(embeddings =FBM_altas_seu[['pca']]@cell.embeddings[,1:20],loadings = FBM_altas_seu[['pca']]@feature.loadings[,1:20],stdev = FBM_altas_seu[['pca']]@stdev)
  BM_altas_seu[['pca']]=merge(FBM_pca_dim20,c(fetal_bm_seu[['ref.pca']],ABM_altas_seu[['ref.pca']]))
  BM_altas_seu[['umap']]=merge(FBM_altas_seu[['umap1']],c(fetal_bm_seu[['ref.umap']],ABM_altas_seu[['ref.umap']]))
  
  
  
  #BM_altas_seu$subcelltype=factor(BM_altas_seu$subcelltype,levels = c(subcelltype_levels[1:6],"YS_BFUE/CFUE_like",subcelltype_levels[7:length(subcelltype_levels)]))
  BM_altas_seu$stage='ABM'
  BM_altas_seu$stage[BM_altas_seu$resource=='GSE143753' ]='EBM'
  BM_altas_seu$stage[BM_altas_seu$resource=="E-MTAB-11343"]='FBM' 
  BM_altas_seu$stage=factor(BM_altas_seu$stage,levels = c('EBM','FBM','ABM'))
  BM_altas_seu$predicted.celltype[BM_altas_seu$stage=='FBM']=BM_altas_seu$subcelltype[BM_altas_seu$stage=='FBM']
  table(BM_altas_seu$stage)
  
  BM_altas_seu$age=as.character(BM_altas_seu$age)
  BM_altas_seu$age[BM_altas_seu$age=='CS22b']='CS22'
  BM_altas_seu$age[BM_altas_seu$age=='CS15b']='CS15'
  BM_altas_seu$age[is.na(BM_altas_seu$age)]=BM_altas_seu$Age[is.na(BM_altas_seu$age)]
  
  BM_altas_seu$age[is.na(BM_altas_seu$age)]='unknown'
  BM_altas_seu$age[! BM_altas_seu$age %in% c( 'CS13','CS15','CS20','CS22','unknown')]=paste0(BM_altas_seu$age[! BM_altas_seu$age %in% c( 'CS13','CS15','CS20','CS22','unknown')],'y')
  ABM_altas_seu$age[BM_altas_seu$stage=='FBM']=gsub(pattern ='y' ,replacement = 'PCW',BM_altas_seu$age[BM_altas_seu$stage=='FBM'])
  BM_altas_seu$age[BM_altas_seu$resource %in% c('GSE150774','GSE133181','GSE169426','GSE181989')]='adult'
  table(BM_altas_seu$age)
  
  
  
  rm(FBM_altas_seu,fetal_bm_seu,ABM_altas_seu);gc()
  
  
  age_levels= c("CS13","CS15","CS20" ,"CS22", "12PCW","14PCW","15PCW","16PCW","17PCW", '2y','25y','31y','41y','44y',"52y","54y","57y","59y","64y","65y","66y","71y","72y","73y","74y",'adult')
  sort(unique(BM_altas_seu$age)[!unique(BM_altas_seu$age) %in% age_levels])
  BM_altas_seu$age=factor(BM_altas_seu$age,levels =age_levels)
  
  table(BM_altas_seu@meta.data[is.na(BM_altas_seu@meta.data$sample),c('resource')])
  BM_altas_seu$resource[ BM_altas_seu$stage=='FBM']='E-MTAB-11343'
  BM_altas_seu$sample[is.na(BM_altas_seu$sample)]=BM_altas_seu$donor[is.na(BM_altas_seu$sample)]
  
  # HBA1 极致高表达，以及网致红细胞数目多，GYPA表达降低，红细胞破裂后HBA1污染更加严重
  douletes_cells=WhichCells(BM_altas_seu,expression = GYPA >1 & IGKC >1 );length(douletes_cells) # 剔除杂细胞，针对红细胞,可能是谱系不忠或者谱系泄露,ABM_altas_NRBC_seu忠IGKC >1 :2028 sample, 主要来自bas和BFUE/CFUE，可能是一类特殊亚类
  douletes_BM_altas_seu=subset(BM_altas_seu,cells=douletes_cells)
  table(douletes_BM_altas_seu@meta.data[,c('resource','subcelltype')])# 主要来在GSE253355：ProE:78/Bas:494/Poly:15/Orth:59 
  table(douletes_BM_altas_seu@meta.data[,c('sample','subcelltype')]) # 主要来自H14_MACS样本：Bas： 240
  IGKC_pos_Ery_cellsid=rownames(douletes_BM_altas_seu@meta.data[douletes_BM_altas_seu$subcelltype %in% c('BFUE/CFUE','ProE','Bas','Poly','Orth'),])
  douletes_BM_altas_seu@meta.data[IGKC_pos_Ery_cellsid,'subcelltype']=paste0('IGKC+',douletes_BM_altas_seu@meta.data[IGKC_pos_Ery_cellsid,'subcelltype'])
  
  VlnPlot(douletes_BM_altas_seu,features = c('IGKC','HBA1','IGHA1'),stack = T)
  BM_altas_seu=subset(BM_altas_seu,cells=colnames(BM_altas_seu)[!colnames(BM_altas_seu) %in% douletes_cells])
  
   
  table(BM_altas_seu@meta.data[is.na(BM_altas_seu$subcelltype),'stage'])
  BM_altas_seu$subcelltype[BM_altas_seu$stage=='EBM']=BM_altas_seu$celltype[BM_altas_seu$stage=='EBM']
  
  p=DimPlot(BM_altas_seu,group.by = c( 'stage','resource','sample'),raster =F,cols = col);p
  ggsave(p,filename = 'NRBC_BM_altas/res_pic/BM_altas_stage_resource_sample_umap.pdf',width =16,height = 8,dpi = 300)
  
  
  DimPlot(BM_altas_seu,group.by = c( 'subcelltype'),raster = T)
  
  BM_altas_seu$predicted.celltype=as.character(BM_altas_seu$predicted.celltype)
  table(is.na(BM_altas_seu$predicted.celltype))
  BM_altas_seu$predicted.celltype[is.na(BM_altas_seu$predicted.celltype)]=BM_altas_seu$subcelltype[is.na(BM_altas_seu$predicted.celltype)]
  
  table(BM_altas_seu@meta.data[,c('subcelltype','resource')])
  table(BM_altas_seu@meta.data[,c('subcelltype','stage')])
  table(BM_altas_seu@meta.data[is.na(BM_altas_seu$subcelltype),'stage'])
  
  head(BM_altas_seu@meta.data[is.na(BM_altas_seu$subcelltype),])

  #BM_altas_seu$subcelltype[is.na(BM_altas_seu$subcelltype)]=BM_altas_seu$celltype1[is.na(BM_altas_seu$subcelltype)]
  #sort(unique(BM_altas_seu$subcelltype))
  
  BM_altas_seu$new_celltype=BM_altas_seu$subcelltype
  sort(unique( BM_altas_seu$new_celltype))
  BM_altas_seu$new_celltype[BM_altas_seu$subcelltype %in% c('Macrophages','MACROPHAGE')]='MACROPHAGE'
  BM_altas_seu$new_celltype[BM_altas_seu$subcelltype %in% c('AEC','SEC','ENDOTHELIUM')]='ENDOTHELIUM'
  BM_altas_seu$new_celltype[BM_altas_seu$subcelltype %in% c("B cell","B CELLs","B CELLS","Pre-B","Pre-Pro B" ,"Pro-B" ,"Mature B" )] ="B CELLs"
  BM_altas_seu$new_celltype[BM_altas_seu$subcelltype %in% c("MONOCYTES","monocyte",'MONOCYTE','Monocyte') ] ='MONOCYTE'
  BM_altas_seu$new_celltype[BM_altas_seu$subcelltype %in% c("Proliferation myeloid cell","Granulocyte","Early Myeloid Progenitor","Late Myeloid") ] ="MYELOCYTE" 
  
  BM_altas_seu$new_celltype[BM_altas_seu$subcelltype %in% c("Platelet",'MEGAKARYOCYTE','Meg',"Megakaryocyte") ] ="MEGAKARYOCYTE" 
  BM_altas_seu$new_celltype[BM_altas_seu$subcelltype %in% c("Osteoblast","OSTEOBLAST") ] ="OSTEOBLAST" 
  BM_altas_seu$new_celltype[BM_altas_seu$subcelltype %in% c("MACROPHAGE_ERY","Mac_Ery") ] ="Mac_Ery" 
  BM_altas_seu$new_celltype[BM_altas_seu$subcelltype %in% c("Neutrophil","NEUTROPHIL") ] ="NEUTROPHIL" 
  BM_altas_seu$new_celltype[BM_altas_seu$subcelltype %in% c("Eosinophil","Ba/Eo/Ma" ) ] ="EO/BASO/MAST"
  BM_altas_seu$new_celltype[BM_altas_seu$subcelltype %in% c("CD4+ T-Cell" ) ] ="CD4 T"
  BM_altas_seu$new_celltype[BM_altas_seu$subcelltype %in% c("CD8+ T-Cell" ) ] ="CD8 T"
  BM_altas_seu$new_celltype[BM_altas_seu$subcelltype %in% c("pDC","Cycling DCs",'DC' ) ] ="DCs"
  BM_altas_seu$new_celltype[BM_altas_seu$new_celltype %in% c("CYCLING_MPP","HSC_MPP",'HSC','MPP',"Cycling HSPC") ] ="HSC_MPP" 
  BM_altas_seu$new_celltype[BM_altas_seu$new_celltype %in% c("CYCLING_MEMP",'MEMP') ] ="MEMP"
  BM_altas_seu$new_celltype[BM_altas_seu$new_celltype %in% c("Plasma Cell"  ) ] ="Plasma cell"
  BM_altas_seu$new_celltype[BM_altas_seu$new_celltype=="YS_Bas"]="YS_Bas/Poly"
  
  celltype_level=c("HSC_MPP","CMP","MEMP" ,"GMP" ,"LMPP_MLP" ,'MEP',"BFUE/CFUE","ProE","Bas","Poly","Orth","YS_Bas/Poly","YS_Orth",
                   "MEGAKARYOCYTE","EO/BASO/MAST","MYELOCYTE","NEUTROPHIL",
                   "MOP","MONOCYTE","DCs","MACROPHAGE","Mac_Ery","ILC","CLP","B CELLs","Plasma cell","CD4 T","CD8 T","Treg","NK/T CELLS",'NK',"Myoprogenitor" ,"MYOCYTE",
                   "Osteoprogenitors","CHONDROBLAST","CHONDROCYTE", "OSTEOBLAST","OSTEOCLAST" ,"VSMC",'PMSC1','PMSC2','OCPs',"BMSC1","BMSC2" ,"Schwan cell",
                   "ENDOTHELIUM", "Adipo-MSC","APOD+ MSC","Fibro-MSC","Osteo-MSC","THY1+ MSC","RNAlo MSC","FIBROBLASTS","SKELETAL MUSCLE","EPITHELIUM" )
  table(unique(BM_altas_seu$new_celltype) %in%  celltype_level)
  unique(BM_altas_seu$new_celltype)[! unique(BM_altas_seu$new_celltype) %in%  celltype_level]
  BM_altas_seu=subset(BM_altas_seu,new_celltype %in% celltype_level)
  
  BM_altas_seu$new_celltype=factor(BM_altas_seu$new_celltype,levels =celltype_level )
  
  p=DimPlot(BM_altas_seu,group.by = c('new_celltype','stage'),cols = cols,raster = T);p
  ggsave(p,file='NRBC_BM_altas/res_pic/BM_altas_celltype_stage_umap.pdf',width = 16,height = 6)
  
  p=DimPlot(BM_altas_seu,group.by = c('new_celltype','predicted.celltype'),cols = cols,raster = T);p
  ggsave(p,file='NRBC_BM_altas/res_pic/BM_altas_celltype_predicted.celltypee_umap.pdf',width = 16,height = 6)
  
  saveRDS(BM_altas_seu,file = 'NRBC_BM_altas/BM_altas_seu_v2.rds')
  saveRDS(BM_altas_seu@meta.data ,file = 'NRBC_BM_altas/BM_altas_seu.meta.data_v2.rds')

  
  BM_altas_seu_meta=BM_altas_seu@meta.data
  BM_altas_seu_meta$id=paste(BM_altas_seu_meta$age,BM_altas_seu_meta$sample,sep = "_")
  BM_celltype_Fre_df=data.frame(table(BM_altas_seu_meta[,c('id','new_celltype')]))
  BM_celltype_Fre_df= BM_celltype_Fre_df[BM_celltype_Fre_df$Freq >0,]
  BM_celltype_Fre_df$id=factor(BM_celltype_Fre_df$id,levels = id_levels)
  
  BM_celltype_age_df=data.frame(table(BM_altas_seu_meta$id));colnames(BM_celltype_age_df)=c('id','count')
  
  p=ggplot(BM_celltype_Fre_df,aes(x=id,fill=new_celltype,y=Freq))+geom_bar(stat ='identity',position = 'fill' )+theme_classic()+scale_fill_manual(values =  col)+
    theme(axis.text.x =  element_text(angle = 30,vjust = 0.85,hjust = 0.75),axis.text = element_text(face = 'bold'))+
    geom_label(data = BM_celltype_age_df,aes(x=id,y=1.06,label=count,fill=NULL),vjust ='top',show.legend = F,label.size = NA)+scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1))
  p
  ggsave(p,filename='res_pic/BM_NRBC_Freq_age_sample_barplot.pdf',width=18,height = 6)
  
}else{
  BM_altas_seu=readRDS('NRBC_BM_altas/BM_altas_seu_v2.rds')#
}



####################################################################################################################################
#-------------------------------------------------- NRBC of BM  analysis -------------------#
####################################################################################################################################
table(BM_altas_seu_meta$new_celltype)
Ery_subcelltype=c("BFUE/CFUE","ProE","Bas","Poly","Orth","YS_Bas/Poly", "YS_Orth")

BM_altas_seu=subset( BM_altas_seu,new_celltype %in% Ery_subcelltype)
DimPlot(BM_NRBC_altas_seu,group.by = c('new_celltype'),cols = col,raster=FALSE,split.by = 'stage')
BM_NRBC_altas_seu$subcelltype=BM_NRBC_altas_seu$new_celltype
BM_NRBC_altas_seu$subcelltype=factor(BM_NRBC_altas_seu$subcelltype,exclude = levels(BM_NRBC_altas_seu$subcelltype)[ !levels(BM_NRBC_altas_seu$subcelltype) %in% unique(BM_NRBC_altas_seu$subcelltype)])
  
BM_NRBC_altas_seu=singleR_analysis_func(refdata =nrbc_ref_se,test =as.matrix(BM_NRBC_altas_seu@assays$RNA$data)[VariableFeatures(BM_NRBC_altas_seu),],outdata = BM_NRBC_altas_seu,an_type1 ='celltype', an_type2 = 'Pre_celltype')
BM_NRBC_altas_seu=BM_NRBC_altas_seu[[1]]

rm(BM_altas_seu);gc();rm(FBM_altas_seu);gc()

# 重新处理效果一般
BM_NRBC_altas_seu=RunPCA(BM_NRBC_altas_seu)
BM_NRBC_altas_seu$donor[is.na(BM_NRBC_altas_seu$donor)]=BM_NRBC_altas_seu$sample[is.na(BM_NRBC_altas_seu$donor)]
BM_NRBC_altas_seu=RunHarmony(BM_NRBC_altas_seu,group.by.vars=c('resource','sample'))
BM_NRBC_altas_seu=RunUMAP(BM_NRBC_altas_seu,reduction = 'harmony',dims = 1:6,reduction.name = 'rumap')

p=DimPlot(BM_NRBC_altas_seu,group.by = c('new_celltype'),cols = col,reduction = 'rumap',pt.size=0.8);p
#ggsave(p,filename = 'NRBC_BM_altas/res_pic/BM_NRBC_altas_celltype_umap.pdf',width = 8,height = 8,dpi = 300)


#BM中，不能有效区分出YS_NRBC 与def——NRBC，所以BM NRBC 应该mapping 到 FL NRBC中查看看
sub_F61_Ery_seu =readRDS('NRBC_FL_altas/tmp_sub_F61_Ery_seu_umap_model.rds')
DimPlot(sub_F61_Ery_seu,reduction = 'umap1',group.by = 'subcelltype',cols = cols)+ggtitle('FL_CS18_F61')+FeaturePlot(sub_F61_Ery_seu,features = c('HBE1'),reduction = 'umap1')

#BM_NRBC_altas_seu=NormalizeData(BM_NRBC_altas_seu) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

anchors=FindTransferAnchors(reference =sub_F61_Ery_seu,query = BM_NRBC_altas_seu,reference.reduction = 'pca' )
BM_NRBC_altas_seu=MapQuery(anchorset =anchors,query =BM_NRBC_altas_seu,reference = sub_F61_Ery_seu, reference.reduction = 'pca',reduction.model = 'umap1',refdata = list(celltype='subcelltype') )
rm(anchors)


p=DimPlot(BM_NRBC_altas_seu,reduction = 'ref.umap',group.by = c('subcelltype'),cols = cols)+VlnPlot(BM_NRBC_altas_seu,features = c('KIT','TFRC','GYPA','HBE1','CD63'),group.by = 'subcelltype',stack = T)+NoLegend();p
ggsave(p,filename='NRBC_BM_altas/res_pic/BM_umap_celltype_vlnplot_feature.pdf',width = 12,height = 6,dpi = 300)

p=DimPlot(BM_NRBC_altas_seu,reduction = 'ref.umap',group.by = c('stage','sample'),cols = cols)+FeaturePlot(BM_NRBC_altas_seu,features = c('HBE1'),reduction = 'ref.umap')
p
ggsave(p,filename='NRBC_BM_altas/res_pic/BM_mapping_FL_F61_refumap.pdf',width = 24,height = 8,dpi = 300)

p=DimPlot(BM_NRBC_altas_seu,reduction = 'ref.umap',group.by = c('subcelltype'),cols = cols,ncol = 1,split.by = 'stage')
p
ggsave(p,filename='NRBC_BM_altas/res_pic/stage_BM_mapping_FL_F61_refumap.pdf',width = 8,height = 24,dpi = 300)

p1=DimPlot(subset(BM_NRBC_altas_seu,stage=='FBM'),reduction = 'ref.umap',group.by = c('subcelltype'),cols = cols)+ggtitle('FBM:18,365')
p2=DimPlot(subset(BM_NRBC_altas_seu,stage=='ABM'),reduction = 'ref.umap',group.by = c('subcelltype'),cols = cols)+ggtitle('ABM:16,089')
p=p1/p2;p
ggsave(p,filename='NRBC_BM_altas/res_pic/BM_celltype_umap.pdf',width = 6,height = 10,dpi = 300)


unique(BM_NRBC_altas_seu$sample)# 

#unique(BM_NRBC_altas_seu@meta.data[,c('sample','resource')])
#id=c( "GSM3901485", "GSM3901486" ,"GSM3901487", "GSM3901488","GSM3901492", "GSM5202203", "GSM5202204", "GSM3993352","GSM3993353" ,"GSM3993354","GSM3993355",'GSM5047372','GSM5515743','GSM5515744' )
#names(id)=c('CD34+BM1','CD34+BM2','CD34+BM3','CD34+BM4','CD34+C1','CD34+Ctrl1','CD34+Ctrl3','CD34+H1','CD34+H2','CD34+H3','CD34+H4','BM_HD_2253','BMH1','BMH2')
  
BM_NRBC_altas_seu$id=paste(BM_NRBC_altas_seu$age,BM_NRBC_altas_seu$sample,sep='_')
id_levels=c(sort(unique(BM_NRBC_altas_seu$id))[38:43],sort(unique(BM_NRBC_altas_seu$id))[c(1:8,10,9,11:25)],sort(unique(BM_NRBC_altas_seu$id))[29:37],sort(unique(BM_NRBC_altas_seu$id))[26:28])
BM_NRBC_altas_seu$id[BM_NRBC_altas_seu$id %in% id_levels[15:40]]=paste0(BM_NRBC_altas_seu$id[BM_NRBC_altas_seu$id %in% id_levels[15:40]],'_CD34+')
BM_NRBC_altas_seu$id[BM_NRBC_altas_seu$id %in% id_levels[41:43]]=paste0(BM_NRBC_altas_seu$id[BM_NRBC_altas_seu$id %in% id_levels[41:43]],'_GYPA+')
id_levels[id_levels %in% id_levels[15:40]]=paste0(id_levels[id_levels%in% id_levels[15:40]],'_CD34+')
id_levels[id_levels %in% id_levels[41:43]]=paste0(id_levels[id_levels%in% id_levels[41:43]],'_GYPA+')
BM_NRBC_altas_seu$id=factor(BM_NRBC_altas_seu$id,id_levels)

BM_NRBC_altas_seu=subset(BM_NRBC_altas_seu,subset=IGKC <=1 ) #因为要比较的是健康青壮年，与老人的差异显著，暂时不做此分析研究比较，专门做一个论文课题分析


Ery_age_celltype_df=data.frame(table(BM_NRBC_altas_seu@meta.data[,c('id','new_celltype')]))
Ery_age_celltype_df=Ery_age_celltype_df[Ery_age_celltype_df$Freq>0,]
Ery_age_celltype_df$id=factor(Ery_age_celltype_df$id,levels =id_levels )
celltype_age_df=data.frame(table(BM_NRBC_altas_seu$id));colnames(celltype_age_df)=c('id','count')
celltype_age_df=celltype_age_df[celltype_age_df$count>0,]

p=ggplot(Ery_age_celltype_df,aes(x=id,fill=new_celltype,y=Freq))+geom_bar(stat ='identity',position = 'fill' )+theme_classic()+scale_fill_manual(values =  col)+
  theme(axis.text.x =  element_text(angle = 30,vjust = 0.85,hjust = 0.75),axis.text = element_text(face = 'bold'))+
  geom_label(data = celltype_age_df,aes(x=id,y=1.06,label=count,fill=NULL),vjust ='top',show.legend = F,label.size = NA)+scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1))
p

ggsave(p,filename='NRBC_BM_altas/res_pic/BM_NRBC_Freq_age_sample_barplot.pdf',width=24,height = 6)

temp_df=data.frame(table(BM_NRBC_altas_seu@meta.data[BM_NRBC_altas_seu$stage=='FBM',c('id','Sort_id')]))
temp_df=temp_df[temp_df$Freq >0,]
p=ggplot(temp_df,aes(x=id,y=Freq,color=Sort_id,group=Sort_id))+geom_point()+geom_line(linewidth=0.8)+theme_classic()+
  RotatedAxis()+scale_color_manual(values = cols[-2])+theme(text = element_text(face = 'bold')) +scale_linetype_manual(values = c(1,2))#  log2(Freq) 更能反应波动
p
ggsave(p,filename='NRBC_BM_altas/res_pic/FBM_NRBC_Freq_dotline.pdf',width=6,height = 6,dpi = 300)

BM_NRBC_altas_seu=CellCycleScoring(BM_NRBC_altas_seu,s.features = cc.genes$s.genes,g2m.features = cc.genes$g2m.genes,set.ident = F)
DimPlot(BM_NRBC_altas_seu,group.by = 'Phase',reduction = 'ref.umap',cols = cols,split.by = 'stage')

temp_df=data.frame(table(BM_NRBC_altas_seu@meta.data[BM_NRBC_altas_seu$Sort_id=='CD45P',c('id','new_celltype','Sort_id')]))
temp_df=temp_df[temp_df$Freq >0,]
celltype_count_df=data.frame(table(BM_NRBC_altas_seu@meta.data[BM_NRBC_altas_seu$Sort_id=='CD45P',c('id')]))
celltype_count_df=celltype_count_df[celltype_count_df$Freq >0,];colnames(celltype_count_df)=c('id','count')

p=ggplot(temp_df,aes(x=id,fill=new_celltype,y=Freq))+geom_bar(stat ='identity',position = 'fill' )+theme(text = element_text(face = 'bold'))+
  theme_classic()+scale_fill_manual(values =  col)+RotatedAxis()+ggtitle('FBM nRBC: CD45P sorting')+scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1))+
  geom_label(data = celltype_count_df,aes(x=id,y=1.16,label=count,fill=NULL),vjust ='top',show.legend = F,label.size = NA)
p
ggsave(p,filename='NRBC_BM_altas/res_pic/FBM_NRBC_CD45P_barplot.pdf',width=8,height = 6,dpi = 300)


saveRDS(BM_NRBC_altas_seu,file='NRBC_BM_altas/res_data/BM_NRBC_altas_seu.rds')

rm(list=ls());gc()





############################################################################################################################################################################
#-----------------------------------------cellchat analysis--------------------------------#
############################################################################################################################################################################

BM_altas_seu=readRDS('NRBC_BM_altas/BM_altas_seu_v2.rds')# GSE253355 少了大半细胞,IGKC+NRBC 占比很高，都是老年人股骨来源细胞
ABM_altas_seu=subset(BM_altas_seu,stage=='ABM');rm(BM_altas_seu);gc()
ABM_altas_seu$orig.ident[ABM_altas_seu$orig.ident=='SeuratProject']=ABM_altas_seu$donor[ABM_altas_seu$orig.ident=='SeuratProject']
del_cellids=subset(ABM_altas_seu, new_celltype %in% c('BFUE/CFUE','ProE','Bas','Poly','Orth') )
del_cellids=colnames(del_cellids)[!colnames(del_cellids) %in% rownames(BM_NRBC_altas_seu@meta.data)]
ABM_altas_seu=subset(ABM_altas_seu,cells= colnames(ABM_altas_seu)[!colnames(ABM_altas_seu) %in% del_cellids ])
saveRDS(ABM_altas_seu,file = 'NRBC_BM_altas/filt_ABM_altas_seu_20251127.rds')
rm(BM_NRBC_altas_seu);gc()


future.seed=TRUE
library(CellChat)
library(Seurat)
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB # simply use the default CellChatDB,
options(future.globals.maxSize=80000 * 1024 ^ 2 )

table(ABM_altas_seu$resource[is.na(ABM_altas_seu$ct)])
cho_cells=sample(rownames(ABM_altas_seu@meta.data),size=50000)# 12 个样本，抽取3/12=1/4
ABM_altas_seu=subset(ABM_altas_seu,cells=cho_cells)
ABM_altas_seu$new_celltype=droplevels(ABM_altas_seu$new_celltype, exclude = setdiff(levels(ABM_altas_seu$new_celltype),unique(ABM_altas_seu$new_celltype)))
cellchat=createCellChat(ABM_altas_seu,assay = 'RNA',group.by='new_celltype')
rm(ABM_altas_seu);gc()

CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat@DB=CellChatDB.use
cellchat=subsetData(cellchat)
future::plan("multisession", workers = 16) # do parallel
cellchat=identifyOverExpressedGenes(cellchat)
cellchat=identifyOverExpressedInteractions(cellchat)
cellchat=projectData(cellchat,PPI.human)
cellchat= computeCommunProb(cellchat)
cellchat=filterCommunication(cellchat, min.cells = 20)
cellchat=computeCommunProbPathway(cellchat)
cellchat=aggregateNet(cellchat)
cellchat=netAnalysis_computeCentrality(cellchat)

saveRDS(cellchat,file='NRBC_BM_altas/ABM_subcelltype_cellchat.rds')

