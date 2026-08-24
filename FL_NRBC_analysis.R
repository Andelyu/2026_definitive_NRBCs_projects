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
library(limma)
library(harmony)

cols=c(brewer.pal(12,"Set3"),brewer.pal(9,"Set1"),brewer.pal(6,"PiYG"),brewer.pal(6,"BrBG"),brewer.pal(8,"Set2"),
       brewer.pal(8,"Pastel2"),brewer.pal(9,"Pastel1"),brewer.pal(8,"Accent"))
col=unique(cols)[-14]

dir.create('NRBC_FL_altas/res/res_pic',recursive = T)
dir.create('NRBC_FL_altas/res/res_data',recursive = T)

all_shared_ensembl_id_info=read.csv('ref_data/all_shared_ensembl_id_info.csv',header=T);all_shared_ensembl_id_info=all_shared_ensembl_id_info[,-1]
rownames(all_shared_ensembl_id_info)=all_shared_ensembl_id_info$X
all_shared_ensembl_id_info=all_shared_ensembl_id_info[,-1]
table(duplicated(all_shared_ensembl_id_info$ref_symbol1)) #497 symbol replicated 

####################################################################################################################################
#----------------------------------6-18WPC,fetus,E-MTAB-7407, 与Teichmann Lab 数据存在交集---------------------------------------------------#
####################################################################################################################################

if(F){
  #E-MTAB-7407
  # FL_altas_seu1: 13 human fetal livers (6-18 PCW)
  FL_altas_seu1 <- LoadH5Seurat(file = "/home/gibh/2021_NRBC_chlyu/ref_data/ref_scRNAseq_data/E-MTAB-7407_hema/fetal_liver_alladata_Copy.h5seurat")
  
  FL_altas_seu1=FindVariableFeatures(FL_altas_seu1,nfeatures = 2500) %>%ScaleData() %>%RunPCA() %>% RunHarmony(group.by.vars=c('orig.ident','sort.ids','stages','batch')) %>% RunUMAP(reductioin='harmony',reduction.name='umap1',dims=1:30)
  DimPlot(FL_altas_seu1,group.by = 'cell.labels',reduction = 'umap1',cols = cols,raster=T)+DimPlot(FL_altas_seu1,group.by = 'cell.labels',reduction = 'umap',cols = cols,raster=T)
  
  FL_altas_seu1=RunUMAP(FL_altas_seu1,reductioin='harmony',reduction.name='umap1',dims=1:20,return.model = T)
  DimPlot(FL_altas_seu1,group.by = 'cell.labels',reduction = 'umap1',cols = cols,raster=T)+DimPlot(FL_altas_seu1,group.by = 'cell.labels',reduction = 'umap',cols = cols,raster=T)
  #unique(FL_altas_seu1$cell.labels)
  cell_levels=c('HSC_MPP','MEMP','Mast cell','Megakaryocyte','Early Erythroid','Mid Erythroid','Late Erythroid',
                'Neutrophil-myeloid progenitor','DC precursor','pDC precursor','DC1','DC2','Monocyte precursor','Monocyte','Mono-Mac','Kupffer Cell','VCAM1+ EI macrophage',
                'Pre pro B cell',"pro-B cell",'pre-B cell','B cell','ILC precursor','Early lymphoid_T lymphocyte','NK',"Endothelial cell","Fibroblast","Hepatocyte")
  FL_altas_seu1$cell.labels=factor(FL_altas_seu1$cell.labels,levels = cell_levels)
  
  DimPlot(FL_altas_seu1,group.by = 'cell.labels',reduction = 'umap1',cols = cols[-8],raster=T)+DimPlot(FL_altas_seu1,group.by = 'cell.labels',reduction = 'umap',cols = cols[-8],raster=T)
  
  
  FL_altas_seu1=RunUMAP(FL_altas_seu1,reductioin='harmony',reduction.name='umap2',dims=1:25)
  DimPlot(FL_altas_seu1,group.by = 'cell.labels',reduction = 'umap2',cols = cols[-8],raster=T,label = T)+DimPlot(FL_altas_seu1,group.by = 'cell.labels',reduction = 'umap',cols = cols[-8],raster=T)
  
  DimPlot(FL_altas_seu1,group.by = 'cell.labels',reduction = 'umap2',cols = cols[-8],raster=T,label = T,split.by = 'stages')
  FeaturePlot(FL_altas_seu1,features = c('HBE1','GYPA','HBB'),reduction = 'umap1',split.by = 'stages',cols = c('gray','firebrick3'))
  unique(FL_altas_seu1$orig.ident)[!unique(FL_altas_seu1$orig.ident) %in% ys_fl_fbm_MEM_Ery_seu$donor]
  unique(FL_altas_seu1$fetal.ids)
  #计算对应的天数
  day_list=c(114,80,119,97,112,101,68,84,55,55,68,57,64,57)
  names(day_list)=as.character(unique(FL_altas_seu1@meta.data[,c('fetal.ids')]))
  FL_altas_seu1@meta.data$day=day_list[as.character(FL_altas_seu1@meta.data[,c('fetal.ids')])]
  saveRDS(FL_altas_seu1,file = '/home/gibh/2021_NRBC_chlyu/ref_data/ref_scRNAseq_data/E-MTAB-7407_hema/fetal_liver_alladata_Copy.rds')
  
}else{
  FL_altas_seu1=readRDS('/home/gibh/2021_NRBC_chlyu/ref_data/ref_scRNAseq_data/E-MTAB-7407_hema/fetal_liver_alladata_Copy.rds')
  DimPlot(FL_altas_seu1,group.by = 'cell.labels',reduction = 'umap',cols = cols,raster=T)
  dim(FL_altas_seu1) #  27080 113063
  table(rownames(FL_altas_seu1) %in% all_shared_ensembl_id_info$ref_symbol1)#FALSE/TRUE : 28/27052
  table(rownames(FL_altas_seu1) %in% all_shared_ensembl_id_info$ref_symbol)#FALSE/TRUE : 7110/19970
  
}

############################################################################################################################################
# -------------------Teichmann Lab: Mapping the developing human immune system across organs----------------#
############################################################################################################################################
# Teichmann Lab:FL_altas_seu2 比FL_altas_seu1 多一些cells, immune paper 新增加了部分FL 样本,
if(F){
  ys_fl_fbm_seu=readRDS('ref_data/ref_scRNAseq_data/ys_fl_fbm_seu.rds') #newly generated sequencing libraries: E-MTAB-11343,  包含 Decoding human fetal liver haematopoiesis：E-MTAB-7407
  FL_altas_seu2=subset(ys_fl_fbm_seu,organ=='LI');rm(ys_fl_fbm_seu);gc()
  head(HIS_genes)
  DimPlot(FL_altas_seu2,group.by = 'subcelltype')
  new_smapleid=unique(FL_altas_seu1$orig.ident)[ !unique(FL_altas_seu1$orig.ident) %in% unique( FL_altas_seu2$donor)] # F17 F16,两个样本不在FL_altas_seu2中，可以添加进来
  table(FL_altas_seu1$orig.ident) # F16/F17:1135/1992
  
  head(colnames(FL_altas_seu2));head(colnames(FL_altas_seu1));
  FL_altas_seu1_cellnames=gsub(pattern = '_',replacement = '-',paste0('FCA',as.character(t(data.frame(strsplit(colnames(FL_altas_seu1),split = '_FCA')))[,2])))
  names(FL_altas_seu1_cellnames)=colnames(FL_altas_seu1)
  table(as.character(FL_altas_seu1_cellnames) %in% colnames(FL_altas_seu2))
  #FALSE:   TRUE 
  #11766:  101297
  table(colnames(FL_altas_seu2) %in% as.character(FL_altas_seu1_cellnames))
  #FALSE   TRUE 
  #111304 101287 
  
  
  new_cellids=FL_altas_seu1_cellnames[!FL_altas_seu1_cellnames %in% colnames(FL_altas_seu2)];length(new_cellids) #E-MTAB-7407 部分cells 不在immune paper FL中
  F17_F16_cellids=rownames(FL_altas_seu1@meta.data[FL_altas_seu1$orig.ident %in% c('F16','F17'),]);length(F17_F16_cellids) # 3127
  new_cellids=new_cellids[ names(new_cellids)[! names(new_cellids) %in% F17_F16_cellids] ]
  new_FL_altas_seu1=subset(FL_altas_seu1,cells=c(names(new_cellids),F17_F16_cellids))
  
  length(F17_F16_cellids[grep(pattern = '4834STDY7',F17_F16_cellids)]) # 3127
  F17_F16_cellids=gsub(pattern = '_',replacement = '-',paste0('4834STDY7',as.character(t(data.frame(strsplit(F17_F16_cellids,split = '_4834STDY7')))[,2])))
  names(F17_F16_cellids)=rownames(FL_altas_seu1@meta.data[FL_altas_seu1$orig.ident %in% c('F16','F17'),])
  
  new_FL_altas_seu1=RenameCells(new_FL_altas_seu1,new.names =c(as.character(new_cellids),as.character(F17_F16_cellids)),old.names=c(names(new_cellids),names(F17_F16_cellids)))
  DimPlot(new_FL_altas_seu1,group.by = 'predicted_doublets')
  table(new_FL_altas_seu1$donor)
  
  rm(FL_altas_seu1)
  
  FL_altas_seu2_meta=FL_altas_seu2@meta.data
  FL_altas_seu2_meta=FL_altas_seu2_meta[,c('nCount_RNA','nFeature_RNA','mito', 'doublet_scores', 'predicted_doublets','Sample.lanes','Sort_id','age', 'method', 'donor', 'sex','uniform_label', 'AnnatomicalPart','anno_lvl_1_final_clean','anno_lvl_2_final_clean','anno_final_celltype2',  'subcelltype' )]
  FL_altas_seu2@meta.data=FL_altas_seu2_meta
  FL_altas_seu2$predicted_doublets[FL_altas_seu2$predicted_doublets=='False']='Singlet'
  FL_altas_seu2=subset(FL_altas_seu2,predicted_doublets=='Singlet')
  FL_altas_seu2$resource='immune_paper'
  rm(FL_altas_seu2_meta)
  
  new_FL_altas_seu1_meta=new_FL_altas_seu1@meta.data
  colnames(new_FL_altas_seu1_meta)=gsub(pattern = 'percent.mito',replacement = 'mito',colnames(new_FL_altas_seu1_meta))
  colnames(new_FL_altas_seu1_meta)=gsub(pattern = 'orig.ident',replacement = 'donor',colnames(new_FL_altas_seu1_meta))
  colnames(new_FL_altas_seu1_meta)=gsub(pattern = 'doublets',replacement = 'predicted_doublets',colnames(new_FL_altas_seu1_meta))
  colnames(new_FL_altas_seu1_meta)=gsub(pattern = 'gender',replacement = 'sex',colnames(new_FL_altas_seu1_meta))
  colnames(new_FL_altas_seu1_meta)=gsub(pattern = 'lanes',replacement = 'Sample.lanes',colnames(new_FL_altas_seu1_meta))
  colnames(new_FL_altas_seu1_meta)=gsub(pattern = 'sort.ids',replacement = 'Sort_id',colnames(new_FL_altas_seu1_meta))
  
  new_FL_altas_seu1_meta$age=paste0(round(new_FL_altas_seu1_meta$day/7,1),'PCW')
  new_FL_altas_seu1_meta$age[new_FL_altas_seu1_meta$donor=='F32']='CS22'
  new_FL_altas_seu1_meta$age[new_FL_altas_seu1_meta$donor=='F34']='CS23'
  new_FL_altas_seu1_meta$age[new_FL_altas_seu1_meta$donor=='F35']='CS22'
  new_FL_altas_seu1_meta$resource='Decoding_FL'
  
  new_FL_altas_seu1@meta.data=new_FL_altas_seu1_meta[,c('nCount_RNA', 'nFeature_RNA','mito','donor','Sort_id','Sample.lanes','age','sex','AnnatomicalPart','cell.labels','resource')]
  saveRDS(new_FL_altas_seu1,file = 'NRBC_FL_altas/new_FL_altas_seu1.rds')
  
  
  donor_age_list= unique(new_FL_altas_seu1_meta[,c('donor','age')])[,2]
  names(donor_age_list)=  as.character(unique(new_FL_altas_seu1_meta[,c('donor','age')])[,1])
  unique(FL_altas_seu2$donor)[!unique(FL_altas_seu2$donor) %in% unique(new_FL_altas_seu1_meta$donor)] # F61,F19
  donor_age_list[['F61']]='CS18'
  
  #  F19,只有84个细胞，去除该样本
  FL_altas_seu2=subset(FL_altas_seu2,donor!='F19')
  FL_altas_seu2$age=as.character(donor_age_list[FL_altas_seu2$donor])
  
  
  
  
  # 先整理不在本身注释基因组信息中，因为添加了后缀
  table(rownames(FL_altas_seu2) %in% all_shared_ensembl_id_info$ref_symbol1) #FALSE/TRUE :10870/ 22668  
  table(rownames(FL_altas_seu2) %in% all_shared_ensembl_id_info$ref_symbol)# FALSE/TRUE :24/33514 
  
  symbol_addpostfix=rownames(FL_altas_seu2)[!rownames(FL_altas_seu2) %in% all_shared_ensembl_id_info$ref_symbol]
  symbol_addpostfix=data.frame(strsplit2(symbol_addpostfix,split = "-",perl=F,fix=T))[,1] 
  symbol_addpostfix_inf=all_shared_ensembl_id_info[all_shared_ensembl_id_info$ref_symbol %in% symbol_addpostfix,]
  symbol_addpostfix_inf=symbol_addpostfix_inf[order(symbol_addpostfix_inf$ref_symbol),]
  symbol_addpostfix_inf$number=1:dim(symbol_addpostfix_inf)[1]
  du_symbol_addpostfix=symbol_addpostfix_inf[match(symbol_addpostfix,symbol_addpostfix_inf$ref_symbol)+1,'ref_symbol1']
  names(du_symbol_addpostfix)=paste(symbol_addpostfix,'1',sep = '-')
  
  # 再处理在范围内的symbol, 得到总list,
  uniq_symbols=rownames(FL_altas_seu2)[rownames(FL_altas_seu2) %in% all_shared_ensembl_id_info$ref_symbol]
  new_symbol_list=all_shared_ensembl_id_info[match(uniq_symbols,all_shared_ensembl_id_info$ref_symbol),'ref_symbol1']
  names(new_symbol_list)=uniq_symbols
  new_symbol_list=c(du_symbol_addpostfix,new_symbol_list)
  
  # 将symbol进行映射处理
  du_symbol=as.character(new_symbol_list)[duplicated(as.character(new_symbol_list))];length(du_symbol)
  tmp_assay=GetAssayData(FL_altas_seu2,assay = 'RNA',layer = 'counts')
  rownames(tmp_assay)=as.character(new_symbol_list[rownames(tmp_assay)])
  
  du_tmp_assay=tmp_assay[rownames(tmp_assay) %in% du_symbol,]
  tmp_assay=tmp_assay[!rownames(tmp_assay) %in% du_symbol,]
  
  #du_tmp_assay=aggregate(du_tmp_assay,by=list(rownames(du_tmp_assay)), FUN=sum)  # 计算太久了
  #du_tmp_assay=column_to_rownames(du_tmp_assay,'Group.1')
  #du_tmp_assay=as(as.matrix(du_tmp_assay),'dgCMatrix')
  
  du_tmp_assay_tmp=colSums(du_tmp_assay[rownames(du_tmp_assay) ==du_symbol[1],])
  
  for(symbol in du_symbol[-1]){
    temp=du_tmp_assay[rownames(du_tmp_assay) == symbol,]
    temp=colSums(temp)
    du_tmp_assay_tmp=rbind(du_tmp_assay_tmp,temp)
  }
  rownames(du_tmp_assay_tmp)=du_symbol
  du_tmp_assay_tmp=as(as.matrix(du_tmp_assay_tmp),'dgCMatrix')
  
  #合并得到最后symbol处理矩阵
  tmp_assay=rbind(tmp_assay,du_tmp_assay_tmp)
  FL_altas_seu2=CreateSeuratObject(counts = tmp_assay,min.cells = 10,min.features = 200,meta.data =FL_altas_seu2@meta.data )
  rm(tmp_assay,du_tmp_assay,du_tmp_assay_tmp) ;gc()
  
        
  
  FL_altas_seu2=NormalizeData(FL_altas_seu2) %>% FindVariableFeatures(assay = 'RNA',nfeatures = 3000) %>% ScaleData() %>% RunPCA() %>% RunHarmony(group.by.vars=c('method','donor','Sort_id')) 
  FL_altas_seu2=RunUMAP(FL_altas_seu2,dim=1:30,reduction='harmony',reduction.name='umap1',return.model = T)
  DimPlot(FL_altas_seu2,reduction = 'umap1',group.by = c('anno_lvl_1_final_clean','subcelltype'),cols = cols,raster = T)
  
  p=DimPlot(FL_altas_seu2,reduction = 'umap1',group.by = 'subcelltype',cols = cols,raster = T)+ggtitle('Teichmann_Lab');p
  ggsave(p,filename='NRBC_FL_altas/res/res_pic/FL_Teichmann_Lab_refdata_celltype_information.pdf',width = 8,height = 6,dpi = 300)
  
  FL_altas_seu2=RunUMAP(FL_altas_seu2,dim=1:20,reduction='harmony',reduction.name='umap2')
  DimPlot(FL_altas_seu2,reduction = 'umap2',group.by = c('anno_lvl_1_final_clean','subcelltype'),cols = cols)
  
  saveRDS(FL_altas_seu2,file = 'NRBC_FL_altas/FL_altas_seu2.rds')
  
  table(rownames(new_FL_altas_seu1) %in% all_shared_ensembl_id_info$ref_symbol1) #FALSE/TRUE :28/27052 
  table(rownames(new_FL_altas_seu1) %in% all_shared_ensembl_id_info$ref_symbol)# FALSE/TRUE :27110 /19970 ,不需要调整
  
  new_FL_altas_seu1=NormalizeData(new_FL_altas_seu1) %>% FindVariableFeatures(nfeatures =3000 ) %>% ScaleData() %>% RunPCA()
  anchors=FindTransferAnchors(reference =FL_altas_seu2,query = new_FL_altas_seu1,dims = 1:30,reference.reduction = 'pca')
  new_FL_altas_seu1=MapQuery(anchorset =anchors,query =new_FL_altas_seu1 ,reference = FL_altas_seu2,refdata =  list(celltype = "subcelltype"), reference.reduction = "pca", reduction.model = "umap1" )
  predictions=TransferData(anchorset =anchors,refdata = FL_altas_seu2$anno_lvl_2_final_clean,dims = 1:30 )
  new_FL_altas_seu1$anno_lvl_2_final_clean=predictions$predicted.id
  colnames(new_FL_altas_seu1@meta.data)=gsub(pattern='predicted.celltype',replacement = 'subcelltype',colnames(new_FL_altas_seu1@meta.data))
  
  DimPlot(new_FL_altas_seu1,reduction = 'ref.umap',group.by = 'subcelltype',cols = cols)
  rm(anchors)
  p=DimPlot(new_FL_altas_seu1,reduction = 'umap',cols = cols,group.by = 'donor')+ggtitle('E-MTAB-7407');p
  ggsave(p,filename='NRBC_FL_altas/res/res_pic/FL_E-MTAB-7407_sample_information.pdf',width = 6,height = 6)
  
  new_FL_altas_seu1[['umap1']]=new_FL_altas_seu1[['ref.umap']]
  new_FL_altas_seu1[['ref.pca']]=NULL
  new_FL_altas_seu1[['harmony']]=NULL
  new_FL_altas_seu1[['umap2']]=NULL
  new_FL_altas_seu1[['ref.umap']]=NULL
  new_FL_altas_seu1$method='unknnow'
  saveRDS(new_FL_altas_seu1,file = 'NRBC_FL_altas/new_FL_altas_seu1.rds')
  
  FL_altas_seu2[['umap2']]=NULL
  new_FL_altas_seu1@meta.data=new_FL_altas_seu1@meta.data[,-grep(pattern = 'subcelltype.score', colnames(new_FL_altas_seu1@meta.data))]
}else{
  FL_altas_seu2=readRDS('NRBC_FL_altas/FL_altas_seu2.rds')
  new_FL_altas_seu1=readRDS('NRBC_FL_altas/new_FL_altas_seu1.rds') # E-MTAB-7407中new 部分，需要注意的是 E-MTAB-7407采用的注释基因组和Teichmann lab采用的基因组注释不同
  table(rownames(FL_altas_seu2) %in% rownames(new_FL_altas_seu1))
  table(rownames(FL_altas_seu2) %in% all_shared_ensembl_id_info$ref_symbol1)# FALSE/TRUE : 6/25333 
  table(rownames(FL_altas_seu2) %in% all_shared_ensembl_id_info$ref_symbol) #FALSE/TRUE :6113/19226 
  table(rownames(new_FL_altas_seu1) %in% all_shared_ensembl_id_info$ref_symbol1) #FALSE/TRUE :28/27052 
  table(rownames(new_FL_altas_seu1) %in% all_shared_ensembl_id_info$ref_symbol)# FALSE/TRUE :27110 /19970 
  
  
}

####################################################################################################################################
#----------------------------------CS14-15WPC,fetus,GSE162950---------------------------------------------------#
####################################################################################################################################

if(F){
  # 分选方式CD34+, 无CD45+或者CD45—
  FL_df=read.csv(gzfile('ref_data/ref_scRNAseq_data/GSE162950/GSM4968839_Liver-11wk-569.csv.gz'),header = T); rownames(FL_df)=FL_df$X;FL_df=FL_df[,-1]
  FL_df1=read.csv(gzfile('ref_data/ref_scRNAseq_data/GSE162950/GSM4968840_Liver-15wk-101.csv.gz'),header = T);rownames(FL_df1)=FL_df1$X;FL_df1=FL_df1[,-1]
  FL_df2=read.csv(gzfile('ref_data/ref_scRNAseq_data/GSE162950/GSM4968841_Liver-4wk-658.csv.gz'),header = T); rownames(FL_df2)=FL_df2$X;FL_df2=FL_df2[,-1]
  FL_df3=read.csv(gzfile('ref_data/ref_scRNAseq_data/GSE162950/GSM4968842_Liver-5wk-575.csv.gz'),header = T); rownames(FL_df3)=FL_df3$X;FL_df3=FL_df3[,-1]
  FL_df4=read.csv(gzfile('ref_data/ref_scRNAseq_data/GSE162950/GSM4968843_Liver-6wk-563.csv.gz'),header = T); rownames(FL_df4)=FL_df4$X;FL_df4=FL_df4[,-1]
  FL_df5=read.csv(gzfile('ref_data/ref_scRNAseq_data/GSE162950/GSM4968844_Liver-8wk-553.csv.gz'),header = T); rownames(FL_df5)=FL_df5$X;FL_df5=FL_df5[,-1]
  
  seu=list()
  seu[[1]]=CreateSeuratObject(counts = FL_df,assay = 'RNA',project = 'FL_11wk', min.features=200,min.cells=10)
  seu[[2]]=CreateSeuratObject(counts = FL_df1,assay = 'RNA',project = 'FL_15wk',min.features=200,min.cells=10)
  seu[[3]]=CreateSeuratObject(counts = FL_df2,assay = 'RNA',project = 'FL_4wk',min.features=200 ,min.cells=10)
  seu[[4]]=CreateSeuratObject(counts = FL_df3,assay = 'RNA',project = 'FL_5wk',min.features=200 ,min.cells=10)
  seu[[5]]=CreateSeuratObject(counts = FL_df4,assay = 'RNA',project = 'FL_6wk',min.features=200 ,min.cells=10)
  seu[[6]]=CreateSeuratObject(counts = FL_df5,assay = 'RNA',project = 'FL_8wk',min.features=200 ,min.cells=10)
  
  rm(FL_df,FL_df1,FL_df2,FL_df3,FL_df4,FL_df5)
  
  FL_seu1=merge(seu[[1]],c(seu[[2]],seu[[3]],seu[[4]],seu[[5]],seu[[6]]),projects='FL')
  FL_seu1$stage=rep(c('11PCW','15PCW','CS14_4PCW','CS15_5PCW','CS17_6PCW','8PCW'),c(dim(seu[[1]])[2],dim(seu[[2]])[2],dim(seu[[3]])[2],dim(seu[[4]])[2],dim(seu[[5]])[2],dim(seu[[6]])[2]))
  FL_seu1$stage=factor(FL_seu1$stage,levels = c('CS14_4PCW','CS15_5PCW','CS17_6PCW','8PCW','11PCW','15PCW'))
  FL_seu1$orig.ident=factor(FL_seu1$orig.ident,levels = c("FL_4wk","FL_5wk","FL_6wk", "FL_8wk","FL_11wk","FL_15wk"))
  rm(seu)
  
  
  #-------------------------- annotation the ensembl info -----------------#
   
  
  head(rownames(FL_seu1))
  table(rownames(FL_seu1) %in% rownames(all_shared_ensembl_id_info))# TRUE:20288 
  
  FL_seu1[['RNA']]=JoinLayers(FL_seu1[['RNA']])
  FL_seu_df=GetAssayData(FL_seu1,assay = 'RNA',layer = 'counts')
  rownames(FL_seu_df)=all_shared_ensembl_id_info[rownames(FL_seu_df),'ref_symbol1']
  table(duplicated(rownames(FL_seu_df)))
  rownames(FL_seu_df)[duplicated(rownames(FL_seu_df))]
  
  du_tmp_assay=FL_seu_df[rownames(FL_seu_df) %in% rownames(FL_seu_df)[duplicated(rownames(FL_seu_df))],]
  tmp_assay=FL_seu_df[!rownames(FL_seu_df) %in% unique(rownames(du_tmp_assay)),]
  du_tmp_assay=aggregate(du_tmp_assay,by=list(rownames(du_tmp_assay)), FUN=sum) 
  library(tibble)
  du_tmp_assay=column_to_rownames(du_tmp_assay,'Group.1')
  du_tmp_assay=as(as.matrix(du_tmp_assay),'dgCMatrix')
  FL_seu_df=rbind(tmp_assay,du_tmp_assay);rm(tmp_assay,du_tmp_assay)
  FL_seu1=CreateSeuratObject(counts = FL_seu_df,assay = 'RNA',project = 'FL',min.features=200 ,min.cells=10,meta.data =FL_seu1@meta.data )
  rm(FL_seu_df,gene_barcode_inf,tmp_assay,du_tmp_assay)
  
  
  table(rownames(FL_seu1) %in% all_shared_ensembl_id_info$ref_symbol1) #FALSE/TRUE : 3854/16405
  
  #--------------------------QC-----------------------------------------#
  
  Idents(FL_seu1)='orig.ident'
  
  all_shared_ensembl_id_info$ref_symbol[grep(pattern = '^MT-',all_shared_ensembl_id_info$ref_symbol1)]
  mt_genes=all_shared_ensembl_id_info$ref_symbol1[grep(pattern = '^MT-',all_shared_ensembl_id_info$ref_symbol1)]
  FL_seu1[["percent.mt"]] <- PercentageFeatureSet(FL_seu1, features = mt_genes)
  # Visualize QC metrics as a violin plot
  VlnPlot(FL_seu1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
  plot1 <- FeatureScatter(FL_seu1, feature1 = "nCount_RNA", feature2 = "percent.mt",split.by = 'orig.ident')
  plot2 <- FeatureScatter(FL_seu1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",split.by = 'orig.ident')
  plot1 + plot2
  
  FL_seu1=subset(FL_seu1,subset = nFeature_RNA >200 & nFeature_RNA <7500 & nCount_RNA <75000 & percent.mt < 50)
  cho_cells=rownames(FL_seu1@meta.data[FL_seu1$orig.ident!='FL_4wk' & FL_seu1$percent.mt<20,])
  cho_cells1=rownames(FL_seu1@meta.data[FL_seu1$orig.ident=='FL_4wk',])
  cho_cells=c(cho_cells1,cho_cells)
  FL_seu1=subset(FL_seu1,cells=cho_cells)
  Idents(FL_seu1)='orig.ident'
  VlnPlot(FL_seu1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
  FL_seu1=NormalizeData(FL_seu1) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData() %>% RunPCA() %>% RunHarmony(group.by.vars='orig.ident')
  FL_seu1<- CellCycleScoring(FL_seu1,s.features =cc.genes$s.genes,g2m.features = cc.genes$g2m.genes,set.ident = T )
  ElbowPlot(FL_seu1,ndims = 50)
  FL_seu1=RunUMAP(FL_seu1,dims = 1:30,reduction.name = 'raw_umap')
  FL_seu1=RunUMAP(FL_seu1,dims = 1:20,reduction ='harmony',reduction.name = 'umap' )
  DimPlot(FL_seu1,group.by = 'orig.ident',cols = cols,reduction = 'raw_umap')+DimPlot(FL_seu1,group.by = 'orig.ident',cols = cols,reduction = 'umap')
  
  
  
  # -------------------------将数据mapping 到FL数据上-------------------------------#

  anchors=FindTransferAnchors(reference =FL_altas_seu1,query = FL_seu1,dims = 1:30,reference.reduction = 'pca')
  FL_seu1=MapQuery(anchorset =anchors,query = FL_seu1,reference = FL_altas_seu1,refdata =  list(celltype = "cell.labels"), reference.reduction = "pca", reduction.model = "umap1" )
  rm(anchors,FL_altas_seu1)
  
  cell_levels[! cell_levels %in% unique(FL_seu1$predicted.celltype)]
  FL_seu1$predicted.celltype=factor(FL_seu1$predicted.celltype,levels = cell_levels[!cell_levels=='VCAM1+ EI macrophage'])
  DimPlot(FL_seu1,group.by = 'predicted.celltype',cols = cols,reduction = 'raw_umap')+DimPlot(FL_seu1,group.by = 'predicted.celltype',cols = cols,reduction = 'ref.umap',label = F)+NoLegend()

  p=DimPlot(FL_seu1,group.by = c( 'orig.ident'),cols = cols,reduction = 'raw_umap')+ggtitle('GSE162950:FL');p
  ggsave(p,filename='NRBC_FL_altas/res/res_pic/FL_GSE162950_donor_iformation_umap.pdf',width = 6,height = 6,dpi = 300)
  
  anchors=FindTransferAnchors(reference =FL_altas_seu2,query = FL_seu1,dims = 1:30,reference.reduction = 'pca')
  FL_seu1=MapQuery(anchorset =anchors,query =  FL_seu1,reference =FL_altas_seu2,refdata = list(celltype='subcelltype'),reference.reduction = 'pca',reduction.model = 'umap1' )
  table(as.character( FL_seu1$predicted.celltype))
  FL_seu1$predicted.celltype=as.character(FL_seu1$predicted.celltype)
  colnames(FL_seu1@meta.data)=gsub(pattern = 'predicted.celltype',replacement = 'subcelltype',colnames(FL_seu1@meta.data))
  predictions=TransferData(anchorset = anchors,refdata =FL_altas_seu2$anno_lvl_2_final_clean,dims = 1:30 )
  FL_seu1$anno_lvl_2_final_clean=predictions$predicted.id
  colnames(FL_seu1@meta.data)=gsub(pattern = 'orig.ident',replacement ='donor' ,colnames(FL_seu1@meta.data))
  FL_seu1[['prediction.score.celltype']]=NULL
  colnames(FL_seu1@meta.data)=gsub(pattern = 'percent.mt',replacement ='mito' ,colnames(FL_seu1@meta.data))
  colnames(FL_seu1@meta.data)=gsub(pattern = 'stage',replacement ='age' ,colnames(FL_seu1@meta.data))
  FL_seu1$resource='embryo_FL'
  FL_seu1[['umap1']]=FL_seu1[['ref.umap']]
  FL_seu1[['ref.umap']]=NULL
  saveRDS(FL_seu1,file = 'NRBC_FL_altas/FL_seu1_embryo.rds')
  
}else{
  FL_seu1=readRDS('NRBC_FL_altas/FL_seu1_embryo.rds')
  table(rownames(FL_seu1) %in% all_shared_ensembl_id_info$ref_symbol1) #FALSE/TRUE :4/27052 
  table(rownames(FL_seu1) %in% all_shared_ensembl_id_info$ref_symbol)# FALSE/TRUE :27110 /19970 
  
}





#############################################################################################################################################
# --------------------------------------merged FL_altas -------------------------------#
############################################################################################################################################



FL_altas_seu=merge(FL_altas_seu2,c(new_FL_altas_seu1,FL_seu1))
FL_altas_seu[['RNA']]=JoinLayers(FL_altas_seu[['RNA']])
FL_altas_seu[['umap1']]=merge(FL_altas_seu2[['umap1']],c(new_FL_altas_seu1[['umap1']],FL_seu1[['umap1']]))
FL_altas_seu@meta.data=FL_altas_seu@meta.data[,1:19]

sort(table(FL_altas_seu$subcelltype)) # 25,1747个细胞
table(FL_altas_seu$anno_lvl_2_final_clean[FL_altas_seu$subcelltype=='PROGENITORS'])
FL_altas_seu$subcelltype[FL_altas_seu$subcelltype=='PROGENITORS']=FL_altas_seu$anno_lvl_2_final_clean[FL_altas_seu$subcelltype=='PROGENITORS']


# CHONDROCYT(软骨)E:1,OSTEOCLAST(破骨细胞):2,TISSUE STROMA: 3,NEUTROPHIL:7,MYE_FIBROBLAST: 21,NEPHRON:27, EPITHELIUM:28,GLIAL :37
# 踢出少于20的细胞以及NEPHRON，CC分析时候，剔除少于200的细胞
table(is.na(FL_altas_seu$subcelltype))# FALSE:251747
table(FL_altas_seu$subcelltype)
FL_altas_seu$subcelltype[FL_altas_seu$subcelltype %in% c('CYCLING_DC','DC2','DC_PROGENITOR')]='DC'
FL_altas_seu$subcelltype[FL_altas_seu$subcelltype %in% c('EOSINOPHIL_BASOPHIL','MAST_CELL')]='EO/BASO/MAST'
FL_altas_seu$subcelltype[FL_altas_seu$subcelltype %in% c('EARLY_MK','LATE_MK')]='MEGAKARYOCYTE'


FL_altas_seu=subset(FL_altas_seu,subcelltype %in% names(table(FL_altas_seu$subcelltype)[table(FL_altas_seu$subcelltype) >50]) )
FL_altas_seu$subcelltype=gsub(pattern = 'YS_Bas',replacement = 'YS_Bas/Poly',FL_altas_seu$subcelltype)
sort(table(FL_altas_seu$subcelltype))

PROGENITORS=c("HSC_MPP", "CMP","CYCLING_MPP","MEMP","CYCLING_MEMP","MEP","GMP", "LMPP_MLP")
subcelltype_levels=c(PROGENITORS,"BFUE/CFUE","ProE","Bas","Poly","Orth","YS_Bas/Poly", "YS_Orth", "MEGAKARYOCYTE","EO/BASO/MAST",
                     "MOP","MONOCYTE","DC","MACROPHAGE","MACROPHAGE_ERY","B CELLS", "ILC","NK/T CELLS","MYELOCYTE",
                     "HEPATOCYTE","ENDOTHELIUM","FIBROBLASTS","SMOOTH MUSCLE","SKELETAL MUSCLE","MYE_FIBROBLAST","MESOTHELIUM","NEURON"
);table(unique(FL_altas_seu$subcelltype) %in%  subcelltype_levels)
FL_altas_seu$subcelltype=factor(FL_altas_seu$subcelltype,levels =subcelltype_levels )

rm(FL_seu1,FL_altas_seu2,new_FL_altas_seu1);gc()

FL_NRBC_altas_metadata=FL_altas_seu@meta.data


#cols=c(brewer.pal(12,"Set3"),brewer.pal(12,"Paired")[-3],brewer.pal(9,"Set1")[-5],brewer.pal(9,"Dark2")[-5])
#cols=unique(cols)
cell_cols=c("#8DD3C7","#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F",
       "#A6CEE3", "#1F78B4", "#33A02C", "#FB9A99", "#E31A1C", "#FF7F00", "#CAB2D6", "#6A3D9A", "#B15928",
       "#FFFF33",  "#F781BF", "#999999",  "#D95F02", "#7570B3", "#E7298A", "#E6AB02", "#A6761D")
DimPlot(FL_altas_seu,group.by = c('resource','donor','subcelltype'),reduction = 'umap1',cols = cols,raster =F)

age_levels=c("CS14_4PCW","CS15_5PCW","CS17_6PCW","CS18", "CS22","CS23",'8PCW', "8.1PCW","9.1PCW","9.7PCW","11PCW","11.4PCW","12PCW","13.9PCW","14.4PCW","15PCW","16PCW","16.3PCW","17PCW")
FL_altas_seu$age=factor(FL_altas_seu$age,levels =age_levels )
FL_altas_seu=CellCycleScoring(FL_altas_seu, s.features =cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = F)

FL_altas_seu$method[FL_altas_seu$resource=='embryo_FL']='3GEX'
FL_altas_seu$method[is.na(FL_altas_seu$method)]='unkown'
FL_altas_seu$Sample.lanes[is.na(FL_altas_seu$Sample.lanes)]=FL_altas_seu$donor[is.na(FL_altas_seu$Sample.lanes)]
saveRDS(FL_altas_seu,file = 'NRBC_FL_altas/tmp_FL_altas_seu.rds')


DimPlot(subset(FL_altas_seu,age %in% c("CS14_4PCW","CS15_5PCW","CS17_6PCW")),reduction='umap1',group.by = 'subcelltype',cols = cell_cols,raster =F)+ggtitle('before CS18')
DimPlot(subset(FL_altas_seu,age %in% c("CS14_4PCW","CS15_5PCW","CS17_6PCW")),reduction='umap1',group.by = 'subcelltype',cols = cell_cols,raster =F,split.by = 'age')+ggtitle('before CS18')
FL_altas_seu$stage='post_CS18'
FL_altas_seu$stage[FL_altas_seu$age %in% c("CS14_4PCW","CS15_5PCW","CS17_6PCW")]='pre_CS18'
FL_altas_seu$stage=factor(FL_altas_seu$stage,levels = c('pre_CS18','post_CS18'))
DimPlot(FL_altas_seu,split.by = 'stage',group.by = 'subcelltype',cols = cols,raster = F)


age_celltype_df=data.frame(table(FL_NRBC_altas_metadata[,c('age','subcelltype')]))
age_celltype_df=age_celltype_df[age_celltype_df$Freq>0,]
age_celltype_df=age_celltype_df[age_celltype_df$subcelltype %in% subcelltype_levels,]
age_celltype_df$subcelltype=factor(age_celltype_df$subcelltype,levels =subcelltype_levels)
age_celltype_df$age=factor(age_celltype_df$age,levels =age_levels)

ggplot(age_celltype_df,aes(x=age,fill=subcelltype,y=Freq))+geom_bar(stat ='identity',position = 'fill' )+theme_classic()+scale_fill_manual(values =  cols)+
  theme(axis.text.x =  element_text(angle = 30,vjust = 0.85,hjust = 0.75),axis.text = element_text(face = 'bold'))


############################################################################################################################################################################
#----------------------------------------------FL Ery altas--------------------#
############################################################################################################################################################################
FL_altas_Ery_seu=subset(FL_altas_seu, subcelltype %in% c("BFUE/CFUE","ProE","Bas","Poly","Orth","YS_Bas", "YS_Orth"))
FL_altas_Ery_seu$subcelltype=gsub(pattern = 'YS_Bas',replacement = 'YS_Bas/Poly',FL_altas_Ery_seu$subcelltype)
FL_altas_Ery_seu$subcelltype=factor(FL_altas_Ery_seu$subcelltype,levels = c("BFUE/CFUE","ProE","Bas","Poly","Orth","YS_Bas/Poly","YS_Orth"))
FL_altas_Ery_seu$method[is.na(FL_altas_Ery_seu$method)]='unkown'
FL_altas_Ery_seu$method[FL_altas_Ery_seu$resource=='embryo_FL']='3GEX'

FL_altas_Ery_seu= RunPCA(FL_altas_Ery_seu) %>% RunHarmony(group.by.vars=c('resource','method','donor','Sample.lanes')) %>% RunUMAP(reduction = 'harmony',dim=1:15)

DimPlot(FL_altas_Ery_seu,reduction = 'umap',group.by = 'subcelltype',cols = cols)

FL_altas_Ery_seu$id=paste(FL_altas_Ery_seu$age,FL_altas_Ery_seu$donor,sep = '_')
FL_altas_Ery_seu$id[FL_altas_Ery_seu$id=='CS14_4PCW_FL_4wk']='CS14_4PCW_FL1'
FL_altas_Ery_seu$id[FL_altas_Ery_seu$id=='CS15_5PCW_FL_5wk']='CS15_5PCW_FL2'
FL_altas_Ery_seu$id[FL_altas_Ery_seu$id=='CS17_6PCW_FL_6wk']='CS17_6PCW_FL3'
FL_altas_Ery_seu$id[FL_altas_Ery_seu$id=='8PCW_FL_8wk']='8PCW_FL4'
FL_altas_Ery_seu$id[FL_altas_Ery_seu$id=='11PCW_FL_11wk']='11PCW_FL5'
FL_altas_Ery_seu$id[FL_altas_Ery_seu$id=='15PCW_FL_15wk']='15PCW_FL6'
FL_altas_Ery_seu$id=factor(FL_altas_Ery_seu$id,levels = c("CS14_4PCW_FL1", "CS15_5PCW_FL2" ,"CS17_6PCW_FL3","CS18_F61", "8PCW_FL4",
                                                          "CS22_F32" ,"CS22_F35" ,"CS23_F34","8.1PCW_F16","9.1PCW_F17","9.7PCW_F22","9.7PCW_F33",
                                                          "11PCW_FL5","11.4PCW_F23","12PCW_F38","13.9PCW_F45","14.4PCW_F30","15PCW_FL6",
                                                          "16PCW_F41","16.3PCW_F21" ,"17PCW_F29" ))

FL_altas_Ery_seu$Sample.lanes[is.na(FL_altas_Ery_seu$Sample.lanes)]=FL_altas_Ery_seu$id[is.na(FL_altas_Ery_seu$Sample.lanes)]

FeaturePlot(subset(FL_altas_Ery_seu,age %in% c('CS14_4PCW', 'CS15_5PCW', 'CS17_6PCW','CS18')),features = 'HBE1',reduction = 'umap',split.by = 'age')/
  DimPlot(subset(FL_altas_Ery_seu,age %in% c('CS14_4PCW', 'CS15_5PCW', 'CS17_6PCW','CS18')),cols = cell_cols,group.by = 'subcelltype',reduction = 'umap',split.by = 'age')
FL_altas_Ery_seu$subcelltype[FL_altas_Ery_seu$subcelltype %in% c('Poly') & FL_altas_Ery_seu$age %in% c('CS14_4PCW','CS15_5PCW','CS17_6PCW')]='YS_Bas/Poly'


sub_F61_Ery_seu =subset(FL_altas_Ery_seu,id %in% c('CS18_F61'))  # immune_paper,来自同一篇文章 
sub_F61_Ery_seu=NormalizeData(sub_F61_Ery_seu) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData() %>% RunPCA() %>% RunUMAP(reduction.name = 'umap1',dims = 1:15,return.model = T)
DimPlot(sub_F61_Ery_seu,reduction = 'umap1',group.by = c('subcelltype','id'),cols = cols)

celltype_features=c('KIT',"TFRC","GYPA",'HBE1','NCL')
VlnPlot(sub_F61_Ery_seu,group.by = 'subcelltype',features = celltype_features ,stack = T)+NoLegend()

p=DimPlot(sub_F61_Ery_seu,reduction = 'umap1',group.by = 'subcelltype',cols = cols)+ggtitle('FL_F61.celltype')+FeaturePlot(sub_F61_Ery_seu,features = c('HBE1'),reduction = 'umap1')
p
ggsave(p,filename='NRBC_FL_altas/res_pic/ref_FL_F61_celltype_HBE1_umap_info.pdf',width = 10,height = 5)

if(F){
  sub_F61_Ery_seu=FindNeighbors(sub_F61_Ery_seu,dims=1:20);sub_F61_Ery_seu=FindClusters(sub_F61_Ery_seu,resolution = 0.5)
  
  # cluster 9 是否是EMP？ 
  DimPlot(sub_F61_Ery_seu,reduction = 'umap1',group.by = c('seurat_clusters','subcelltype'),cols = cols)/
    VlnPlot(sub_F61_Ery_seu,group.by = 'seurat_clusters',features = celltype_features ,stack = T)+NoLegend()
  
  
  sub_F61_Ery_seu$subcelltype[sub_F61_Ery_seu$seurat_clusters %in% c(8,4,2,7) & sub_F61_Ery_seu$subcelltype=='Bas']='YS_Bas/Poly'
  sub_F61_Ery_seu$subcelltype[sub_F61_Ery_seu$seurat_clusters %in% c(8,4,2,7) & sub_F61_Ery_seu$subcelltype=='Poly']='YS_Bas/Poly'
  sub_F61_Ery_seu$subcelltype[sub_F61_Ery_seu$seurat_clusters %in% c(8,4,2,7) & sub_F61_Ery_seu$subcelltype=='Orth']='YS_Orth'
  
}


saveRDS(sub_F61_Ery_seu,file = 'NRBC_FL_altas/tmp_sub_F61_Ery_seu_umap_model.rds')

sub_FL_Ery_seu =subset(FL_altas_Ery_seu,donor %in% unique(FL_altas_Ery_seu$donor)[-3])                                                         
sub_FL_Ery_seu=NormalizeData(sub_FL_Ery_seu) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% RunUMAP(reduction.name = 'umap1',dims = 1:15,return.model = T)
anchors=FindTransferAnchors(reference =sub_F61_Ery_seu,query = sub_FL_Ery_seu,reference.reduction = 'pca' )
sub_FL_Ery_seu=MapQuery(anchorset =anchors,query =sub_FL_Ery_seu,reference = sub_F61_Ery_seu, reference.reduction = 'pca',reduction.model = 'umap1' )
DimPlot(sub_FL_Ery_seu,reduction = 'ref.umap',group.by = 'subcelltype',cols = cols)+FeaturePlot(sub_FL_Ery_seu,features = c('HBE1','HBZ'),reduction = 'ref.umap')
rm(anchors)

FL_altas_Ery_seu[['umap2']]=merge(sub_F61_Ery_seu[['umap1']],sub_FL_Ery_seu[['ref.umap']])
p=DimPlot(FL_altas_Ery_seu,reduction = 'umap2',group.by = c('subcelltype'),cols = cols)+ggtitle('FL: 63,606');p
ggsave(p,filename='NRBC_FL_altas/res_pic/FL_Ery_celltype_umap.pdf',width=6,height = 6,dpi = 300)

FL_altas_Ery_seu$sourceid=FL_altas_Ery_seu$resource
FL_altas_Ery_seu$sourceid[FL_altas_Ery_seu$resource=='Decoding_FL']='E-MTAB-7407'
FL_altas_Ery_seu$sourceid[FL_altas_Ery_seu$resource=='embryo_FL']='GSE162950'
FL_altas_Ery_seu$sourceid[FL_altas_Ery_seu$resource=='immune_paper']='Teichmann Lab'

p=DimPlot(FL_altas_Ery_seu,reduction = 'umap2',group.by = c('sourceid','id'),cols = cols);p
ggsave(p,filename='NRBC_FL_altas/res_pic/FL_Ery_source_sample_info_umap.pdf',width=14,height = 6,dpi = 300)

FeaturePlot(FL_altas_Ery_seu,reduction = 'umap2',features = 'HBE1',ncol = 23)

Ery_age_celltype_df=data.frame(table(FL_altas_Ery_seu@meta.data[,c('id','subcelltype')]))
Ery_age_celltype_df=Ery_age_celltype_df[Ery_age_celltype_df$Freq>0,]
celltype_age_df=data.frame(table(FL_altas_Ery_seu$id));colnames(celltype_age_df)=c('id','count')

p=ggplot(Ery_age_celltype_df,aes(x=id,fill=subcelltype,y=Freq))+geom_bar(stat ='identity',position = 'fill' )+theme_classic()+scale_fill_manual(values =  cell_cols)+
  theme(axis.text.x =  element_text(angle = 30,vjust = 0.85,hjust = 0.75),axis.text = element_text(face = 'bold'))+
  geom_label(data = celltype_age_df,aes(x=id,y=1.06,label=count,fill=NULL),vjust ='top',show.legend = F,label.size = NA)+scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1))
  
ggsave(p,filename='NRBC_FL_altas/res_pic/FL_celltype_barplot.pdf',width = 12,height = 6)

FL_altas_Ery_seu$orign='definitive'
FL_altas_Ery_seu$orign[grep('YS',FL_altas_Ery_seu$subcelltype)]='primitive'

table(FL_altas_Ery_seu$Sort_id)
FL_altas_Ery_seu$Sort_id[FL_altas_Ery_seu$Sort_id=='CD45-']='CD45N'
FL_altas_Ery_seu$Sort_id[FL_altas_Ery_seu$Sort_id=='CD45+']='CD45P'
FL_altas_Ery_seu$Sort_id[FL_altas_Ery_seu$Sort_id=='TOT']='Total'
FL_altas_Ery_seu$Sort_id[FL_altas_Ery_seu$resource=='embryo_FL']='Total'


temp_df1=data.frame(table(FL_altas_Ery_seu@meta.data[FL_altas_Ery_seu$orign=='primitive',c('id','Sort_id')]))
temp_df1=temp_df1[temp_df1$Freq >0,]
temp_df1$orign='primitive'
temp_df=data.frame(table(FL_altas_Ery_seu@meta.data[FL_altas_Ery_seu$orign=='definitive',c('id','Sort_id')]))
temp_df=temp_df[temp_df$Freq >0,]
temp_df$orign='definitive'
temp_df=rbind(temp_df1,temp_df)
temp_df$group=paste(temp_df$orign,temp_df$Sort_id,sep=':')

p=ggplot(temp_df,aes(x=id,y=Freq,color=Sort_id,linetype=orign,group=group))+geom_point()+geom_line(size=0.8)+theme_classic()+
  RotatedAxis()+scale_color_manual(values = cols[-2])+theme(text = element_text(face = 'bold')) +scale_linetype_manual(values = c(1,2))#  log2(Freq) 更能反应波动
p
ggsave(as.ggplot(p),filename = 'NRBC_FL_altas/res_pic/FL_source_count_poin_line.pdf',width = 6,height = 4,dpi = 300)


rm(sub_F61_Ery_seu,sub_FL_Ery_seu)
saveRDS(FL_altas_Ery_seu,file = 'NRBC_FL_altas/tmp_FL_altas_Ery_seu.rds')

# CS 16之前，HSPC暂未定殖在FL，prim_stage, YS_Ery+ definitive Ery: CS17-CS18, mix_stage,def_stage:CS22- ,definitive Ery;
FL_altas_seu$stage='def_stage'
FL_altas_seu$stage[FL_altas_seu$age %in% c("CS17_6PCW", "CS18")]='mix_stage'
FL_altas_seu$stage[FL_altas_seu$age %in% c("CS14_4PCW","CS15_5PCW")]='prim_stage'
FL_altas_seu$stage=factor(FL_altas_seu$stage,levels = c('prim_stage','mix_stage','def_stage'))

DimPlot(FL_altas_seu,reduction = 'umap1',group.by = 'subcelltype',cols = cell_cols)

cell_levels2=gsub(pattern = 'YS_Bas/Poly/Poly',replacement = 'YS_Bas/Poly',levels(FL_altas_seu@meta.data$subcelltype))
FL_altas_seu@meta.data$subcelltype=as.character(FL_altas_seu@meta.data$subcelltype)
FL_altas_seu@meta.data[rownames(FL_altas_Ery_seu@meta.data[FL_altas_Ery_seu$subcelltype=='YS_Bas/Poly',]),'subcelltype']='YS_Bas/Poly'
FL_altas_seu@meta.data$subcelltype=factor(FL_altas_seu@meta.data$subcelltype,levels =cell_levels2 )
table(FL_altas_seu$subcelltype)

temp_df=data.frame(table(FL_altas_Ery_seu@meta.data[FL_altas_Ery_seu$Sort_id=='CD45P',c('id','subcelltype','Sort_id')]))
temp_df=temp_df[temp_df$Freq >0,]
celltype_count_df=data.frame(table(FL_altas_Ery_seu@meta.data[FL_altas_Ery_seu$Sort_id=='CD45P',c('id')]))
celltype_count_df=celltype_count_df[celltype_count_df$Freq >0,];colnames(celltype_count_df)=c('id','count')

#  可以选择细胞数目大于50
p=ggplot(temp_df,aes(x=id,fill=subcelltype,y=Freq))+geom_bar(stat ='identity',position = 'fill' )+theme(text = element_text(face = 'bold'))+
  theme_classic()+scale_fill_manual(values =  col)+RotatedAxis()+ggtitle('FBM nRBC: CD45P sorting')+scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1))+
  geom_label(data = celltype_count_df,aes(x=id,y=1.16,label=count,fill=NULL),vjust ='top',show.legend = F,label.size = NA)
p
ggsave(p,filename='NRBC_FL_altas/res_pic/FL_NRBC_CD45P_barplot.pdf',width=10,height = 6,dpi = 300)


saveRDS(FL_altas_seu,file = 'NRBC_FL_altas/tmp_FL_altas_seu.rds')

############################################################################################################################################################################
#-----------------------------------------cellchat analysis--------------------------------#
############################################################################################################################################################################
future.seed=TRUE
library(CellChat)
library(Seurat)
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB # simply use the default CellChatDB,
options(future.globals.maxSize=80000 * 1024 ^ 2 )
cho_cells=sample(colnames(FL_altas_seu),size=50000)
FL_altas_seu=subset(FL_altas_seu,cells=cho_cells)
FL_altas_seu=subset(FL_altas_seu, subcelltype %in% names(table(FL_altas_seu$subcelltype))[table(FL_altas_seu$subcelltype) >20])
FL_altas_seu$subcelltype=droplevels(FL_altas_seu$subcelltype, exclude = setdiff(levels(FL_altas_seu$subcelltype),unique(FL_altas_seu$subcelltype)))
FL_altas_seu=NormalizeData(FL_altas_seu)
cellchat=createCellChat(FL_altas_seu,assay = 'RNA',group.by='subcelltype')
rm(FL_altas_seu)

CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat@DB=CellChatDB.use
cellchat=subsetData(cellchat)
future::plan("multisession", workers = 4) # do parallel
cellchat=identifyOverExpressedGenes(cellchat)
cellchat=identifyOverExpressedInteractions(cellchat)
cellchat=projectData(cellchat,PPI.human)
celllchat= computeCommunProb(cellchat, type = "triMean",raw.use = FALSE,population.size = TRUE)
#celllchat=computeCommunProb(cellchat)
cellchat=filterCommunication(cellchat, min.cells = 20)
saveRDS(cellchat,file='NRBC_FL_altas/FL_subcelltype_cellchat.rds')
cellchat=computeCommunProbPathway(cellchat)
cellchat=aggregateNet(cellchat)
cellchat=netAnalysis_computeCentrality(cellchat)

saveRDS(cellchat,file='NRBC_FL_altas/FL_subcelltype_cellchat.rds')






