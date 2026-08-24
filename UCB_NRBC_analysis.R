
source('/home/gibh/2021_NRBC_chlyu/zx_lab_NRBC/scripts/scRNAseq_pipline/scRNAseq_analysis_model.R')
set.seed(123)
# 设置全局变量内存大小，默认是500M，我们数据超过，需要重新设置
options(future.globals.maxSize=50000 * 1024 ^ 2 )
library(future)
#plan(multicore, workers = 4) # not supported in Rstudio because it is considered unstable

library(corrplot)
library(clusterProfiler)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(limma)# 使用strsplit2函数
#  变量 col，颜色
library(RColorBrewer)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(SeuratDisk)
library(tibble)
library(pheatmap)
library(harmony)
library(dplyr)

cols=c(brewer.pal(12,"Set3"),brewer.pal(6,"PiYG"),brewer.pal(6,"BrBG"),brewer.pal(8,"Set2"),
       brewer.pal(12,"Set3"),brewer.pal(8,"Pastel2"),brewer.pal(9,"Pastel1"),brewer.pal(8,"Accent"))
col=unique(cols)[-14]
setwd('/home/gibh/2021_NRBC_chlyu/')


########################################################################################################################
#--------------------------------------------prepare data---------------------------------------------#
########################################################################################################################
# our lab data, 3 samples UBC ~20 week, our lab data 
if(!file.exists('/home/gibh/2021_NRBC_chlyu/results/zxlab_NRBC_vis_NRBC.Rdata')){
  load('/home/1.gibh/2021_NRBC_chlyu/results/zxlab_NRBC_vis_NRBC.Rdata',verbose = T) 
  unique(an_vis_NRBC$orig.ident)
  an_vis_NRBC@meta.data=an_vis_NRBC@meta.data[,c(1:8,15:16)]
  an_vis_NRBC@meta.data$orig.ident=paste0('fetal_',an_vis_NRBC@meta.data$orig.ident)
  an_vis_NRBC$type='fetal_UCB'
  an_vis_NRBC$stage='22WPC'
  an_vis_NRBC$stage[an_vis_NRBC$orig.ident=='fetal_NRBC3']='18WPC'
  an_vis_NRBC$ref_dataid='zx_lab'
  
  colnames(an_vis_NRBC@meta.data)= gsub(colnames(an_vis_NRBC@meta.data),pattern = '^celltype',replacement = 'cell.labels')
  colnames(an_vis_NRBC@meta.data)= gsub(colnames(an_vis_NRBC@meta.data),pattern = 'age_PCW',replacement = 'stage')
  
  save(an_vis_NRBC,file = '/home/gibh/2021_NRBC_chlyu/results/zxlab_NRBC_vis_NRBC.Rdata')
  fetal_UCB=an_vis_NRBC;rm(an_vis_NRBC)
  DimPlot(fetal_UCB,reduction = 'umap',cols = cols)
  
  fetal_UCB$cell.labels=as.character(fetal_UCB$cell.labels)
  fetal_UCB$cell.labels[fetal_UCB$cell.labels %in% c( 'lBas', 'ePoly', 'lPoly','ProE/eBas')]='Poly'
  fetal_UCB$cell.labels[fetal_UCB$cell.labels %in% c( 'CFUE', 'BFUE')]='Bas'
  fetal_UCB$cell.labels=factor(fetal_UCB$cell.labels,levels = c('Bas','Poly',"eOrth","mOrth","lOrth" ))
  Idents(fetal_UCB)='cell.labels'
  
  FeaturePlot(fetal_UCB,reduction = 'umap',features = c('KIT',"GYPA","TFRC","MALAT1",'NCL','CD63'))
  DimPlot(fetal_UCB,reduction = 'umap',cols = cols,group.by =c( 'cell.labels','orig.ident'))
  saveRDS(fetal_UCB,file = 'results/zxlab_NRBC.rds')
  
}else{
  fetal_UCB=readRDS('results/zxlab_NRBC.rds')
}


# Nature immunology : immunology erythroid precursors, 所提供 的 ery 没有分型，不知道是否有MEMP
if(!file.exists('ref_data/primitive_NI_Ery.RDS')){
  if(!file.exists('ref_data/primitive_NI_Ery.RDS')){
    # Yolk sac, data source :immunology erythroid precursors
    CS10_YS=LoadLoom('/home/gibh/2021_NRBC_chlyu/ref_data/ref_scRNAseq_data/new_blood_data/matrix_loom/CS10_YS.loom');CS10_YS$orig.ident='CS10_YS';head(CS10_YS) # 与zx_lab feuature inf 高度重合，
    CS11_YS=LoadLoom('/home/gibh/2021_NRBC_chlyu/ref_data/ref_scRNAseq_data/new_blood_data/matrix_loom/CS11_YS.loom');CS11_YS$orig.ident='CS11_YS';head(CS11_YS)
    CS15_YS=LoadLoom('/home/gibh/2021_NRBC_chlyu/ref_data/ref_scRNAseq_data/new_blood_data/matrix_loom/CS15_YS.loom');CS15_YS$orig.ident='CS15_YS';head(CS15_YS)
    CS11_YS=RenameCells(CS11_YS,new.names = gsub(pattern = 'YS_R2',replacement = 'YS',rownames(CS11_YS@meta.data)))
    CS15_YS=RenameCells(CS15_YS,new.names = gsub(pattern = 'YS_R1',replacement = 'YS',rownames(CS15_YS@meta.data)))
    
    fetal_YS=merge(CS10_YS,c(CS11_YS,CS15_YS));
    fetal_YS$type='YS'
    
    dim(fetal_YS)
    rm(list = c('CS10_YS','CS11_YS','CS15_YS'))
    if(! exists('ref_data/NI_fetal_YS_test.rds')){
      fetal_YS_test=fetal_YS
      fetal_YS_test[['percent.mt']] <- PercentageFeatureSet(fetal_YS_test,pattern = '^MT-' )
      fetal_YS_test$GS_score=fetal_YS_test$G2M.Score-fetal_YS_test$S.Score
      
      fetal_YS_test[['RNA']]=split(fetal_YS_test[['RNA']],f = fetal_YS_test$orig.ident) 
      fetal_YS_test=NormalizeData(fetal_YS_test)
      fetal_YS_test=FindVariableFeatures(fetal_YS_test)
      fetal_YS_test=ScaleData(fetal_YS_test,vars.to.regress = c('percent.mt','nCount_RNA',''))
      fetal_YS_test=RunPCA(fetal_YS_test)
      ElbowPlot(fetal_YS_test,ndims = 50)
      fetal_YS_test=IntegrateLayers(fetal_YS_test,method = CCAIntegration,orig.reduction = 'pca',new.reduction='cca')
      fetal_YS_test=RunUMAP(fetal_YS_test,dims = 1:30,reduction = 'cca',reduction.name = 'cca_umap')
      DimPlot(fetal_YS_test,cols = col,group.by = 'orig.ident',reduction = 'cca_umap')+FeaturePlot(fetal_YS_test,features = 'GYPA',reduction = 'cca_umap')
      
      fetal_YS_test=IntegrateLayers(fetal_YS_test,method = HarmonyIntegration,orig.reduction = 'pca',new.reduction='harmony')
      fetal_YS_test=RunUMAP(fetal_YS_test,dims = 1:30,reduction = 'harmony')
      # 有双细胞存在，需要过滤去除
      
      fetal_YS_test$celltype='others'
      fetal_YS_test@meta.data[rownames(all_NI_Ery_meta)[rownames(all_NI_Ery_meta) %in% rownames(fetal_YS_test@meta.data) ],'celltype']='Ery'
      DimPlot(fetal_YS_test,cols = col,group.by = c('orig.ident','celltype'))+FeaturePlot(fetal_YS_test,features = 'GYPA')
      
      ggplot(data.frame(table(fetal_YS_test@meta.data[,c('orig.ident','celltype')])),aes(x=orig.ident,y=Freq,fill=celltype))+geom_bar(stat = 'identity',position = 'fill',alpha=0.8)+theme_classic()
      saveRDS(fetal_YS_test,file = 'ref_data/NI_fetal_YS_test.rds')
      
      rm(fetal_YS_test)
    }
    
    #[1] 33538 21140,远多于YS_Ery 细胞数据
    
    pre_UBC1=sRNA_qc_cal_func(inpath = 'ref_data/ref_scRNAseq_data/new_blood_data/pre_UCB_raw_matrix/pre_UBC1/filtered_feature_bc_matrix',proname = 'pre_UCB1')
    pre_UBC2=sRNA_qc_cal_func(inpath = 'ref_data/ref_scRNAseq_data/new_blood_data/pre_UCB_raw_matrix/pre_UBC2/filtered_feature_bc_matrix',proname = 'pre_UCB2')
    pre_UBC3=sRNA_qc_cal_func(inpath = 'ref_data/ref_scRNAseq_data/new_blood_data/pre_UCB_raw_matrix/pre_UBC3/filtered_feature_bc_matrix',proname = 'pre_UCB3')
    pre_UBC1[[1]]=RenameCells(pre_UBC1[[1]],add.cell.id = 'pre_UCB1')
    pre_UBC2[[1]]=RenameCells(pre_UBC2[[1]],add.cell.id = 'pre_UCB2')
    pre_UBC3[[1]]=RenameCells(pre_UBC3[[1]],add.cell.id = 'pre_UCB3')
    pre_UBC=merge(pre_UBC1[[1]],c(pre_UBC2[[1]],pre_UBC3[[1]]))
    pre_UBC[['RNA']]=JoinLayers(pre_UBC[['RNA']])
    pre_UBC$type='pre_UCB'
    
    pre_UBC=RenameCells(pre_UBC,new.names = gsub(pattern = '-',replacement = '_',rownames(pre_UBC@meta.data)))
    dim(pre_UBC) ;rm(list = c('pre_UBC1','pre_UBC2','pre_UBC3'))
    #[1] 19073 16282,feature 个数明显偏少了
    # 细胞数目远远多于6640, 理论上不可能过滤掉大半，存在其他类型的细胞
    
    # ~ 40 week UBC NRBC sample, read the raw data ,data source : immunology erythroid precursor
    UCB1=LoadLoom('/home/gibh/2021_NRBC_chlyu/ref_data/ref_scRNAseq_data/new_blood_data/matrix_loom/UCB_FACS_R1.loom');UCB1$orig.ident='term_UCB1'
    UCB2=LoadLoom('/home/gibh/2021_NRBC_chlyu/ref_data/ref_scRNAseq_data/new_blood_data/matrix_loom/UCB_FACS_R2.loom');UCB2$orig.ident='term_UCB2'
    UCB3=LoadLoom('/home/gibh/2021_NRBC_chlyu/ref_data/ref_scRNAseq_data/new_blood_data/matrix_loom/UCB_FACS_R3.loom');UCB3$orig.ident='term_UCB3'
    term_UCB=merge(UCB1,c(UCB2,UCB3))
    term_UCB$type='term_UCB'
    dim(term_UCB);rm(list = c('UCB1','UCB2','UCB3'))
    #[1] 33538 24057
    primitive_NI=merge(fetal_YS,c(pre_UBC,term_UCB));primitive_NI@meta.data=primitive_NI@meta.data[,c('orig.ident', 'nCount_RNA', 'nFeature_RNA', 'Clusters','type')]
    dim(primitive_NI)
    primitive_NI[['RNA']]=JoinLayers(primitive_NI[['RNA']])
    
    
    primitive_NI_Ery[['percent.mt']] <- PercentageFeatureSet(primitive_NI_Ery,pattern = '^MT-' )
    ribo <- c(grep(pattern = "^RPS", x = rownames(primitive_NI_Ery), value = T),
              grep(pattern = "^RPL", rownames(primitive_NI_Ery), value = T))
    primitive_NI_Ery[['percent.rb']] <- PercentageFeatureSet(primitive_NI_Ery,features = ribo)
    
    
    saveRDS(primitive_NI,file = 'ref_data/primitive_NI.RDS')
  }else{
    primitive_NI=readRDS('ref_data/primitive_NI.RDS')
  }
  
  rm(list = c('fetal_YS','pre_UBC','term_UCB'))
  # -------------------------------------------meta info -------------------------------------------#
  if(!file.exists('ref_data/all_NI_Ery_meta.Rdata')){
    load('ref_data/ref_scRNAseq_data/new_blood_data/Ery.integrated.RData',verbose = T);rm(Ery.integrated_counts_matrix)
    table(Ery.integrated_meta_info$sample)
    NI_Ery_meta=Ery.integrated_meta_info
    dim(NI_Ery_meta)
    #[1] 31121     6
    NI_Ery_meta$stage=NA
    NI_Ery_meta$stage[NI_Ery_meta$sample=='CS10_YS']='4WPC'
    NI_Ery_meta$stage[NI_Ery_meta$sample=='CS11_YS']='4WPC'
    NI_Ery_meta$stage[NI_Ery_meta$sample=='CS15_YS']='5WPC'
    NI_Ery_meta$stage[NI_Ery_meta$sample=='UCB_FACS_R1']='~40WPC'
    NI_Ery_meta$stage[NI_Ery_meta$sample=='UCB_FACS_R2']='~40WPC'
    NI_Ery_meta$stage[NI_Ery_meta$sample=='UCB_FACS_R3']='~40WPC'
    NI_Ery_meta$stage[NI_Ery_meta$sample=='FL_10WPC']='10WPC'
    NI_Ery_meta$stage[NI_Ery_meta$sample=='FL_10.5WPC']='10.5WPC'
    NI_Ery_meta$stage[NI_Ery_meta$sample=='FL_12.5WPC']='12.5WPC'
    rownames(NI_Ery_meta)=paste(NI_Ery_meta$sample,rownames(NI_Ery_meta),sep = ':')
    rownames(NI_Ery_meta)=gsub(pattern = 'YS_',replacement = '',rownames(NI_Ery_meta))
    rownames(NI_Ery_meta)=gsub(pattern = ':UCB_',replacement = ':',rownames(NI_Ery_meta))
    rownames(NI_Ery_meta)=gsub(pattern = '_.$',replacement = 'x',rownames(NI_Ery_meta))
    load('ref_data/ref_scRNAseq_data/new_blood_data/pre_UCB_Ery.Rdata',verbose = T)
    unique(Ery.integrated_meta_info[,3:4])
    #samples stage
    #AAACCCATCGACGACC_1 pre_UCB_FACS_R1 35WPC
    #AAACGCTAGTGGTGGT_2 pre_UCB_FACS_R2 33WPC
    #AAACCCAAGTCGAGGT_3 pre_UCB_FACS_R3 32WPC
    pre_UCB_Ery_meta=Ery.integrated_meta_info
    dim(pre_UCB_Ery_meta)
    #[1] 6640    6
    colnames(pre_UCB_Ery_meta)=c(colnames(pre_UCB_Ery_meta)[1:2],'sample',colnames(pre_UCB_Ery_meta)[4:5],'cluster')
    pre_UCB_Ery_meta$site='UCB'
    pre_UCB_Ery_meta$sample=gsub(pre_UCB_Ery_meta$sample,pattern = '_FACS_R',replacement = '')
    rownames(pre_UCB_Ery_meta)=paste(pre_UCB_Ery_meta$sample,rownames(pre_UCB_Ery_meta),sep = "_")
    rownames(pre_UCB_Ery_meta)=gsub(pattern = '_2',replacement = '_1',rownames(pre_UCB_Ery_meta))
    rownames(pre_UCB_Ery_meta)=gsub(pattern = '_3',replacement = '_1',rownames(pre_UCB_Ery_meta))
    all_NI_Ery_meta=rbind(NI_Ery_meta,pre_UCB_Ery_meta)
    save(all_NI_Ery_meta,file='ref_data/all_NI_Ery_meta.Rdata')
    rm(list = c('NI_Ery_meta','pre_UCB_Ery_meta'))
  }else{load(file = 'ref_data/all_NI_Ery_meta.Rdata',verbose = T)}
  
  primitive_NI_Ery=subset(primitive_NI,cells=rownames(all_NI_Ery_meta[! all_NI_Ery_meta$sample %in% c('FL_10WPC','FL_10.5WPC','FL_12.5WPC'),]));dim(primitive_NI_Ery)
  primitive_NI_Ery@meta.data[,c('Clusters','stage')]=all_NI_Ery_meta[rownames(primitive_NI_Ery@meta.data),c('cluster','stage')]
  primitive_NI_Ery$ref_dataid='Shi_Lilong_Lab'
  primitive_NI_Ery[['RNA']]=JoinLayers(primitive_NI_Ery[['RNA']])
  primitive_NI_Ery[["percent.mt"]] <- PercentageFeatureSet(primitive_NI_Ery, pattern = "^MT-")
  primitive_NI_Ery[["percent.rb"]] <- PercentageFeatureSet(primitive_NI_Ery, pattern = "^RP[SL]")
  
  saveRDS(primitive_NI_Ery,file = 'ref_data/primitive_NI_Ery.RDS')
  rm(primitive_NI);gc()
  
}else{
  primitive_NI_Ery=readRDS('ref_data/primitive_NI_Ery.RDS')
  
} 


pre_UBC=subset(primitive_NI_Ery,type=='pre_UCB')
term_UCB=subset(primitive_NI_Ery,type=='term_UCB')
dim(pre_UBC)
dim(term_UCB)
pre_UCB_qc_pare=c(200,6000,70000,10,60)
pre_UBC=subset(pre_UBC, subset = nFeature_RNA > pre_UCB_qc_pare[1] & nFeature_RNA < pre_UCB_qc_pare[2] & nCount_RNA <pre_UCB_qc_pare[3] & percent.mt < pre_UCB_qc_pare[4] & percent.rb < pre_UCB_qc_pare[5])
dim(pre_UBC)
term_UCB_qc_pare=c(200,7000,70000,25,60)
term_UCB=subset(term_UCB, subset = nFeature_RNA > term_UCB_qc_pare[1] & nFeature_RNA < term_UCB_qc_pare[2] & nCount_RNA <term_UCB_qc_pare[3] & percent.mt < term_UCB_qc_pare[4] & percent.rb < term_UCB_qc_pare[5])
dim(term_UCB)


sapply(colnames(all_shared_ensembl_id_info), function(x){table(rownames(fetal_UCB) %in% all_shared_ensembl_id_info[,x])})
sapply(colnames(all_shared_ensembl_id_info), function(x){table(rownames(pre_UBC) %in% all_shared_ensembl_id_info[,x])})
sapply(colnames(all_shared_ensembl_id_info), function(x){table(rownames(term_UCB) %in% all_shared_ensembl_id_info[,x])}) 

# 属于相同symbol注释结果信息：zxlab_symbol:ref_symbol

UCB_NRBC_altas=merge(pre_UBC,term_UCB)
UCB_NRBC_altas[['RNA']]=JoinLayers(UCB_NRBC_altas[['RNA']])

UCB_NRBC_altas=merge(UCB_NRBC_altas,fetal_UCB)
UCB_NRBC_altas[['RNA']]=JoinLayers(UCB_NRBC_altas[['RNA']])

UCB_NRBC_altas[['prediction.score.celltype.l1']]=NULL
UCB_NRBC_altas[['integrated']]=NULL

rm(fetal_UCB,pre_UCB,term_UCB);gc()

saveRDS(UCB_NRBC_altas,file = 'NRBC_UCB_altas/UCB_NRBC_20260515.rds')

dir.create('NRBC_UCB_altas/res_pic')

p=DimPlot(UCB_NRBC_altas,reduction = 'umap',group.by = c('celltype','type'),cols = cols,raster = F);p
ggsave(p,file='NRBC_UCB_altas/res_pic/UCB_nRBC_celltype_stage_umap.pdf',width = 10,height = 4)


# 先整理不在本身注释基因组信息中，因为添加了后缀
symbol_addpostfix=rownames(UCB_NRBC_altas)[!rownames(UCB_NRBC_altas) %in% all_shared_ensembl_id_info$ref_symbol]
symbol_addpostfix=data.frame(strsplit2(symbol_addpostfix,split = '.',fix=T))[,1] 
symbol_addpostfix_inf=all_shared_ensembl_id_info[all_shared_ensembl_id_info$ref_symbol %in% symbol_addpostfix,]
symbol_addpostfix_inf=symbol_addpostfix_inf[order(symbol_addpostfix_inf$ref_symbol),]
symbol_addpostfix_inf$number=1:dim(symbol_addpostfix_inf)[1]
du_symbol_addpostfix=symbol_addpostfix_inf[match(symbol_addpostfix,symbol_addpostfix_inf$ref_symbol)+1,'ref_symbol1']
names(du_symbol_addpostfix)=paste(symbol_addpostfix,'1',sep = '.')

# 再处理在范围内的symbol
uniq_symbols=rownames(UCB_NRBC_altas)[rownames(UCB_NRBC_altas) %in% all_shared_ensembl_id_info$ref_symbol]
new_symbol_list=all_shared_ensembl_id_info[match(uniq_symbols,all_shared_ensembl_id_info$ref_symbol),'ref_symbol1']
names(new_symbol_list)=uniq_symbols
new_symbol_list=c(du_symbol_addpostfix,new_symbol_list)

# 得到总list
du_symbol=as.character(new_symbol_list)[duplicated(as.character(new_symbol_list))];length(du_symbol)
tmp_assay=GetAssayData(UCB_NRBC_altas,assay = 'RNA',layer = 'counts')
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
new_UCB_NRBC_altas=CreateSeuratObject(counts = tmp_assay,min.cells = 10,min.features = 200,meta.data =UCB_NRBC_altas@meta.data )
rm(tmp_assay,du_tmp_assay)
saveRDS(new_UCB_NRBC_altas,file = 'NRBC_UCB_altas/dealt_symbol_UCB_NRBC.rds')

new_UCB_NRBC_altas=NormalizeData(new_UCB_NRBC_altas) %>%FindVariableFeatures() %>%RunPCA()

############################################################################################################################################################
#---------------------term NRBC as reference------------------#
############################################################################################################################################################
term_UCB_seu=subset(new_UCB_NRBC_altas,type=='term_UCB')
term_UCB_seu=NormalizeData(term_UCB_seu) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
term_UCB_seu=RunHarmony(term_UCB_seu,group.by.vars='orig.ident')
colnames(term_UCB_seu@meta.data)
colnames(term_UCB_seu@meta.data)=
  
  
load('/home/gibh/2021_NRBC_chlyu/ref_data/NRBC_ref_bullk_RNAseq_se.Rdata',verbose = T)  
term_UCB_seu=singleR_analysis_func(refdata =nrbc_ref_se2,test =GetAssayData(term_UCB_seu,layer = 'data')[VariableFeatures(term_UCB_seu),],outdata = term_UCB_seu,an_type1 ='celltype', an_type2 = 'Pre_celltype')
term_UCB_seu=term_UCB_seu[[1]]
table(is.na(term_UCB_seu$Pre_celltype)) # F:5860
table(term_UCB_seu$Pre_celltype)
term_UCB_seu$Pre_celltype=factor(term_UCB_seu$Pre_celltype,levels =c('BFUE','CFUE','ProE','eBas','lBas','Poly','Orth') )


term_UCB_seu=RunUMAP(term_UCB_seu,dims = 1:8,reduction = 'harmony',return.model = T)
term_UCB_seu=RunUMAP(term_UCB_seu,dims = 1:8,reduction = 'pca',reduction.name = 'pumap')

DimPlot(term_UCB_seu,group.by = 'orig.ident',cols = cols,reduction = 'pumap')/
  DimPlot(term_UCB_seu,group.by = 'celltype',cols = cols,reduction = 'umap')

DimPlot(term_UCB_seu,group.by = 'celltype',cols = cols,reduction = 'umap')
VlnPlot(term_UCB_seu,features = c('KIT','TFRC','GYPA','NCL'),group.by = 'celltype',stack = T)
term_UCB_seu=FindNeighbors(term_UCB_seu,reduction = 'harmony',dims = 1:30)
term_UCB_seu=FindClusters(term_UCB_seu,resolution = 0.2)
DimPlot(term_UCB_seu,group.by = c('RNA_snn_res.0.2','Pre_celltype'),cols = cols)
term_UCB_seu=FindSubCluster(term_UCB_seu,cluster = '4',resolution = 0.25,graph.name = 'RNA_snn')
DimPlot(term_UCB_seu,group.by = c('sub.cluster','Pre_celltype'),cols = cols)

term_UCB_seu$final_celltype=as.character(term_UCB_seu$celltype)
term_UCB_seu$final_celltype[term_UCB_seu$RNA_snn_res.0.2 %in% c('5')]='BFUE/CFUE'
term_UCB_seu$final_celltype[term_UCB_seu$RNA_snn_res.0.2 %in% c('2')]='ProE'
term_UCB_seu$final_celltype[term_UCB_seu$RNA_snn_res.0.2 %in% c('3','4')]='Baso'
term_UCB_seu$final_celltype[term_UCB_seu$RNA_snn_res.0.2 %in% c('1')]='Poly'
term_UCB_seu$final_celltype[term_UCB_seu$RNA_snn_res.0.2 %in% c('0','6')]='Orth'
term_UCB_seu$final_celltype=factor(term_UCB_seu$final_celltype,levels = c('BFUE/CFUE','ProE','Baso','Poly','Orth'))
saveRDS(term_UCB_seu,file = 'NRBC_UCB_altas/ref_term_UCB_NRBC.rds')

DimPlot(term_UCB_seu,group.by = c('final_celltype'),cols = cols)+
  VlnPlot(term_UCB_seu,features = c('KIT','TFRC','GYPA','NCL','CD63'),group.by = 'final_celltype',stack = T)+NoLegend()



other_UCB_NRBC_altas=subset(new_UCB_NRBC_altas,type %in% c('fetal_UCB','pre_UCB'))
other_UCB_NRBC_altas=NormalizeData(other_UCB_NRBC_altas) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
anchors=FindTransferAnchors(reference =term_UCB_seu,query = other_UCB_NRBC_altas,dims = 1:30,reference.reduction = 'pca')
other_UCB_NRBC_altas=MapQuery(anchorset =anchors,query =other_UCB_NRBC_altas  ,reference =term_UCB_seu ,refdata =  list(celltype = "final_celltype"), reference.reduction = "pca", reduction.model = "umap" )
DimPlot(other_UCB_NRBC_altas,reduction = 'ref.umap',group.by = c('predicted.celltype','type'),cols = cols)
other_UCB_NRBC_altas$final_celltype=other_UCB_NRBC_altas$predicted.celltype



new_UCB_NRBC_altas=merge(other_UCB_NRBC_altas,term_UCB_seu)
new_UCB_NRBC_altas[['RNA']]=JoinLayers(new_UCB_NRBC_altas[['RNA']])
new_UCB_NRBC_altas[['prediction.score.celltype']]=NULL
new_UCB_NRBC_altas[['ref.pca']]=NULL
new_UCB_NRBC_altas[['ref.umap']]=NULL
new_UCB_NRBC_altas[['ref.pca']]=merge(other_UCB_NRBC_altas[['ref.pca']],term_UCB_seu[['pca']])
new_UCB_NRBC_altas[['ref.umap']]=merge(other_UCB_NRBC_altas[['ref.umap']],term_UCB_seu[['umap']])
new_UCB_NRBC_altas@meta.data=new_UCB_NRBC_altas@meta.data[,c('orig.ident', 'nCount_RNA', 'nFeature_RNA', 'type', 'stage','ref_dataid',  'percent.mt',
                                                     'percent.rb',  'S.Score', 'G2M.Score', 'Phase' ,'final_celltype')]



mast_cell_genes <- c("HPGDS", "CPA3", "FCER1A")
mast_expr <- GetAssayData(new_UCB_NRBC_altas, assay = "RNA", layer = "data")[mast_cell_genes, ]
# 标记阳性细胞：任一基因表达 > 1
cells_to_remove <- colnames(mast_expr)[which(colSums(mast_expr > 0.2) > 0)];length(cells_to_remove)
new_UCB_NRBC_altas <- subset(new_UCB_NRBC_altas, cells = setdiff(colnames(new_UCB_NRBC_altas), cells_to_remove))
VlnPlot(new_UCB_NRBC_altas,features = c('HPGDS','CPA3','FCER1A' ),stack = T,group.by = 'final_celltype')
saveRDS(new_UCB_NRBC_altas,file = 'NRBC_UCB_altas/dealt_symbol_UCB_NRBC.rds')


p=DimPlot(new_UCB_NRBC_altas,reduction = 'ref.umap',cols = cols,group.by =c('type','final_celltype'),label = F )
ggsave(p,file='NRBC_UCB_altas/res_pic/UCB_nRBC_niche_type_subcelltype_refumap.pdf',width = 12,height = 5)


defintive_markers=readRDS('Protein_NRBC_marker/res_data/main_figure3/defintive_markers.rds')
new_UCB_NRBC_altas$source_celltype=paste(new_UCB_NRBC_altas$type,new_UCB_NRBC_altas$final_celltype)
new_UCB_NRBC_altas$source_celltype=factor(new_UCB_NRBC_altas$source_celltype,levels =c("fetal_UCB BFUE/CFUE","fetal_UCB ProE","fetal_UCB Baso","fetal_UCB Poly", "fetal_UCB Orth" ,
                                                                               "pre_UCB BFUE/CFUE","pre_UCB ProE","pre_UCB Baso","pre_UCB Poly", "pre_UCB Orth" ,
                                                                               "term_UCB BFUE/CFUE","term_UCB ProE","term_UCB Baso","term_UCB Poly", "term_UCB Orth" ) )

VlnPlot(new_UCB_NRBC_altas,group.by = 'source_celltype',features = defintive_markers[['FL_top10_markers']],stack = T,split.by = 'type',cols = alpha(cols[-2],alpha = 0.7))

p=VlnPlot(new_UCB_NRBC_altas,group.by = 'source_celltype',features = c('HBE1','MT1F','GDF15','MINPP1','STOM',defintive_markers[['top_fetal_unique_marker_genes']],defintive_markers[['ABM_cho_topmarkers']]),stack = T,split.by = 'type',cols = alpha(cols[-2],alpha = 0.7))
p
ggsave(p,file='NRBC_UCB_altas/res_pic/UCB_nRBC_niche_subcelltype_niche_marker_expression_vlnplot.pdf',width = 12,height = 6)

library(UCell)

geneset_list=list( 'fetal_signature'=defintive_markers[['top_fetal_unique_marker_genes']],
                   'adult_signature'=defintive_markers[['ABM_cho_topmarkers']])

new_UCB_NRBC_altas <- AddModuleScore_UCell(new_UCB_NRBC_altas, features =geneset_list ,ncores = 6) # 计算速度很快

p=VlnPlot(new_UCB_NRBC_altas,features = c('fetal_signature_UCell','adult_signature_UCell'),group.by = 'source_celltype',stack = T)
ggsave(p,file='NRBC_UCB_altas/res_pic/UCB_fetal_adult_signature_vlntplot.pdf',width = 6.5,height = 6)



filt_NBRC_altas_seu=readRDS('20251125_filt_NBRC_altas_seu.rds')
filt_NBRC_altas_seu=subset(filt_NBRC_altas_seu,tissue_stage!='YS')
filt_NBRC_altas_seu$source_celltype=factor(filt_NBRC_altas_seu$source_celltype,levels = levels(filt_NBRC_altas_seu$source_celltype)[-1:-4])


sample_mexp=AverageExpression(new_UCB_NRBC_altas,features =  rownames(new_UCB_NRBC_altas),group.by = 'source_celltype')$RNA
niche_mexp=AverageExpression(filt_NBRC_altas_seu,features =   rownames(filt_NBRC_altas_seu),group.by = 'source_celltype')$RNA
shared_symbol=rownames(sample_mexp)[rownames(sample_mexp) %in% rownames(niche_mexp)];length(shared_symbol)
all_niche_mexp2=cbind(niche_mexp[shared_symbol,],sample_mexp[shared_symbol,])

cho_gene2=rownames((all_niche_mexp2))[rowSums((all_niche_mexp2)) >1];length(cho_gene2)
p=corrplot(corr = cor(as.matrix(all_niche_mexp2[cho_gene2,])),type = 'lower',method = 'ellipse',diag = F,order = 'hclust',col=alpha(COL2('RdBu',10)[10:1],0.8))

#ggsave(p,file='NRBC_UCB_altas/res_pic/UCB_nRBC_subcelltype_cor_niche_nRBC.pdf',width = 10,height = 10)


sample_mexp=AverageExpression(UCB_NRBC_altas,features =  rownames(UCB_NRBC_altas),group.by = 'type')$RNA
niche_mexp=AverageExpression(filt_NBRC_altas_seu,features =   rownames(filt_NBRC_altas_seu),group.by = 'tissue_stage')$RNA
shared_symbol=rownames(sample_mexp)[rownames(sample_mexp) %in% rownames(niche_mexp)];length(shared_symbol)
all_niche_mexp=cbind(niche_mexp[shared_symbol,],sample_mexp[shared_symbol,])
cho_gene=rownames(all_niche_mexp)[rowSums(all_niche_mexp) >0.5];length(cho_gene)
corrplot.mixed(corr = cor(as.matrix(all_niche_mexp[cho_gene,])),lower = 'number',upper = 'ellipse',order = 'hclust',upper.col=alpha(COL2('RdBu',10)[10:1],0.8))
#ggsave(p,file='NRBC_UCB_altas/res_pic/UCB_nRBC_whole_level_cor_niche_nRBC.pdf',width = 8,height = 6)


saveRDS(new_UCB_NRBC_altas,file = 'NRBC_UCB_altas/dealt_symbol_UCB_NRBC.rds')



