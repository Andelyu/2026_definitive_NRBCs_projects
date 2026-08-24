# 加载所需要的包
# ------------------- require libraries for scRNA-seq-------------------------------#
{
options(stringsAsFactors = F)
.libPaths('/usr/local/lib/R/site-library')
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
set.seed(42)        
library(RColorBrewer)
cols=c(brewer.pal(12,"Set3"),brewer.pal(6,"PiYG"),brewer.pal(6,"BrBG"),brewer.pal(8,"Set2"),brewer.pal(12,"Set3"),brewer.pal(8,"Pastel2"),brewer.pal(9,"Pastel1"),brewer.pal(8,"Accent"))
col=unique(cols)

}

# cal qc func of 10x 

sRNA_qc_cal_func<-function(inpath,proname,seurate_sRNA_raw=NA){
	if(is.na(seurate_sRNA_raw)){
		  sRNA_raw <- Read10X(data.dir  = inpath,gene.column = 2)
      seurate_sRNA_raw <- CreateSeuratObject(counts =sRNA_raw,project = proname,min.cells = 3,min.features = 200)
    }else{Idents(seurate_sRNA_raw)=proname}
    #Idents(seurate_sRNA_raw)<- proname
    #QC and selecting cells for further analysis
    {
	      #The number of unique genes detected in each cell.
	      #Low-quality cells or empty droplets will often have very few genes
	      #Cell doublets or multiplets may exhibit an aberrantly high gene count
	      #Similarly, the total number of molecules detected within a cell (correlates strongly with unique genes)
	      #The percentage of reads that map to the mitochondrial genome
	      #Low-quality / dying cells often exhibit extensive mitochondrial contamination
	      #We calculate mitochondrial QC metrics with the PercentageFeatureSet() function, which calculates the percentage of counts originating from a set of features
	      #We use the set of all genes starting with MT- as a set of mitochondrial genes
	      # #  [[ operator can add columns to object metadata. This is a great place to stash QC stats
     }
      
     seurate_sRNA_raw[['percent.mt']] <- PercentageFeatureSet(seurate_sRNA_raw,pattern = '^MT-' )
     ribo <- c(grep(pattern = "^RPS", x = rownames(seurate_sRNA_raw), value = T),
		 grep(pattern = "^RPL", rownames(seurate_sRNA_raw), value = T))
      seurate_sRNA_raw[['percent.rb']] <- PercentageFeatureSet(seurate_sRNA_raw,features = ribo)
	  
      colnames(seurate_sRNA_raw@meta.data)
      # Visualize QC metrics as a violin plot
      p1=VlnPlot(seurate_sRNA_raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol = 4)
	    
      # FeatureScatter is typically used to visualize feature-feature relationships,but can be used for anything calculated 
      # by the object, i.e. columns in object metadata, PC scores etc.
      p21 = FeatureScatter(seurate_sRNA_raw,feature1 = 'nCount_RNA',feature2 = 'percent.mt')
      p22 = FeatureScatter(seurate_sRNA_raw,feature1 = 'nCount_RNA',feature2 = 'percent.rb')
      p3 = FeatureScatter(seurate_sRNA_raw,feature1 = 'nCount_RNA',feature2 = 'nFeature_RNA')
	      
      # 关于拼图：https://www.jianshu.com/p/73057774b4cb
      p= p1/(p21+p22+p3)+plot_annotation(title = paste(proname,'QC'))
      qc_pare= c(200,as.numeric(quantile(seurate_sRNA_raw$nFeature_RNA,probs = 0.99)),as.numeric(quantile(seurate_sRNA_raw$nCount_RNA,probs = 0.99)),
		  as.numeric(quantile(seurate_sRNA_raw$percent.mt,probs = 0.99)),as.numeric(quantile(seurate_sRNA_raw$percent.rb,probs = 0.99)))
      return(list(seurate_sRNA_raw,p,qc_pare))   
}

# qc and Normalization 
qc_normalization_varfeature_func=function(seu,qc_pare,var_nfeatures=2000,var_method='vst',top_nvar=10){
    qc_seu_sRNA <- subset(seu, subset = nFeature_RNA > qc_pare[1] & nFeature_RNA < qc_pare[2] & nCount_RNA <qc_pare[3] & percent.mt < qc_pare[4] & percent.rb < qc_pare[5])
    # normalize and identify variable feature for each datasets independently #
    nr_qc_seu_sRNA<-NormalizeData(qc_seu_sRNA)
    nr_qc_seu_sRNA=CellCycleScoring(nr_qc_seu_sRNA,s.features =cc.genes$s.genes,g2m.features = cc.genes$g2m.genes,set.ident = T )
    nr_qc_seu_sRNA<-FindVariableFeatures(nr_qc_seu_sRNA,nfeatures =var_nfeatures,selection.method=var_method)
    
    #展示top_nvar 高变异的feature信息
    topn <- head(VariableFeatures(nr_qc_seu_sRNA), top_nvar)
    p1 <- VariableFeaturePlot(nr_qc_seu_sRNA)
    p2 <- LabelPoints(plot = p1, points = topn, repel = TRUE)  
    return(list(nr_qc_seu_sRNA,p2))
}

  
# integrate_multiple_samples ,注意整合时选择的FindVariableFeatures个数最好一致，不能超过max(VariableFeatures)
integrate_sRNAs_func<-function(seu_list,method='LogNormalize',nfeatures=2000){
	# select features that are repeatly vaiable across datasets for integration
	# Prepare an object list normalized with sctransform for integration.
	if(method %in% 'SCT'){ 
	 # seu_list <- lapply(X = seu_list, FUN = SCTransform)
	       for(i in 1:length(seu_list)){
	      	    # 默认排除细胞周期的影响，后续可调试比较不进行细胞周期的变量排除	
              seu_list[[i]] <- SCTransform(seu_list[[i]], verbose = F, vars.to.regress = c('nCount_RNA','percent.mt','S.Score'))
         }
          features<- SelectIntegrationFeatures(object.list = seu_list,nfeatures =nfeatures)
          seu_list <- PrepSCTIntegration(object.list = seu_list,anchor.features = features,verbose = F)
    }else{features<- SelectIntegrationFeatures(object.list = seu_list,nfeatures =nfeatures )  }
  
    #perform integration 
    # find feature anchors according to features from selected features
    find_anchors <- FindIntegrationAnchors(seu_list,anchor.features = features,normalization.method  = method)
    integrate_sRNA <- IntegrateData(anchorset =find_anchors,new.assay.name = 'integrated')
    DefaultAssay(integrate_sRNA) <- "integrated"
    return(integrate_sRNA)
}


# 预留参数接口，目前暂只是用PCA线性降维的方法，
visual_umap_func<-function(sRNA_count,npca=30,resolution=0.4,reduction1='pca',reduction2='umap',Integ=F,n.neighbors=30L,vars.to.regress=NULL){
        # 计算细胞周期
       sRNA_count=CellCycleScoring(sRNA_count,s.features =cc.genes$s.genes,g2m.features = cc.genes$g2m.genes,set.ident = T )
       sRNA_count <-ScaleData(sRNA_count,vars.to.regress =vars.to.regress )
       sRNA_count<- RunPCA(sRNA_count,npcs=50)
       
      
      sRNA_count <- FindNeighbors(sRNA_count,dims=1:npca)
      sRNA_count <- FindClusters(sRNA_count,resolution=resolution,method = 'igraph')
      p0=ElbowPlot(sRNA_count,ndims = 30)+ DimPlot(sRNA_count,reduction = 'pca',label = T,repel = T)
        
      sRNA_count <- RunUMAP(sRNA_count,reduction=reduction1,dims=1:npca,n.neighbors =n.neighbors)
      sRNA_count <- RunTSNE(sRNA_count,reduction=reduction1,dims=1:npca,check_duplicates=F)

        
      # visualization 
      cols=c(brewer.pal(12,"Set3"),brewer.pal(6,"PiYG"),brewer.pal(6,"BrBG"),brewer.pal(8,"Set2"),brewer.pal(12,"Set3"),brewer.pal(8,"Pastel2"),brewer.pal(9,"Pastel1"),brewer.pal(8,"Accent"))
      col=unique(cols)
      
      p1 <- DimPlot(sRNA_count,reduction = 'umap',label = T,repel = T,cols = col)+DimPlot(sRNA_count,reduction = 'tsne',label = T,repel = T,cols = col)
	if(Integ){
		      p2 <- DimPlot(sRNA_count,reduction = reduction2,group.by = 'orig.ident',cols = col)
	        p3 <- DimPlot(sRNA_count,reduction = reduction2,split.by= 'orig.ident',label = T,repel = T,cols = col)
	        p1 <- p1/p2/p3
	 }
	  
	 p4=FeaturePlot(sRNA_count,features = c('nCount_RNA','nFeature_RNA','percent.mt','percent.rb'),label = T)
	    # 细胞周期结果
	 p5=DimPlot(sRNA_count,group.by = 'Phase')
	    
	 p1<- p1+plot_annotation(paste0('ParM: ','npca=',npca,' resolution=',resolution))
	 p2=p4/p5
	 return(list(sRNA_count,p0,p1,p2))
}

# 单独测试 聚类结果
cluster_pca_func<-function(sRNA_count,npca=30,resolution=c(0.3,0.4,0.5,0.6,0.8,1)){
  sRNA_count <- FindNeighbors(sRNA_count,dims=1:npca)
  sRNA_count <- FindClusters(sRNA_count,resolution=resolution,method = 'igraph')
  if(grep('SCT_snn_res',colnames(sRNA_count@meta.data))>0){
    p=clustree(sRNA_count@meta.data,prefix = 'SCT_snn_res.',node_label_aggr = "median")
  }else if(grep('integrated_snn_res',colnames(sRNA_count@meta.data))>0){
    p=clustree(sRNA_count@meta.data,prefix = 'integrated_snn_res',node_label_aggr = "median")
  }else{
    p=clustree(sRNA_count@meta.data,prefix = 'snn_res',node_label_aggr = "median")
  }
  return(p)
}

# 抽取亚类再进行分析聚类分析




# 载入参考数据
# 在容器id：8414c426f27c初始产生这些数据的脚本,存放记录而已
if(F){
  fetal_kiedny_mexp=read.csv('/home/data/fetal_celltype_nr_mexp.csv',header = T)
  rownames(fetal_kiedny_mexp)=fetal_kiedny_mexp$X;fetal_kiedny_mexp=fetal_kiedny_mexp[,-1]
  colData=data.frame('celltype'=colnames(fetal_kiedny_mexp))
  fetal_kiedney_se <- SummarizedExperiment(assays=list(counts=fetal_kiedny_mexp),colData=colData)
  saveRDS(fetal_kiedney_se,file = '~/fetal_kiedney_mexp_se.rds')
  
  library(celldex)
  hpca.se <- HumanPrimaryCellAtlasData()
  saveRDS(hpca.se,file = '~/HumanPrimaryCellAtlasData.rds')
  
}

singleR_analysis_func<-function(test,refdata,outdata,an_type1='celltype',an_type2='anno_type',test_assay=1,ref_assay=1,clusters=NULL){
     library(SingleR );library(ggplotify)
    #for single-cell references
	  #testdata: A numeric matrix of (usually log-transformed) expression values from a reference dataset, 
	  #refdata: A numeric matrix or SummarizedExperiment  对象
	  #outdata : 被注释seu类型数据
	  #an_type1 ：refdata中的标签colname；an_type2： 被注释标签colname
    #适用于的单细胞数据 
	  #if(!is.matrix(test)){test=as.matrix(test)}
          pred_singler_res <- SingleR(test =test , ref = refdata, assay.type.test=test_assay,labels = refdata@colData[,an_type1],assay.type.ref =ref_assay ,clusters = clusters)
          p1=as.ggplot(plotScoreHeatmap(pred_singler_res))
          outdata@meta.data[,an_type2] <- pred_singler_res[rownames(outdata@meta.data),'pruned.labels'] 
          p2=DimPlot(outdata,label = T)/DimPlot(outdata,group.by = an_type2,label = T,cols = col) 
          #print(p2)
          #print(data.frame(table( outdata@meta.data[,an_type2])))
         return(list(outdata,p1,p2))
}





# --------------Constructing single-cell trajectories------------ #
#   monocle3 安装
if(!require(monocle3)){
  library(BiocManager);
  if(!require(monocle)){BiocManager::install('monocle')} # leidenbase依赖此包
  BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                         'limma', 'S4Vectors', 'SingleCellExperiment',
                         'SummarizedExperiment', 'batchelor', 'Matrix.utils'))
  if(!require(devtools)){install.packages("devtools")}
  devtools::install_github('cole-trapnell-lab/leidenbase')
  devtools::install_github('cole-trapnell-lab/monocle3')
}else{library(monocle3)}

monocle3_trajectories_func<-function(read_cho=T,seu=NULL,seura_path,mono_date_cds=NULL,preprocess_method='PCA',reduction_method='UMAP',pc_num_dim=30,
                                     resolution=1e-4,color_type='celltype',use_genes =NULL,align=F,alignment_group=NULL,res_str=NULL,plotmain='',assay='RNA'){
  # Step 1: load  10X dataset or use the mono_date_cds and pre_porcess the cds
  if(read_cho){ 
    mono_date_cds <- load_cellranger_data(seura_path)
  }else if(!is.null(seu)){
    data <- GetAssayData(seu, assay = assay, slot = 'counts')
    cell_metadata <- seu@meta.data
    gene_annotation <- data.frame(gene_short_name = rownames(data))
    rownames(gene_annotation) <- rownames(data)
    mono_date_cds <- new_cell_data_set(expression_data = as.matrix(data),cell_metadata = cell_metadata,gene_metadata = gene_annotation)
    rm(data)
   # if(is.null(seu@assays$integrated@var.features) & is.null(use_genes)){use_genes=seu@assays$integrated@var.features}
    if( is.null(use_genes)){
      Idents(seu)=color_type
      use_genes=FindAllMarkers(seu) %>% group_by(cluster) %>% filter(p_val_adj <0.05)  %>% top_n(n =20,wt = avg_log2FC  )
    }
  }
  #mono_date_cds= estimate_size_factors(mono_date_cds)
  # Step 2:normalize 
  mono_date_cds <- preprocess_cds(mono_date_cds,method = preprocess_method,num_dim = pc_num_dim,use_genes  =use_genes )
  p1=plot_pc_variance_explained(mono_date_cds)
  #去除批次效应和一些变量的影响
  if(align){
    print('---------------align-------------')
    if(is.null(res_str)){res_str=c('~percent.mt+nFeature_RNA+S.Score+G2M.Score')}
    mono_date_cds=align_cds(cds =mono_date_cds,alignment_group = alignment_group ,residual_model_formula_str = res_str)
  }
  # Step 3:reduce_dimension
  mono_date_cds <- reduce_dimension(mono_date_cds, preprocess_method=preprocess_method,reduction_method=reduction_method,max_components = 10)
  
  # Step 4:cluster
  mono_date_cds <- cluster_cells(mono_date_cds,resolution =resolution)
  p2=plot_cells(mono_date_cds, show_trajectory_graph = F,color_cells_by =color_type,label_groups_by_cluster = F)
  p3=plot_cells(x = 1,y = 3,mono_date_cds, show_trajectory_graph = F)
  
  # plot_cells(mono_date_cds, show_trajectory_graph = F,color_cells_by =  'celltype')
  ## Step 5: Learn a graph
  mono_date_cds <- learn_graph(mono_date_cds)
  #plot_cells(mono_date_cds)
  p4=plot_cells(mono_date_cds,show_trajectory_graph = T,group_label_size = 5,cell_size = 0.6,color_cells_by = color_type,label_groups_by_cluster = F)
  
  ## Step 6: Order cells
  mono_date_cds <- order_cells(mono_date_cds)
  
  p=p1+p2+p3+p4 + plot_annotation(title = paste(plotmain,preprocess_method,reduction_method,pc_num_dim,resolution,sep = ": "))
  return(list(mono_date_cds,p))
}











