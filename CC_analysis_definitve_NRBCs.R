
library(Seurat)
library(CellChat)
library(nichenetr)
library(clusterProfiler)
library(reshape)
library(ggplot2)
library(patchwork )
library(pheatmap)
library(network)
library(ggnetwork)

library(RColorBrewer )
cols=c(brewer.pal(12,"Set3"),brewer.pal(6,"PiYG"),brewer.pal(6,"BrBG"),brewer.pal(8,"Set2"),
       brewer.pal(12,"Set3"),brewer.pal(8,"Pastel2"),brewer.pal(9,"Pastel1"),brewer.pal(8,"Accent"))

signaling_pic_func=function(ligands_oi,targets_oi,weighted_networks=weighted_networks,ligand_tf_matrix=ligand_tf_matrix,top_n_regulators=4){
  
  active_signaling_network <- get_ligand_signaling_path(ligands_all = ligands_oi,targets_all = targets_oi,
                                                        weighted_networks = weighted_networks,ligand_tf_matrix = ligand_tf_matrix,
                                                        top_n_regulators = top_n_regulators,minmax_scaling = TRUE) 
  
  graph_min_max <- diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network,
                                                     ligands_all = ligands_oi, targets_all = targets_oi,
                                                     sig_color = "indianred", gr_color = "steelblue")
  
  
  #构建 network 对象
  # 使用 edges_df 中的 from 和 to 列构建有向图
  edges_df <- graph_min_max$edges_df
  nodes_df <- graph_min_max$nodes_df
  net_obj <- network(as.matrix(edges_df[, c("from", "to")]), directed = TRUE)
  
  # 将节点属性合并到 network 对象中
  # 这一步至关重要，它让 ggnetwork 能读取节点的颜色、类型等信息
  # 确保节点名称匹配
  net_vertices <- network.vertex.names(net_obj)
  
  # 将属性赋值给 network 对象的顶点
  # 使用 match 确保顺序一致
  net_obj %v% "label" <- nodes_df$label[match(net_vertices, nodes_df$id)]
  net_obj %v% "type" <- nodes_df$type[match(net_vertices, nodes_df$id)] # 节点类型（配体/靶点/中间分子）
  net_obj %v% "fillcolor" <- nodes_df$fillcolor[match(net_vertices, nodes_df$id)] # 节点颜色
  
  # 为了映射边的粗细，我们需要将边的权重信息也加入到 network 对象中
  # 注意：network 对象构建时边的顺序可能与 edges_df 不完全一致，需要重新匹配
  # 获取 network 对象中的边列表
  net_edges <- as.data.frame(net_obj, edge.names = TRUE)
  # 匹配 from 和 to 来合并权重
  net_edges$weight <- edges_df$penwidth[match(paste(net_edges$.tail, net_edges$.head, sep = "-"), 
                                              paste(edges_df$from, edges_df$to, sep = "-"))]
  # 将权重赋值给 network 对象的边属性
  net_obj %e% "weight" <- net_edges$weight
  
  # 开始绘图
  p <- ggplot(net_obj, aes(x = x, y = y, xend = xend, yend = yend))+
    # 绘制边 使用 weight 映射线条粗细
    geom_edges(aes(size = weight),alpha = 0.4,  color = "grey70", arrow = arrow(length = unit(1, "mm"), type = "closed")) +
    # 绘制节点 使用 fillcolor 映射填充颜色
    geom_nodes(aes(color = fillcolor), size = 5) +
    # 添加节点标签
    geom_nodetext(aes(label = label),   fontface = "bold", color = "black", size = 3, vjust = 1.5) + # vjust 调整标签位置
    # 主题调整, 移除背景网格
    theme_blank() + labs(title = "NicheNet Signaling Pathway") +
    # 手动设置颜色（保持 NicheNet 原有的配色）
    scale_color_identity() + scale_size_continuous(range = c(0.5, 2)) # 调整线条粗细范围
  
  print(unique(net_obj$label))
  return(p)
}


NRBC_subcelltype=c("BFUE/CFUE","ProE","Bas","Poly","Orth" )

NRBC_altas_LR_df=read.csv('NRBC_altas_CC/res_data/filt_NRBC_altas_LR_df.csv',sep="\t")
NRBC_altas_LR_df=NRBC_altas_LR_df[NRBC_altas_LR_df$stage!='YS',]


all_NRBC_receptor_genes=CellChatDB.human$interaction[CellChatDB.human$interaction$interaction_name %in% unique(NRBC_altas_LR_df$interaction_name[NRBC_altas_LR_df$target_type!='Ery2Other']),'receptor.symbol'] # 取receptor gene
all_NRBC_ligand_genes=CellChatDB.human$interaction[CellChatDB.human$interaction$interaction_name %in%unique(NRBC_altas_LR_df$interaction_name[NRBC_altas_LR_df$target_type=='Ery2Other']),'ligand.symbol']# 取ligand gene
all_NRBC_receptor_genes=sort(unique(unlist(strsplit(all_NRBC_receptor_genes,split=','))))
all_NRBC_ligand_genes=sort(unique(unlist(strsplit(all_NRBC_ligand_genes,split=','))))
length(all_NRBC_receptor_genes);length(all_NRBC_ligand_genes)# 49， 77



filt_NBRC_altas_seu=readRDS('20251125_filt_NBRC_altas_seu.rds')
NRBC_altas_LR_df$celltype=factor(NRBC_altas_LR_df$celltype, levels = levels(filt_NBRC_altas_seu$final_celltype))

############################################################################################################################################
#-----------the specific interactions ---------------------------#
############################################################################################################################################
fetal_adult_NRBC_whole_marker=readRDS('Protein_NRBC_marker/res_data/HSPC_derived_nRBC_wholelevel_RNA_markers.rds')
sub_fetal_adult_all_Ery_tissue_markers=read.csv('Protein_NRBC_marker/DE_marker/fetal_adult_all_Ery_RNA_markers.csv')
sub_fetal_adult_all_Ery_tissue_markers=sub_fetal_adult_all_Ery_tissue_markers[,-1]

Ery_recepted_LR_df=NRBC_altas_LR_df[NRBC_altas_LR_df$target_type!='Ery2Other',];
Ery_outgoing_LR_df=NRBC_altas_LR_df[NRBC_altas_LR_df$target_type=='Ery2Other',]

# recepted LR interaction
fetal_nRBC_recepted_lr_list=unique(Ery_recepted_LR_df[Ery_recepted_LR_df$stage %in% c('FL','FBM'),'interaction_name']);length(fetal_nRBC_recepted_lr_list)
adult_nRBC_recepted_lr_list=unique(Ery_recepted_LR_df[Ery_recepted_LR_df$stage %in% c('ABM'),'interaction_name']);length(adult_nRBC_recepted_lr_list)

specific_adult_nRBC_recepted_lr_list=adult_nRBC_recepted_lr_list[!adult_nRBC_recepted_lr_list %in% fetal_nRBC_recepted_lr_list];length(specific_adult_nRBC_recepted_lr_list)
specific_fetal_nRBC_recepted_lr_list=fetal_nRBC_recepted_lr_list[!fetal_nRBC_recepted_lr_list %in% adult_nRBC_recepted_lr_list];length(specific_fetal_nRBC_recepted_lr_list)

cho_Ery_recepted_LR_df=Ery_recepted_LR_df[Ery_recepted_LR_df$interaction_name %in% c(specific_fetal_nRBC_recepted_lr_list,specific_adult_nRBC_recepted_lr_list),]
cho_Ery_recepted_LR_df$stage=factor(cho_Ery_recepted_LR_df$stage,levels = c('FL','FBM','ABM'))
cho_Ery_recepted_LR_df$type='fetal';cho_Ery_recepted_LR_df$type[cho_Ery_recepted_LR_df$interaction_name %in% specific_adult_nRBC_recepted_lr_list]='adult';
cho_Ery_recepted_LR_df=cho_Ery_recepted_LR_df[order(cho_Ery_recepted_LR_df$type,cho_Ery_recepted_LR_df$annotation,cho_Ery_recepted_LR_df$receptor),]
cho_Ery_recepted_LR_df$interaction_name=factor(cho_Ery_recepted_LR_df$interaction_name,levels = unique(cho_Ery_recepted_LR_df$interaction_name))

LR_max_prob_func=function(LR_df){
  temp_max_prob_LRs_df=data.frame() # interaction,stage,celltype,target_type
  for(LR in unique(LR_df$interaction_name) ){
    temp_df=LR_df[LR_df$interaction_name ==LR ,]
    temp_df$stage_celltype=paste(as.character(temp_df$stage),as.character(temp_df$celltype),sep='_')
    if(dim(temp_df)[1]>1){
      for (celltype  in unique(temp_df$stage_celltype)) {
        temp_df1=temp_df[temp_df$stage_celltype==celltype,]
        if(length(unique(temp_df1$target_type)) >1){
          for( target_type in unique(temp_df1$target_type) ){
            temp_df2=temp_df1[temp_df1$target_type==target_type,]
            temp_df2=temp_df2[which.max(temp_df2$prob),]
            temp_max_prob_LRs_df=rbind(temp_max_prob_LRs_df,temp_df2)
          }
        }else{
          temp_df2=temp_df1[which.max(temp_df1$prob),]
          temp_max_prob_LRs_df=rbind(temp_max_prob_LRs_df,temp_df2)
          
        }
      }
    }else{
      temp_df2=temp_df[which.max(temp_df$prob),]
      temp_max_prob_LRs_df=rbind(temp_max_prob_LRs_df,temp_df2)
      
    }
    
  }
  colnames(temp_max_prob_LRs_df)=gsub(pattern ='prob' ,replacement = 'max.prob',colnames(temp_max_prob_LRs_df))
  return(temp_max_prob_LRs_df)
  
}

max_cho_Ery_recepted_LR_df=LR_max_prob_func(LR_df =cho_Ery_recepted_LR_df )

p1=ggplot(max_cho_Ery_recepted_LR_df,aes(x=celltype,y=interaction_name,size=max.prob,shape=target_type,color=annotation))+geom_point()+theme_classic()+scale_shape_manual(values = c(0:2,5))+
  scale_fill_manual(values = cols[1:4])+theme(axis.text.x = element_text(angle = 45,hjust = 1),text = element_text(face ='bold'))+ggtitle(label = 'NRBCs recepted specific interactions')+
  facet_grid(~stage )+scale_color_manual(values = cols[-2])
p1

setwd('definitive_nRBC_marker_project')
dir.create('res_pic/main_figure3')
ggsave(p1,file='res_pic/main_figure3/stage_specific_recepted_interaction.pdf',width =8 ,height = 12)


# outgoing LR interaction
fetal_nRBC_outgoing_lr_list=unique(Ery_outgoing_LR_df[Ery_outgoing_LR_df$stage %in% c('FL','FBM'),'interaction_name']);length(fetal_nRBC_outgoing_lr_list)
adult_nRBC_outgoing_lr_list=unique(Ery_outgoing_LR_df[Ery_outgoing_LR_df$stage %in% c('ABM'),'interaction_name']);length(adult_nRBC_outgoing_lr_list)

specific_adult_nRBC_outgoing_lr_list=adult_nRBC_outgoing_lr_list[!adult_nRBC_outgoing_lr_list %in% fetal_nRBC_outgoing_lr_list];length(specific_adult_nRBC_outgoing_lr_list)
specific_fetal_nRBC_outgoing_lr_list=fetal_nRBC_outgoing_lr_list[!fetal_nRBC_outgoing_lr_list %in% adult_nRBC_outgoing_lr_list];length(specific_fetal_nRBC_outgoing_lr_list)

cho_Ery_outgoing_LR_df=Ery_outgoing_LR_df[Ery_outgoing_LR_df$interaction_name %in% c(specific_fetal_nRBC_outgoing_lr_list,specific_adult_nRBC_outgoing_lr_list),]
cho_Ery_outgoing_LR_df$stage=factor(cho_Ery_outgoing_LR_df$stage,levels = c('FL','FBM','ABM'))
cho_Ery_outgoing_LR_df$type='fetal';cho_Ery_outgoing_LR_df$type[cho_Ery_outgoing_LR_df$interaction_name %in% specific_adult_nRBC_outgoing_lr_list]='adult';
cho_Ery_outgoing_LR_df=cho_Ery_outgoing_LR_df[order(cho_Ery_outgoing_LR_df$type,cho_Ery_outgoing_LR_df$annotation,cho_Ery_outgoing_LR_df$receptor),]
cho_Ery_outgoing_LR_df$interaction_name=factor(cho_Ery_outgoing_LR_df$interaction_name,levels = unique(cho_Ery_outgoing_LR_df$interaction_name))

max_cho_Ery_outgoing_LR_df=LR_max_prob_func(LR_df =cho_Ery_outgoing_LR_df )
p2=ggplot(max_cho_Ery_outgoing_LR_df,aes(x=celltype,y=interaction_name,size=max.prob,color=annotation))+geom_point()+theme_classic()+scale_shape_manual(values = c(0:2,5))+
  scale_fill_manual(values = cols[1:4])+theme(axis.text.x = element_text(angle = 45,hjust = 1),text = element_text(face ='bold'))+ggtitle(label = 'NRBCs recepted specific interactions')+
  facet_grid(~stage )+scale_color_manual(values = cols[-2])
p2

ggsave(p1,file='res_pic/main_figure3/stage_specific_outgoing_interaction.pdf',width =8 ,height = 16)





############################################################################################################################################################
#-----------------------------using the nichnet to analysis the key LR to target genes-------------------------#
############################################################################################################################################################
library(nichenetr)
library(ggplot2)
library(dplyr)

lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
ligand_target_matrix <- readRDS("NRBC_altas_CC/ligand_target_matrix_nsga2r_final.rds")
weighted_networks <- readRDS("NRBC_altas_CC/weighted_networks_nsga2r_final.rds")
signaling_pic_func=function(ligands_oi,targets_oi,weighted_networks=weighted_networks,ligand_tf_matrix=ligand_tf_matrix,top_n_regulators=4){
  
  active_signaling_network <- get_ligand_signaling_path(ligands_all = ligands_oi,targets_all = targets_oi,
                                                        weighted_networks = weighted_networks,ligand_tf_matrix = ligand_tf_matrix,
                                                        top_n_regulators = top_n_regulators,minmax_scaling = TRUE) 
  
  graph_min_max <- diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network,
                                                     ligands_all = ligands_oi, targets_all = targets_oi,
                                                     sig_color = "indianred", gr_color = "steelblue")
  
  
  #构建 network 对象
  # 使用 edges_df 中的 from 和 to 列构建有向图
  edges_df <- graph_min_max$edges_df
  nodes_df <- graph_min_max$nodes_df
  net_obj <- network(as.matrix(edges_df[, c("from", "to")]), directed = TRUE)
  
  # 将节点属性合并到 network 对象中
  # 这一步至关重要，它让 ggnetwork 能读取节点的颜色、类型等信息
  # 确保节点名称匹配
  net_vertices <- network.vertex.names(net_obj)
  
  # 将属性赋值给 network 对象的顶点
  # 使用 match 确保顺序一致
  net_obj %v% "label" <- nodes_df$label[match(net_vertices, nodes_df$id)]
  net_obj %v% "type" <- nodes_df$type[match(net_vertices, nodes_df$id)] # 节点类型（配体/靶点/中间分子）
  net_obj %v% "fillcolor" <- nodes_df$fillcolor[match(net_vertices, nodes_df$id)] # 节点颜色
  
  # 为了映射边的粗细，我们需要将边的权重信息也加入到 network 对象中
  # 注意：network 对象构建时边的顺序可能与 edges_df 不完全一致，需要重新匹配
  # 获取 network 对象中的边列表
  net_edges <- as.data.frame(net_obj, edge.names = TRUE)
  # 匹配 from 和 to 来合并权重
  net_edges$weight <- edges_df$penwidth[match(paste(net_edges$.tail, net_edges$.head, sep = "-"), 
                                              paste(edges_df$from, edges_df$to, sep = "-"))]
  # 将权重赋值给 network 对象的边属性
  net_obj %e% "weight" <- net_edges$weight
  
  # 开始绘图
  p <- ggplot(net_obj, aes(x = x, y = y, xend = xend, yend = yend))+
    # 绘制边 使用 weight 映射线条粗细
    geom_edges(aes(size = weight),alpha = 0.4,  color = "grey70", arrow = arrow(length = unit(1, "mm"), type = "closed")) +
    # 绘制节点 使用 fillcolor 映射填充颜色
    geom_nodes(aes(color = fillcolor), size = 5) +
    # 添加节点标签
    geom_nodetext(aes(label = label),   fontface = "bold", color = "black", size = 3, vjust = 1.5) + # vjust 调整标签位置
    # 主题调整, 移除背景网格
    theme_blank() + labs(title = "NicheNet Signaling Pathway") +
    # 手动设置颜色（保持 NicheNet 原有的配色）
    scale_color_identity() + scale_size_continuous(range = c(0.5, 2)) # 调整线条粗细范围
  
  print(unique(net_obj$label))
  return(p)
}

draw_gonetcwork_pic_func=function(res=fa_fetal_LR_targetgene_enrichGO_res,showCategory=20,xlimits = c(-0.2, 2.5)){
  
  res <- as.data.frame(res)
  res <- head(res,showCategory)
  
  # 构建连接
  links <- do.call(rbind, lapply(1:nrow(res), function(i) {
    data.frame(gene = strsplit(res$geneID[i], "/")[[1]], 
               pathway = res$Description[i])}))
  
  # 计算位置
  links$y_gene <- match(links$gene, names(sort(table(links$gene), T))) / (length(unique(links$gene)) + 1)
  links$y_pathway <- match(links$pathway, res$Description[order(res$p.adjust)]) / (nrow(res) + 1)
  
  # 绘图
  p=ggplot() +  geom_segment(aes(x = 0, xend = 0.8, y = y_gene, yend = y_pathway), data = links,   color = "gray70", alpha = 0.4, size = 0.4) +
    geom_point(aes(x = 0, y = y_gene), data = distinct(links, gene, y_gene), size = 2, color = "#4878D0") +
    geom_text(aes(x = -0.05, y = y_gene, label = gene), data = distinct(links, gene, y_gene), hjust = 1, size = 3, color = "black") +
    geom_point(aes(x = 0.8, y = y_pathway), data = distinct(links, pathway, y_pathway), size = 2, color = "#EE854A") +
    geom_text(aes(x = 0.85, y = y_pathway, label = pathway), data = distinct(links, pathway, y_pathway),   hjust = 0, size = 2.5, color = "black") +
    scale_x_continuous(limits =xlimits, breaks = c(0, 1), labels = c("Genes", "Pathways")) +
    scale_y_continuous(breaks = NULL) +
    theme_classic() +
    theme(axis.text.x = element_text(face = "bold", size = 10),
          plot.margin = margin(5, 15, 25, 5))
  print(p)
  return(p)
}


filt_NBRC_altas_seu$fa_celltype2=paste(as.character(filt_NBRC_altas_seu$fa_type),as.character(filt_NBRC_altas_seu$final_celltype),sep = "_")
fa_all_mexp_df=as.matrix(AverageExpression(subset(filt_NBRC_altas_seu,tissue_stage !='YS'),group.by = 'fa_celltype2',features =rownames(filt_NBRC_altas_seu) )$RNA)

known_key_regulator_erythropoiesis_genes=c('GATA1','GATA2','KLF1','TAL1','LMO2','FOG1','BCL11A','EPOR','HIF1A','HIF2A','STAT5A','ALAS2','SLC11A2','FOXO3')



# --------compare the LR on primitive and definitve NRBC -------------#

LR_df=CellChatDB.human$interaction

fetal_NBRC_ligands=unique(NRBC_altas_LR_df[NRBC_altas_LR_df$stage %in% c('FL','FBM') & NRBC_altas_LR_df$target_type!='Ery2Other','interaction_name'])
fetal_NBRC_ligands=unique(LR_df$ligand.symbol[LR_df$interaction_name %in%fetal_NBRC_ligands ])
fetal_NBRC_ligands= unique(as.character(t(data.frame(strsplit(fetal_NBRC_ligands,split = ', ')))[,1]))
fetal_NBRC_ligands[!fetal_NBRC_ligands %in% rownames(filt_NBRC_altas_seu)]  

adult_NBRC_ligands=unique(NRBC_altas_LR_df[NRBC_altas_LR_df$stage %in% c('ABM') & NRBC_altas_LR_df$target_type!='Ery2Other','interaction_name'])
adult_NBRC_ligands=unique(LR_df$ligand.symbol[LR_df$interaction_name %in%adult_NBRC_ligands ])
adult_NBRC_ligands= unique(as.character(t(data.frame(strsplit(adult_NBRC_ligands,split = ', ')))[,1]))
adult_NBRC_ligands[!adult_NBRC_ligands %in% rownames(filt_NBRC_altas_seu)]  
adult_NBRC_ligands=c(adult_NBRC_ligands,'EPO')

# --------- HSPC-derived : fetal vs adult ------------------#
top_hspc_degs=fetal_adult_NRBC_whole_marker[fetal_adult_NRBC_whole_marker$avg_log2FC >1 & fetal_adult_NRBC_whole_marker$pct.1 >0.05 & fetal_adult_NRBC_whole_marker$pct.2 <0.3, ] %>% group_by(cluster) %>% top_n(wt =avg_log2FC,n = 500 )
top_hspc_sub_degs=sub_fetal_adult_all_Ery_tissue_markers[sub_fetal_adult_all_Ery_tissue_markers$avg_log2FC >1 & sub_fetal_adult_all_Ery_tissue_markers$pct.1>0.1,]  %>% group_by(cluster,celltype) %>% do(head(.,n=300))
table(top_hspc_sub_degs$cluster);table(top_hspc_degs$cluster);

fetal_target_genes=unique(c(top_hspc_degs$gene[top_hspc_degs$cluster=='fetal'],top_hspc_sub_degs$gene[top_hspc_sub_degs$cluster=='fetal']))
adult_target_genes=unique(c(top_hspc_degs$gene[top_hspc_degs$cluster=='adult'],top_hspc_sub_degs$gene[top_hspc_sub_degs$cluster=='adult']))
length(fetal_target_genes);length(adult_target_genes)


target_genes=fetal_target_genes;
target_genes=target_genes[rowMax(fa_all_mexp_df[target_genes,grep('fetal',colnames(fa_all_mexp_df))])>0.5]
length(target_genes)
potential_ligands=fetal_NBRC_ligands[fetal_NBRC_ligands %in% colnames(ligand_target_matrix)]

fetal_ligand_target_res_list=nichenet_predict_func(legend_title = 'fetal NRBC ligand to target gene score',potential_ligands =potential_ligands,target_genes = target_genes,ligand_target_matrix = ligand_target_matrix,background_expressed_genes = background_expressed_genes )
fetal_ligand_target_res_list[[2]]
ggsave(fetal_ligand_target_res_list[[2]],filename = 'res_pic/main_figure4/fetal_recept_ligand_target_gene_heatmap_score.pdf',width = 10,height = 10)

saveRDS(fetal_ligand_target_res_list,file = 'res_data/fa_fetal_ligand_target_res_list.rds')

#
an_df=data.frame(fetal_ligand_target_res_list[[1]]$data[fetal_ligand_target_res_list[[1]]$data$y %in% c('EPO','IGF1','IFNG','IGF2','TNF','LTB','LTA'),])
an_df=data.frame(row.names = an_df$y,aupr_score=an_df$score)

p=pheatmap(t(fetal_ligand_target_res_list[[3]][,c('EPO','IGF1','IFNG','IGF2','TNF','LTB','LTA')]),annotation_row = an_df,border_color = 'white',cluster_rows = F,cluster_cols =F ,color =colorRampPalette(colors =  c("white" ,"#EF6548"))(100))
ggsave(as.ggplot(p),file='res_pic/main_figure3/fetal_key_ligand_target_genes_heatmap.pdf',width = 18,height = 3,dpi = 300)

fa_fetal_LR_targetgene_enrichGO_res=enrichGO(gene =rownames(fetal_ligand_target_res_list[[3]]),OrgDb = org.Hs.eg.db,keyType = 'SYMBOL',ont = 'BP' )
cnetplot(fa_fetal_LR_targetgene_enrichGO_res,showCategory =20,   color_category = "#E41A1C", color_gene = "#377EB8", cex_label_category = 0.6,cex_label_gene = 0.5, node_label = "all",  layout = "fr")
p=draw_gonetcwork_pic_func(res =fa_fetal_LR_targetgene_enrichGO_res,showCategory = 20 )
ggsave(p,filename='res_pic/main_figure3/fetal_target_gene_enrichGO_res.pdf',width =6,height = 6)
saveRDS(fa_fetal_LR_targetgene_enrichGO_res,file = 'res_data/fa_fetal_LR_targetgene_enrichGO_res.rds')

cho_fetal_ligand_target_df=fetal_ligand_target_res_list[[3]][,c('EPO','IGF1','IFNG','IGF2','TNF','LTB','LTA')]
cho_fetal_ligand_target_list=list()
cho_fetal_ligand_target_list[['EPO']]=names(cho_fetal_ligand_target_df[,'IGF1'])[cho_fetal_ligand_target_df[,'EPO'] >0.04]
cho_fetal_ligand_target_list[['IGF1']]=names(cho_fetal_ligand_target_df[,'IGF1'])[cho_fetal_ligand_target_df[,'IGF1'] >0.04]
cho_fetal_ligand_target_list[['IFNG']]=names(cho_fetal_ligand_target_df[,'IGF1'])[cho_fetal_ligand_target_df[,'IFNG'] >0.04]
cho_fetal_ligand_target_list[['TNF']]=names(cho_fetal_ligand_target_df[,'TNF'])[cho_fetal_ligand_target_df[,'TNF'] >0.04]
cho_fetal_ligand_target_list[['LTB']]=names(cho_fetal_ligand_target_df[,'EPO'])[cho_fetal_ligand_target_df[,'LTB'] >0.04]
cho_fetal_ligand_target_list[['LTB']]=names(cho_fetal_ligand_target_df[,'LTA'])[cho_fetal_ligand_target_df[,'LTA'] >0.04]


#----------------------adult ------------------#
target_genes=adult_target_genes;
target_genes=target_genes[rowMax(fa_all_mexp_df[target_genes,grep('adult',colnames(fa_all_mexp_df))])>0.5]
length(target_genes)
potential_ligands=adult_NBRC_ligands[adult_NBRC_ligands %in% colnames(ligand_target_matrix)]
adult_ligand_target_res_list=nichenet_predict_func(legend_title = 'adult NRBC ligand to target gene score',potential_ligands =potential_ligands,target_genes = target_genes,ligand_target_matrix = ligand_target_matrix,background_expressed_genes = background_expressed_genes )
p=adult_ligand_target_res_list[[2]]
p
ggsave(p,file='res_pic/main_figure3/adult_recept_ligand_target_gene_heatmap_score.pdf',height =10 ,width = 6)
saveRDS(adult_ligand_target_res_list,file = 'res_data/fa_adult_ligand_target_res_list.rds')

#
an_df=data.frame(adult_ligand_target_res_list[[1]]$data[adult_ligand_target_res_list[[1]]$data$y %in% c('EPO','CXCL12'),])
an_df=data.frame(row.names = an_df$y,aupr_score=an_df$score)
p=pheatmap(adult_ligand_target_res_list[[3]][,c('EPO','CXCL12')],annotation_col = an_df,border_color = 'white',cluster_rows = F,cluster_cols =F ,color =colorRampPalette(colors =  c("white" ,"#EF6548"))(100))
ggsave(as.ggplot(p),file='res_pic/main_figure3/adult_key_ligand_target_genes_heatmap.pdf',width = 4,height = 8,dpi = 300)


fa_adult_LR_targetgene_enrichGO_res=enrichGO(gene =rownames(adult_ligand_target_res_list[[3]]),OrgDb = org.Hs.eg.db,keyType = 'SYMBOL',ont = 'BP' )
saveRDS(fa_adult_LR_targetgene_enrichGO_res,file = 'res_data/fa_adult_LR_targetgene_enrichGO_res.rds')
cnetplot(fa_adult_LR_targetgene_enrichGO_res,20)
p=draw_gonetcwork_pic_func(res =fa_adult_LR_targetgene_enrichGO_res,showCategory = 20 ,xlimits = c(-0.5,2.5))
ggsave(p,filename='res_pic/main_figure3/adult_target_gene_enrichGO_res.pdf',width =6,height = 6)




cho_adult_ligand_target_df=adult_ligand_target_res_list[[3]][,c('EPO','CXCL12')]
cho_adult_ligand_target_list=list()
cho_adult_ligand_target_list[['EPO']]=names(cho_adult_ligand_target_df[,'EPO'])[cho_adult_ligand_target_df[,'EPO'] >0.04]
cho_adult_ligand_target_list[['CXCL12']]=names(cho_adult_ligand_target_df[,'CXCL12'])[cho_adult_ligand_target_df[,'CXCL12'] >0.04]
key_fa_adult_LR_targetgene_enrichGO_res=enrichGO(gene =unique(unlist(cho_adult_ligand_target_list)),OrgDb = org.Hs.eg.db,keyType = 'SYMBOL',ont = 'BP' )
p=draw_gonetcwork_pic_func(res =key_fa_adult_LR_targetgene_enrichGO_res,showCategory = 20 )
ggsave(p,filename='res_pic/main_figure5/key_fa_fetal_hDEG_targetgene_expression.pdf',width = 6,height = 6)



###############################################################################################################################
# ---------------obtain the specific transmembrane genes involved in cc communication network----------------#
###############################################################################################################################

dir.create('res_pic/main_figure5')

# hDEGs the receptors or ligands(hDEG-LRs) from specific interaction；
receptor_genes=unique(unlist(strsplit(unique(cho_Ery_recepted_LR_df$receptor),split = '_')))
hDEG_receptor=fetal_adult_NRBC_whole_marker[fetal_adult_NRBC_whole_marker$gene %in% receptor_genes & fetal_adult_NRBC_whole_marker$avg_log2FC >1 & fetal_adult_NRBC_whole_marker$pct.1 >0.05 & fetal_adult_NRBC_whole_marker$pct.2<0.1,]
hDEG_receptor_early_df=sub_fetal_adult_all_Ery_tissue_markers[sub_fetal_adult_all_Ery_tissue_markers$gene %in% receptor_genes & sub_fetal_adult_all_Ery_tissue_markers$avg_log2FC >1 & sub_fetal_adult_all_Ery_tissue_markers$pct.2 <0.1 & sub_fetal_adult_all_Ery_tissue_markers$pct.1 >0.1,]
hDEG_receptor_early_df=hDEG_receptor_early_df[hDEG_receptor_early_df$celltype=='early_Ery',]
hDEG_receptor_early_df=hDEG_receptor_early_df[hDEG_receptor_early_df$gene %in% hDEG_receptor$gene,]


ligand_genes=unique(unlist(strsplit(unique(cho_Ery_outgoing_LR_df$ligand),split = '_')))
hDEG_ligand_df=fetal_adult_NRBC_whole_marker[fetal_adult_NRBC_whole_marker$gene %in% ligand_genes & fetal_adult_NRBC_whole_marker$avg_log2FC >1 & fetal_adult_NRBC_whole_marker$pct.2 <0.1 & fetal_adult_NRBC_whole_marker$pct.1 >0.05,]
hDEG_ligand_early_df=sub_fetal_adult_all_Ery_tissue_markers[sub_fetal_adult_all_Ery_tissue_markers$gene %in% ligand_genes & sub_fetal_adult_all_Ery_tissue_markers$avg_log2FC >1 & sub_fetal_adult_all_Ery_tissue_markers$pct.2 <0.1 & sub_fetal_adult_all_Ery_tissue_markers$pct.1 >0.1,]
hDEG_ligand_early_df=hDEG_ligand_early_df[hDEG_ligand_early_df$celltype=='early_Ery',]
hDEG_ligand_early_df=hDEG_ligand_early_df[hDEG_ligand_early_df$gene %in% hDEG_ligand_df$gene,]

# 主要在早期，采用early stage DE 结果
cho_hDEGs_LR_early_df=rbind(hDEG_receptor_early_df,hDEG_ligand_early_df)
cho_hDEGs_LR_early_df=cho_hDEGs_LR_early_df[!duplicated(cho_hDEGs_LR_early_df$gene),]
cho_hDEGs_LR_early_df=cho_hDEGs_LR_early_df[order(cho_hDEGs_LR_early_df$cluster,cho_hDEGs_LR_early_df$avg_log2FC,decreasing = T),]
cho_hDEGs_LR_early_df$gene=factor(cho_hDEGs_LR_early_df$gene,levels = cho_hDEGs_LR_early_df$gene)



# 获得specific fetal target genes
hDEGs_fetal_target_df=fetal_adult_NRBC_whole_marker[fetal_adult_NRBC_whole_marker$gene %in% rownames(fetal_ligand_target_res_list[[3]]) & fetal_adult_NRBC_whole_marker$avg_log2FC>1 & fetal_adult_NRBC_whole_marker$pct.1>0.05 & fetal_adult_NRBC_whole_marker$pct.2<0.1,]
hDEGs_fetal_target_early_df=sub_fetal_adult_all_Ery_tissue_markers[sub_fetal_adult_all_Ery_tissue_markers$gene %in% rownames(fetal_ligand_target_res_list[[3]]) & sub_fetal_adult_all_Ery_tissue_markers$avg_log2FC >1 & sub_fetal_adult_all_Ery_tissue_markers$pct.2 <0.1 & sub_fetal_adult_all_Ery_tissue_markers$pct.1 >0.1,]
hDEGs_fetal_target_early_df=hDEGs_fetal_target_early_df[hDEGs_fetal_target_early_df$celltype=='early_Ery',]
hDEGs_fetal_target_early_df=hDEGs_fetal_target_early_df[hDEGs_fetal_target_early_df$gene %in% hDEGs_fetal_target_df$gene,]

hDEGs_adult_target_df=fetal_adult_NRBC_whole_marker[fetal_adult_NRBC_whole_marker$gene %in% rownames(adult_ligand_target_res_list[[3]]) & fetal_adult_NRBC_whole_marker$avg_log2FC>1 & fetal_adult_NRBC_whole_marker$pct.1>0.05 & fetal_adult_NRBC_whole_marker$pct.2<0.1,]
hDEGs_adult_target_early_df=sub_fetal_adult_all_Ery_tissue_markers[sub_fetal_adult_all_Ery_tissue_markers$gene %in% rownames(adult_ligand_target_res_list[[3]]) & sub_fetal_adult_all_Ery_tissue_markers$avg_log2FC >1 & sub_fetal_adult_all_Ery_tissue_markers$pct.2 <0.1 & sub_fetal_adult_all_Ery_tissue_markers$pct.1 >0.1,]
hDEGs_adult_target_early_df=hDEGs_adult_target_early_df[hDEGs_adult_target_early_df$celltype=='early_Ery',]
hDEGs_adult_target_early_df=hDEGs_adult_target_early_df[hDEGs_adult_target_early_df$gene %in% hDEGs_adult_target_df$gene,]
hDEGs_adult_target_early_df

hDEGs_target_gene_df=rbind(hDEGs_fetal_target_early_df,hDEGs_adult_target_early_df)

CC_hDEGs_df=rbind(cho_hDEGs_LR_early_df,hDEGs_target_gene_df)
# DLK1 & ANXA1 同为ligand +target
CC_hDEGs_df=CC_hDEGs_df[!duplicated(CC_hDEGs_df$gene),]
CC_hDEGs_df=CC_hDEGs_df[order(CC_hDEGs_df$cluster,CC_hDEGs_df$avg_log2FC,decreasing = T),]
CC_hDEGs_df$gene=factor(CC_hDEGs_df$gene,levels = CC_hDEGs_df$gene)

CC_hDEGs_df$transmembrane='un'
CC_hDEGs_df[CC_hDEGs_df$gene %in% unique(unlist(strsplit(CellChatDB.human$interaction$receptor.symbol[CellChatDB.human$interaction$receptor.transmembrane],split = ', '))),'transmembrane']='yes'
CC_hDEGs_df[CC_hDEGs_df$gene %in% unique(unlist(strsplit(CellChatDB.human$interaction$ligand.symbol[CellChatDB.human$interaction$ligand.transmembrane],split = ', '))),'transmembrane']='yes'
CC_hDEGs_df[CC_hDEGs_df$gene %in% unique(unlist(strsplit(CellChatDB.human$interaction$receptor.symbol[!CellChatDB.human$interaction$receptor.transmembrane],split = ', '))),'transmembrane']='no'
CC_hDEGs_df[CC_hDEGs_df$gene %in% unique(unlist(strsplit(CellChatDB.human$interaction$ligand.symbol[!CellChatDB.human$interaction$ligand.transmembrane],split = ', '))),'transmembrane']='no'

CC_hDEGs_df[CC_hDEGs_df$transmembrane=='un',]
CC_hDEGs_df[CC_hDEGs_df$transmembrane=='no',]
CC_hDEGs_df[CC_hDEGs_df$gene %in% c('CD69','SLC6A9','HLA-DQB1','HLA-DRB5','CLEC2B'),'transmembrane']='yes'
CC_hDEGs_df[CC_hDEGs_df$gene %in% c('TNFSF13B','SEMA7A','CISH','IFI16'),'transmembrane']='no'
CC_hDEGs_df[CC_hDEGs_df$transmembrane=='no',]
CC_hDEGs_df[CC_hDEGs_df$transmembrane=='yes',]

saveRDS(CC_hDEGs_df,file = 'res_data/CC_hDEGs_df.rds')

p1=ggplot(CC_hDEGs_df[CC_hDEGs_df$transmembrane=='yes',],aes(x=gene,y=avg_log2FC,fill=cluster))+geom_bar(stat = 'identity')+theme_classic()+theme(axis.text.x = element_text(angle = 90,hjust = 1))+NoLegend()+
  geom_hline(yintercept = 1,linetype = "dashed", color = "red", linewidth = 0.4)+NoLegend()

p2=DotPlot(filt_NBRC_altas_seu,group.by = 'source_celltype',features =CC_hDEGs_df[CC_hDEGs_df$transmembrane=='yes','gene'],cols = c('white','firebrick3'))+RotatedAxis()

p=p1+p2+plot_layout(ncol = 1,heights = c(0.6,1.2));p # 后续添加上taget genes 中同类满足条件的基因
ggsave(p,filename='res_pic/main_figure5/specific_transmenbrane_profile_expression.pdf',width=8,height = 8)


########################################################################################################################################################################
###################----------------------------------------------the signaling pathway -------------------------------------###################
########################################################################################################################################################################

sig_network <- readRDS('signaling_network_human_21122021.rds')
gr_network <- readRDS('gr_network_human_21122021.rds')
ligand_tf_matrix <- readRDS('ligand_tf_matrix_nsga2r_final.rds')


#---------------------------------fetal ligand target path--------------------------------------------------#
cho_fetal_ligand_target_df=fetal_ligand_target_res_list[[3]][,c('IGF1','IGF2','IFNG','TNF','LTA','LTB','EPO')]
cho_fetal_ligand_target_list=list()
cho_fetal_ligand_target_list[['IGF1']]=names(cho_fetal_ligand_target_df[,'IGF1'])[cho_fetal_ligand_target_df[,'IGF1'] >0.04]
cho_fetal_ligand_target_list[['IGF2']]=names(cho_fetal_ligand_target_df[,'IGF2'])[cho_fetal_ligand_target_df[,'IGF2'] >0.04]
cho_fetal_ligand_target_list[['IFNG']]=names(cho_fetal_ligand_target_df[,'IGF1'])[cho_fetal_ligand_target_df[,'IFNG'] >0.04]
cho_fetal_ligand_target_list[['TNF']]=names(cho_fetal_ligand_target_df[,'TNF'])[cho_fetal_ligand_target_df[,'TNF'] >0.04]
cho_fetal_ligand_target_list[['LTA']]=names(cho_fetal_ligand_target_df[,'LTA'])[cho_fetal_ligand_target_df[,'LTA'] >0.04]
cho_fetal_ligand_target_list[['LTB']]=names(cho_fetal_ligand_target_df[,'LTB'])[cho_fetal_ligand_target_df[,'LTB'] >0.04]
cho_fetal_ligand_target_list[['EPO']]=names(cho_fetal_ligand_target_df[,'EPO'])[cho_fetal_ligand_target_df[,'EPO'] >0.04]

ligands_oi <- c("EPO")
targets_oi <- cho_fetal_ligand_target_list[['EPO']]
fetal_EPO_screated_network=signaling_pic_func(ligands_oi =ligands_oi,targets_oi =targets_oi,weighted_networks =weighted_networks,ligand_tf_matrix =ligand_tf_matrix,top_n_regulators = 4    )
# EPO_signaling_fetal_network.pdf,5 X 5 

ligands_oi <- c("IGF1") # this can be a list of multiple ligands if required
targets_oi <- unique(c(as.character(unlist(cho_fetal_ligand_target_list[c("IGF1")]))))
IGF1_fetal_screated_network=signaling_pic_func(ligands_oi =ligands_oi,targets_oi =targets_oi,weighted_networks =weighted_networks,ligand_tf_matrix =ligand_tf_matrix,top_n_regulators = 4    )
# IGF1_signaling_network.pdf, 6x6


ligands_oi <- c("IGF2") # this can be a list of multiple ligands if required
targets_oi <- unique(c(as.character(unlist(cho_fetal_ligand_target_list[c("IGF2")]))))
IGF2_fetal_screated_network=signaling_pic_func(ligands_oi =ligands_oi,targets_oi =targets_oi,weighted_networks =weighted_networks,ligand_tf_matrix =ligand_tf_matrix,top_n_regulators = 4    )
# IGF2_signaling_network.pdf, 5x5


ligands_oi <- c("IFNG") # this can be a list of multiple ligands if required
targets_oi <- unique(c(as.character(unlist(cho_fetal_ligand_target_list[c("IFNG")]))))
IFNG_fetal_screated_network=signaling_pic_func(ligands_oi =ligands_oi,targets_oi =targets_oi,weighted_networks =weighted_networks,ligand_tf_matrix =ligand_tf_matrix,top_n_regulators = 4    )
# IFNG_signaling_network.pdf, 6x6

ligands_oi <- c("TNF") # this can be a list of multiple ligands if required
targets_oi <- unique(c(as.character(unlist(cho_fetal_ligand_target_list[c("TNF")]))))
TNF_fetal_screated_network=signaling_pic_func(ligands_oi =ligands_oi,targets_oi =targets_oi,weighted_networks =weighted_networks,ligand_tf_matrix =ligand_tf_matrix,top_n_regulators = 4    )
# TNF_signaling_network.pdf, 6x6


ligands_oi <- c("LTA") # this can be a list of multiple ligands if required
targets_oi <- unique(c(as.character(unlist(cho_fetal_ligand_target_list[c("LTA")]))))
LTA_fetal_screated_network=signaling_pic_func(ligands_oi =ligands_oi,targets_oi =targets_oi,weighted_networks =weighted_networks,ligand_tf_matrix =ligand_tf_matrix,top_n_regulators = 4    )
# LTA_signaling_network.pdf, 5x5

ligands_oi <- c("LTB") # this can be a list of multiple ligands if required
targets_oi <- unique(c(as.character(unlist(cho_fetal_ligand_target_list[c("LTB")]))))
LTB_fetal_screated_network=signaling_pic_func(ligands_oi =ligands_oi,targets_oi =targets_oi,weighted_networks =weighted_networks,ligand_tf_matrix =ligand_tf_matrix,top_n_regulators = 4    )
# LTB_signaling_network.pdf, 5x5


#---------------------------------adult ligand target path--------------------------------------------------#
cho_adult_ligand_target_df=adult_ligand_target_res_list[[3]][,c('EPO','CXCL12')]
cho_adult_ligand_target_list=list()
cho_adult_ligand_target_list[['EPO']]=names(cho_adult_ligand_target_df[,'EPO'])[cho_adult_ligand_target_df[,'EPO'] >0.04]
cho_adult_ligand_target_list[['CXCL12']]=names(cho_adult_ligand_target_df[,'CXCL12'])[cho_adult_ligand_target_df[,'CXCL12'] >0.04]

ligands_oi <- 'EPO' # this can be a list of multiple ligands if required
targets_oi <- cho_adult_ligand_target_list[['EPO']]
adult_EPO_screated_network=signaling_pic_func(ligands_oi =ligands_oi,targets_oi =targets_oi,weighted_networks =weighted_networks,ligand_tf_matrix =ligand_tf_matrix,top_n_regulators = 4    )
# EPO_signaling_adult_network.pdf,5 X 5 


targets_oi <-unique(c( cho_adult_ligand_target_list[['EPO']], cho_fetal_ligand_target_list[['EPO']]))
fetal_adulkt_EPO_screated_network=signaling_pic_func(ligands_oi =ligands_oi,targets_oi =targets_oi,weighted_networks =weighted_networks,ligand_tf_matrix =ligand_tf_matrix,top_n_regulators = 4    )
# CXCL12_signaling_fetal_adult_network.pdf

ligands_oi <- 'CXCL12' # this can be a list of multiple ligands if required
targets_oi <- cho_adult_ligand_target_list[['CXCL12']]
CXCL12_screated_network=signaling_pic_func(ligands_oi =ligands_oi,targets_oi =targets_oi,weighted_networks =weighted_networks,ligand_tf_matrix =ligand_tf_matrix,top_n_regulators = 4    )
# CXCL12_signaling_network.pdf ,5 x5 

# 分析 不同ligand signal的Signaling mediators
Signaling_mediator_genelist=list()
Signaling_mediator_genelist[['CXCL12']]=unique(CXCL12_screated_network$data$label)
Signaling_mediator_genelist[['CXCL12']]=Signaling_mediator_genelist[['CXCL12']][ !Signaling_mediator_genelist[['CXCL12']] %in% c(cho_adult_ligand_target_list[['CXCL12']],'CXCL12') ]
Signaling_mediator_genelist[['EPO']]=unique(fetal_adulkt_EPO_screated_network$data$label)
Signaling_mediator_genelist[['EPO']]=Signaling_mediator_genelist[['EPO']][ !Signaling_mediator_genelist[['EPO']] %in% unique(c( cho_adult_ligand_target_list[['EPO']],'EPO', cho_fetal_ligand_target_list[['EPO']]))]

ligand='IGF1'
Signaling_mediator_genelist[[ligand]]=unique(IGF1_fetal_screated_network$data$label)
Signaling_mediator_genelist[[ligand]]=Signaling_mediator_genelist[[ligand]][ !Signaling_mediator_genelist[[ligand]] %in% c(cho_fetal_ligand_target_list[[ligand]],'IGF1')]
ligand='IGF2'
Signaling_mediator_genelist[[ligand]]=unique(IGF1_fetal_screated_network$data$label)
Signaling_mediator_genelist[[ligand]]=Signaling_mediator_genelist[[ligand]][ !Signaling_mediator_genelist[[ligand]] %in% c(cho_fetal_ligand_target_list[[ligand]],'IGF2')]
ligand='IFNG'
Signaling_mediator_genelist[[ligand]]=unique(IFNG_fetal_screated_network$data$label)
Signaling_mediator_genelist[[ligand]]=Signaling_mediator_genelist[[ligand]][ !Signaling_mediator_genelist[[ligand]] %in% c(cho_fetal_ligand_target_list[[ligand]],'IFNG')]
ligand='TNF'
Signaling_mediator_genelist[[ligand]]=unique(TNF_fetal_screated_network$data$label)
Signaling_mediator_genelist[[ligand]]=Signaling_mediator_genelist[[ligand]][ !Signaling_mediator_genelist[[ligand]] %in% c(cho_fetal_ligand_target_list[[ligand]],'TNF')]
#ligand='LTA'
#Signaling_mediator_genelist[[ligand]]=unique(LTA_fetal_screated_network$data$label)
#Signaling_mediator_genelist[[ligand]]=Signaling_mediator_genelist[[ligand]][ !Signaling_mediator_genelist[[ligand]] %in% c(cho_fetal_ligand_target_list[[ligand]],'LTA)]
ligand='LTB'
Signaling_mediator_genelist[[ligand]]=unique(LTB_fetal_screated_network$data$label)
Signaling_mediator_genelist[[ligand]]=Signaling_mediator_genelist[[ligand]][ !Signaling_mediator_genelist[[ligand]] %in% c(cho_fetal_ligand_target_list[[ligand]],'TLB')]
length(unique(unlist(Signaling_mediator_genelist)))
all_genes=unique(unlist(Signaling_mediator_genelist))
saveRDS(Signaling_mediator_genelist,file = 'Signaling_mediator_genelist.rds')

Reduce(intersect, Signaling_mediator_genelist) # "RELA"  "STAT3", TP53 , 所有的交集

# 转换为二进制矩阵
all_genes <- unique(unlist(Signaling_mediator_genelist))
binary_mat <- sapply(Signaling_mediator_genelist, function(x) as.integer(all_genes %in% x))
rownames(binary_mat) <- all_genes
p=pheatmap(t(as.data.frame(binary_mat)),color = c('white','firebrick3'))
ggsave(as.ggplot(p),file='res_pic/main_figure3/Signaling_mediator_distribution.pdf',width = 15,height = 3)

focused_genes=c('MYC','ESR1','MAPK1','MAPK8','TP53','STAT3','STAT1','STAT5A','STAT5B','RELA','NFKB1','EP300','JUN','JUND','FOS')# ESR1 几乎不表达
p=VlnPlot(filt_NBRC_altas_seu,group.by = 'source_celltype',features = focused_genes,stack = T,cols = cols)+NoLegend();p
ggsave(p,file='res_pic/main_figure3/foucused_Signaling_mediator_expression_vlnplot.pdf',width = 10,height = 6)

p=DotPlot(subset(filt_NBRC_altas_seu,tissue_stage!='YS'),group.by = 'source_celltype',features = focused_genes,scale = F,cols = c('white','firebrick3'))+RotatedAxis()

sub_fetal_adult_all_Ery_tissue_markers[sub_fetal_adult_all_Ery_tissue_markers$gene %in% focused_genes & sub_fetal_adult_all_Ery_tissue_markers$avg_log2FC >0 & sub_fetal_adult_all_Ery_tissue_markers$celltype=='early_Ery',]


VlnPlot(subset(filt_NBRC_altas_seu,tissue_stage!='YS'),group.by = 'source_celltype',features = c('MYC','ESRR1','MAPK1','TP53','RELA','STAT3','NFKB1'),stack = T)+NoLegend()

all_genes=all_genes[!all_genes %in% c('CXCL12','IGF1','IGF2','EPO','LTA','LTB','IFNG','TNF')]
p=DotPlot(object = filt_NBRC_altas_seu,features =unique(unlist(Signaling_mediator_genelist)),group.by = 'source_celltype',scale = F )
temp_df=p$data
temp_df=temp_df[-grep('YS',temp_df$id),]
rownames(temp_df)=NULL



temp_df=temp_df[temp_df$avg.exp >0.1 & temp_df$pct.exp>10,]
length(unique(temp_df$features.plot))
ggplot(temp_df,aes(y=id,x=features.plot,size=pct.exp,color=avg.exp))+geom_point()+theme_classic()+scale_color_gradient(low = 'white',high = 'firebrick3')+RotatedAxis()

order_Signaling_mediator_genes=c(unique(temp_df$features.plot)[unique(temp_df$features.plot) %in% Signaling_mediator_genelist[['EPO']]],
                                 unique(temp_df$features.plot)[unique(temp_df$features.plot) %in%    Signaling_mediator_genelist[['IGF1']]],
                                 unique(temp_df$features.plot)[unique(temp_df$features.plot) %in% Signaling_mediator_genelist[['IFNG']]],
                                 unique(temp_df$features.plot)[unique(temp_df$features.plot) %in% Signaling_mediator_genelist[['TNF']]],
                                 unique(temp_df$features.plot)[unique(temp_df$features.plot) %in% Signaling_mediator_genelist[['IGF2']]],
                                 unique(temp_df$features.plot)[unique(temp_df$features.plot) %in% Signaling_mediator_genelist[['CXCL12']]]
)

order_Signaling_mediator_genes=as.character(unique(order_Signaling_mediator_genes));length(order_Signaling_mediator_genes)
p1=DotPlot(object = filt_NBRC_altas_seu,features =unique(temp_df$features.plot),group.by = 'source_celltype',scale = F,cols = c('gray','firebrick3') )+RotatedAxis()






########################################################################################################################################################################
#---------------------- cell source of EPO、CXCL12、IGF1 、IGF2、IFNG、TNF(LTA\LTB) expression in niches ------------------------------------#
########################################################################################################################################################################

dir.create('res_pic/main_figure4')

unique(NRBC_altas_LR_df[grep('IGF1',NRBC_altas_LR_df$ligand),'receptor']) # IGF1R
unique(NRBC_altas_LR_df[grep('IFNG',NRBC_altas_LR_df$ligand),'receptor']) # IFNGR1_IFNGR2
unique(NRBC_altas_LR_df[grep('TNF',NRBC_altas_LR_df$ligand),'receptor']) #  "TNFRSF1A"  "LTBR"      "TNFRSF21"  "TNFRSF13B" "TNFRSF13C" "TNFRSF17"
unique(NRBC_altas_LR_df[grep('LTA',NRBC_altas_LR_df$ligand),'receptor']) #  "TNFRSF1A"  "LTBR"   
unique(NRBC_altas_LR_df[grep('LTB',NRBC_altas_LR_df$ligand),'receptor']) #  LTB4R
unique(NRBC_altas_LR_df[grep('CXCL12',NRBC_altas_LR_df$ligand),'receptor']) #  CXCR4


p=DotPlot(filt_NBRC_altas_seu,scale=F,features = c('EPOR','CXCR4','IFNGR1','IFNGR2','IGF1R','TNFRSF1A',"LTBR"),group.by = 'source_celltype')+
  RotatedAxis()+scale_color_gradient(low = 'gray',high = 'firebrick3');p
ggsave(p,filename='res_pic/main_figure4/the_key_ligand_receptor_expression_in_definitive.pdf',width =6 ,height = 6)

#-----------------------------------------check the key secreted genes expression in FL/FBM/ABM---------------------#
FL_altas_seu=readRDS('../NRBC_FL_altas/tmp_FL_altas_seu.rds')

cho_celtype=c( "MONOCYTE","MACROPHAGE", "NK/T CELLS",'B CELLS','LMPP_MLP','DC' ,'ILC',"HEPATOCYTE","FIBROBLASTS","SMOOTH MUSCLE","SKELETAL MUSCLE","MESOTHELIUM","NEPHRON" )
cho_feature= c('EPO','CXCL12','IGF1','IGF2','IFNG','TNF','LTA','LTB')
p1=VlnPlot(subset(FL_altas_seu,subcelltype %in% cho_celtype),group.by = 'subcelltype',features =cho_feature,stack = T)+NoLegend()+ggtitle('FL ALTAS')
ggsave(p1,filename='res_pic/main_figure4/FL_key_ligand_expression_celltype_niche_vlnplot.pdf',width = 5,height = 5)


FL_altas_seu$age=factor(FL_altas_seu$age,levels = c("CS14_4PCW" ,"CS15_5PCW", "CS17_6PCW", "CS18","CS22","CS23",  "8PCW","8.1PCW","9.1PCW","9.7PCW",  "11PCW","11.4PCW" ,"12PCW","13.9PCW","14.4PCW","15PCW","16.3PCW","16PCW","17PCW"  ))
p=VlnPlot(subset(FL_altas_seu,subcelltype %in% c('MACROPHAGE')),group.by = 'age',features =cho_feature[-1],stack = T)+NoLegend()+ggtitle('FL  MACROPHAGE');p
ggsave(p,filename='res_pic/main_figure4/FL_Mac_key_ligand_expression_vlnplot.pdf',height = 6,width = 4)

# Mac 类别：IGF1+IGF2+，IGF1+TNF+, IGF1+TNF+CXCL12+,IGF1+CXCL12+
{
  
  macrophages <- subset(FL_altas_seu, subcelltype == "MACROPHAGE")
  
  macrophages$IGF1_pos <- FetchData(macrophages, "IGF1")[,1] > 0
  macrophages$IGF2_pos <- FetchData(macrophages, "IGF2")[,1] > 0
  macrophages$TNF_pos <- FetchData(macrophages, "TNF")[,1] > 0
  macrophages$CXCL12_pos <- FetchData(macrophages, "CXCL12")[,1] > 0
  
  macrophages@meta.data$group <- "other"
  macrophages@meta.data$group[macrophages@meta.data$IGF1_pos & macrophages@meta.data$IGF2_pos & !macrophages@meta.data$CXCL12_pos & !macrophages@meta.data$TNF_pos] <- "IGF1+IGF2+"
  macrophages@meta.data$group[macrophages@meta.data$IGF1_pos & macrophages@meta.data$TNF_pos & macrophages@meta.data$CXCL12_pos] <- "IGF1+TNF+CXCL12+"
  macrophages@meta.data$group[macrophages@meta.data$IGF1_pos & macrophages@meta.data$TNF_pos & !macrophages@meta.data$CXCL12_pos] <- "IGF1+TNF+"
  macrophages@meta.data$group[!macrophages@meta.data$IGF1_pos & macrophages@meta.data$TNF_pos & !macrophages@meta.data$CXCL12_pos] <- "TNF+"
  macrophages@meta.data$group[!macrophages@meta.data$IGF1_pos & macrophages@meta.data$CXCL12_pos & !macrophages@meta.data$TNF_pos ] <- "CXCL12+"
  macrophages@meta.data$group[macrophages@meta.data$IGF1_pos & macrophages@meta.data$CXCL12_pos & !macrophages@meta.data$TNF_pos] <- "IGF1+CXCL12+"
  macrophages@meta.data$group[!macrophages@meta.data$IGF1_pos & macrophages@meta.data$CXCL12_pos & macrophages@meta.data$TNF_pos ] <- "TNF+CXCL12+"
  
  count_df <- macrophages@meta.data %>%
    filter(group != "other") %>%
    group_by(age, group) %>%
    summarise(n = n(), .groups = "drop")
  
  count_df$group=factor(count_df$group,levels = c( "IGF1+IGF2+", "IGF1+TNF+CXCL12+","IGF1+TNF+","IGF1+CXCL12+","TNF+CXCL12+","TNF+",'CXCL12+'))
  p=ggplot(count_df, aes(x = age, y = n, color = group, group = group)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    theme_classic() +scale_color_manual(values = cols[-2])+
    labs(x = "Age", y = "Number of macrophages", color = "Subset")+RotatedAxis()
 p   
 ggsave(p,filename='res_pic/main_figure4/FL_Mac_key_ligand_expression_number_along_age.pdf',height = 8,width = 4)
 

 
 VlnPlot(subset(macrophages,anno_lvl_2_final_clean %in% c('MACROPHAGE_IRON_RECYCLING','MACROPHAGE_KUPFFER_LIKE','MACROPHAGE_LYVE1_HIGH','MACROPHAGE_MHCII_HIGH','MACROPHAGE_PROLIFERATING')),
         group.by = 'age',features = c('IGF1','IGF2','CXCL12','TNF'),stack = T,split.by = 'anno_lvl_2_final_clean')
 
 macrophages$VCAM1_pos <- FetchData(macrophages, "VCAM1")[,1] > 1
 VlnPlot(macrophages,group.by = 'VCAM1_pos',features = c('IGF1','IGF2','IGF2','CXCL12','TNF'),stack = T)

 central_mac_markers <- list( Central_Macrophage = c( "VCAM1", "CD163", "HMOX1", "SLC40A1", "SPIC")) # 
  
 library(UCell) 
 macrophages <- AddModuleScore_UCell(obj = macrophages,features = central_mac_markers,name = "_UCell",ncores = 4)
 summary(macrophages$Central_Macrophage_UCell)
 # 查看分布
 ggplot(macrophages@meta.data, aes(x = Central_Macrophage_UCell)) +
   geom_density(fill = "grey70", alpha = 0.5) +
   geom_rug(alpha = 0.1) +
   theme_classic() +scale_x_continuous(breaks = seq(0, 1, 0.03), limits = c(0, 0.9)) 
   labs(x = "Central macrophage UCell score", y = "Density")
 
 threshold <-0.62
 
 macrophages$Central_Macrophage <- ifelse(macrophages$Central_Macrophage_UCell >= threshold,"Central_Macrophage","Other")
 p=VlnPlot(macrophages,group.by = 'Central_Macrophage',features = c('IGF1','IGF2','CXCL12','TNF'),stack = T,pt.size = 0,flip = T)+NoLegend();p
 ggsave(p,filename='res_pic/main_figure4/FL_Central_Mac_key_ligand_expression.pdf',height = 6,width = 3)
 
 p=VlnPlot(macrophages,group.by = 'age',split.by  = 'Central_Macrophage',features = c('IGF1','IGF2','CXCL12','TNF'),stack = T,pt.size = 0)+NoLegend();p
 ggsave(p,filename='res_pic/main_figure4/FL_Central_Mac_key_ligand_expression_along_age.pdf',height = 8,width = 4)

 
 plot_data <- FetchData(macrophages, vars = c("IGF1",'IGF2', "CXCL12", "TNF", "Central_Macrophage"))
 t.test(IGF1 ~ Central_Macrophage, data = plot_data ) # ***
 t.test(IGF2 ~ Central_Macrophage, data = plot_data ) # ***,
 t.test(CXCL12 ~ Central_Macrophage, data = plot_data )# ***
 t.test(TNF ~ Central_Macrophage, data = plot_data )# na
 summary(plot_data[plot_data$Central_Macrophage=='Other','TNF'])
 summary(plot_data[plot_data$Central_Macrophage=='Central_Macrophage','TNF'])
}

p=VlnPlot(subset(FL_altas_seu,subcelltype %in% "NK/T CELLS"),group.by = 'age',features =c('IFNG'),pt.size = 0)+NoLegend()+ggtitle('FL  NT/T CELLS: IFNG');p
# # IFNG 持续存在
p1=VlnPlot(subset(FL_altas_seu,subcelltype %in% c('MONOCYTE')),group.by = 'age',features =c('TNF'),pt.size = 0)+NoLegend()+ggtitle('FL  MONOCYTE:TNF ');p1
#TNF 持续存在
p=p/p1;p
ggsave(p,filename='res_pic/main_figure4/FL_MONOCYTE_NKT_key_ligand_expression_vlnplot.pdf',height = 6,width = 8)

# all-- > MONOCYTE_III_IL1  TNF+
p=VlnPlot(subset(FL_altas_seu,anno_lvl_2_final_clean %in% c('MONOCYTE_I_CXCR4','MONOCYTE_II_CCR2','MONOCYTE_III_IL1B','PROMONOCYTE')),
          split.by = 'anno_lvl_2_final_clean',group.by = 'age',features =cho_feature[-1],stack = T)+ggtitle('FL  MONOCYTE ')
 
ggsave(p,filename='res_pic/main_figure4/FL_Mono_subcelltype_key_ligand_expression_vlnplot.pdf',height = 8,width = 6)

p=VlnPlot(subset(FL_altas_seu,subcelltype %in% c("HEPATOCYTE")),group.by = 'age',features =c( "CXCL12", "IGF1","IGF2" ,  "TNF"  ,"LTB"   ),stack = T)+NoLegend()+ggtitle('HEPATOCYTE');p
p=VlnPlot(subset(FL_altas_seu,subcelltype %in% c("HEPATOCYTE")),group.by = 'age',features =c('EPO'),pt.size = 0)+NoLegend()+ggtitle('HEPATOCYTE');p
ggsave(p,filename='res_pic/main_figure4/FL_HEPATOCYTE_key_ligand_expression_vlnplot.pdf',height = 3,width =8)
rm(p);gc()

p=VlnPlot(subset(FL_altas_seu,subcelltype %in% c("FIBROBLASTS")),group.by = 'age',features =c( "CXCL12", "IGF1","IGF2" ,  "TNF"  ,"LTB"   ),stack = T)+NoLegend()+ggtitle('FIBROBLASTS');p
p=VlnPlot(subset(FL_altas_seu,subcelltype %in% c("SMOOTH MUSCLE","SKELETAL MUSCLE")),group.by = 'age',features =cho_feature[-1],stack = T)+NoLegend()+ggtitle('MUSCLE');p
p=VlnPlot(subset(FL_altas_seu,subcelltype %in% c("MESOTHELIUM")),group.by = 'age',features =c('CXCL12',"IGF1","IGF2"  , "LTA","LTB" ),stack = T)+NoLegend()+ggtitle('MESOTHELIUM');p

rm(p);gc()
FL_altas_seu$id=FL_altas_seu$donor
FL_altas_seu$id[ !FL_altas_seu$id %in% unique(FL_altas_seu$id)[grep('wk',unique(FL_altas_seu$id))]]=paste(FL_altas_seu$donor[ !FL_altas_seu$id %in% unique(FL_altas_seu$id)[grep('wk',unique(FL_altas_seu$id))]],
                                                                                                          FL_altas_seu$age[ !FL_altas_seu$id %in% unique(FL_altas_seu$id)[grep('wk',unique(FL_altas_seu$id))]],sep='_')
FL_cho_gene_aggregated_exp=AggregateExpression(FL_altas_seu,features = cho_feature,group.by = 'id')$RNA
FL_cho_gene_aggregated_exp=log2(FL_cho_gene_aggregated_exp+1)
colnames(FL_cho_gene_aggregated_exp)=gsub(pattern = 'wk','WPC',colnames(FL_cho_gene_aggregated_exp))

order_sampleid=c( "FL-4WPC" ,"FL-5WPC" ,"FL-6WPC", "F61-CS18","F35-CS22","F32-CS22","F34-CS23", "FL-8WPC"  ,"F16-8.1PCW","F17-9.1PCW","F22-9.7PCW","F33-9.7PCW","FL-11WPC" ,
                  "F23-11.4PCW","F30-14.4PCW","F38-12PCW" , "F45-13.9PCW" ,"F30-14.4PCW", "FL-15WPC"  ,"F41-16PCW","F21-16.3PCW","F29-17PCW")

colnames(FL_cho_gene_aggregated_exp)[!colnames(FL_cho_gene_aggregated_exp) %in% order_sampleid]
FL_cho_gene_aggregated_exp=FL_cho_gene_aggregated_exp[,order_sampleid]

p=pheatmap(FL_cho_gene_aggregated_exp,cluster_rows = F,cluster_cols = F,color = colorRampPalette(colors = c('navy','white','firebrick3'))(100))
ggsave(as.ggplot(p),filename='res_pic/main_figure5/key_ligand_expression_inFL_niche_samples.pdf',width =6,height = 4 )

FL_cho_gene_aggregated_exp_df=t(FL_cho_gene_aggregated_exp)
FL_cho_gene_aggregated_exp_df=data.frame(FL_cho_gene_aggregated_exp_df/FL_cho_gene_aggregated_exp_df[,'IGF1'])
FL_cho_gene_aggregated_exp_df=FL_cho_gene_aggregated_exp_df[,c('IGF1','EPO','CXCL12','IFNG','TNF','IGF2','LTA','LTB')]
pheatmap(FL_cho_gene_aggregated_exp_df,cluster_rows = F)
FL_cho_gene_aggregated_exp_df[FL_cho_gene_aggregated_exp_df==0]=NA
boxplot(FL_cho_gene_aggregated_exp_df,rm.na=T)
round(colMedians(as.matrix(FL_cho_gene_aggregated_exp_df),na.rm = T),digits = 1)
#IGF1    EPO CXCL12   IFNG    TNF   IGF2    LTA    LTB 
#1.0    0.7    1.4    1.2    1.2    1.2    0.9    1.4 
# FL_agrregated_expression_ref_IGF1_boxplot.pdf,6x6 


#-----------------------------------BM---------------------#

BM_altas_seu=readRDS('../NRBC_BM_altas/BM_altas_seu_v2.rds')#
VlnPlot(BM_altas_seu,group.by = 'age',features = cho_feature,stack = T)

table(BM_altas_seu$new_celltype[BM_altas_seu$stage=='EBM'])[table(BM_altas_seu$new_celltype[BM_altas_seu$stage=='EBM']) >0]
EBM_cho_celltype=c("MACROPHAGE", "Myoprogenitor","Osteoprogenitors","CHONDROBLAST",   "CHONDROCYTE","PMSC1","PMSC2","OCPs","BMSC1","BMSC2" )
p0=VlnPlot(subset(subset(BM_altas_seu,stage=='EBM'),new_celltype %in% EBM_cho_celltype ),features = cho_feature[-5],stack = T,group.by = 'new_celltype')+NoLegend()+ggtitle('EBM ALTAS')
ggsave(p0,filename='res_pic/main_figure5/EBM_celltype_key_ligand_expression_vlnplot.pdf',height = 4,width = 6)

levels(BM_altas_seu$new_celltype)

FBM_cho_cells=c("LMPP_MLP","DCs", "MACROPHAGE" ,"ILC","B CELLs" ,"NK/T CELLS", "CHONDROCYTE" , "OSTEOBLAST","FIBROBLASTS","SKELETAL MUSCLE","VSMC","Schwan cell"  )
p2=VlnPlot(subset(subset(BM_altas_seu,stage=='FBM'),new_celltype %in% FBM_cho_cells ),features = cho_feature,stack = T,group.by = 'new_celltype')+NoLegend()+ggtitle('FBM ALTAS')

p=VlnPlot(subset(subset(BM_altas_seu,stage=='FBM'),new_celltype %in% c('MACROPHAGE')),group.by = 'age',features =cho_feature[-1],stack = T)+NoLegend()+ggtitle('FBM  MACROPHAGE');p
ggsave(p,filename='res_pic/main_figure5/FBM_Mac_key_ligand_expression_vlnplot.pdf',height = 6,width =4)
rm(p);gc()

central_mac_markers <- list( Central_Macrophage = c( "VCAM1", "CD163", "HMOX1", "SLC40A1", "SPIC")) # 
macrophages=subset(BM_altas_seu,new_celltype %in% c('MACROPHAGE'))
macrophages <- AddModuleScore_UCell(obj = macrophages,features = central_mac_markers,name = "_UCell",ncores = 4 )
summary(macrophages$Central_Macrophage_UCell)
# 查看分布
ggplot(macrophages@meta.data, aes(x = Central_Macrophage_UCell)) +
  geom_density(fill = "grey70", alpha = 0.5) +
  geom_rug(alpha = 0.1) +
  theme_classic() +scale_x_continuous(breaks = seq(0, 1, 0.03), limits = c(0, 1)) 
labs(x = "Central macrophage UCell score", y = "Density")


threshold <-0.43
macrophages$Central_Macrophage <- ifelse(macrophages$Central_Macrophage_UCell >= threshold,"Central_Macrophage","Other")
table(macrophages@meta.data[,c('Central_Macrophage','stage')])

p=VlnPlot(subset(macrophages,stage=='FBM'),group.by = 'Central_Macrophage',features = c('IGF1','CXCL12'),stack = T,pt.size = 0,flip = T)+NoLegend();p
ggsave(p,filename='res_pic/main_figure4/FBM_Central_Mac_key_ligand_expression.pdf',height = 6,width = 4)

p=VlnPlot(subset(macrophages,stage=='FBM'),group.by = 'age',split.by  = 'Central_Macrophage',features = c('IGF1','CXCL12'),stack = T,pt.size = 0);p
ggsave(p,filename='res_pic/main_figure4/FBM_Central_Mac_key_ligand_expression_along_age.pdf',height = 6,width = 4)

plot_data <- FetchData(subset(macrophages,stage=='FBM'), vars = c("IGF1", "CXCL12", "Central_Macrophage"))
t.test(IGF1 ~ Central_Macrophage, data = plot_data ) # ***
t.test(CXCL12 ~ Central_Macrophage, data = plot_data )# ***


p=VlnPlot(subset(macrophages,stage=='ABM'),group.by = 'Central_Macrophage',features = c('CXCL12','IGF1'),stack = T,pt.size = 0,flip = T)+NoLegend();p
ggsave(p,filename='res_pic/main_figure4/ABM_Central_Mac_key_ligand_expression_number_along_age.pdf',height = 8,width = 4)

plot_data <- FetchData(subset(macrophages,stage=='ABM'), vars = c( "CXCL12", "Central_Macrophage"))
t.test(CXCL12 ~ Central_Macrophage, data = plot_data )# ***



p=VlnPlot(subset(subset(BM_altas_seu,stage=='FBM'),new_celltype %in% c('MONOCYTE')),group.by = 'age',features =cho_feature[-1],stack = T)+NoLegend()+ggtitle('FBM  MONOCYTE');p
ggsave(p,filename='res_pic/main_figure4/FBM_Mono_key_ligand_expression_vlnplot.pdf',height = 6,width =4)

VlnPlot(subset(BM_altas_seu,stage=='FBM'),features = c("IFNG",'TNF','LTA'),stack = T,group.by = 'anno_final_celltype2')+NoLegend()+ggtitle('FBM ALTAS')

# CD16 monocyte 和Proliferation T/NK 表达TNF 
VlnPlot(subset(BM_altas_seu,stage=='ABM'),features = c("IFNG",'TNF','LTA'),stack = T,group.by = 'ct')+NoLegend()+ggtitle('ABM ALTAS')


ABM_cho_cells=c("NEUTROPHIL","MACROPHAGE","Mac_Ery" ,"CLP","B CELLs","Plasma cell", "CD4 T","CD8 T","Treg","NK/T CELLS","OSTEOCLAST",
                "ENDOTHELIUM","Adipo-MSC","APOD+ MSC","Fibro-MSC","Osteo-MSC","THY1+ MSC","RNAlo MSC"  )
p3=VlnPlot(subset(subset(BM_altas_seu,stage=='ABM'),new_celltype %in% ABM_cho_cells),features = cho_feature,stack = T,group.by = 'new_celltype')+NoLegend()+ggtitle('ABM ALTAS')

p=p1+p2+p3;p
ggsave(p,filename='res_pic/main_figure4/key_ligand_expression_celltype_niche_vlnplot.pdf',width = 15,height = 5)


BM_altas_seu$donor[is.na(BM_altas_seu$donor)]=BM_altas_seu$sample[is.na(BM_altas_seu$donor)]
BM_altas_seu$id=paste(BM_altas_seu$donor,BM_altas_seu$age,sep="_")
table(BM_altas_seu$id)
BM_cho_gene_aggregated_exp=AggregateExpression(BM_altas_seu,group.by = 'id',features =cho_feature )$RNA
BM_cho_gene_aggregated_exp=log2(BM_cho_gene_aggregated_exp+1)
p=pheatmap(BM_cho_gene_aggregated_exp[,grep('H|F|CS',colnames(BM_cho_gene_aggregated_exp))],cluster_cols = F,cluster_rows = F)
ggsave(as.ggplot(p),filename='res_pic/main_figure4/key_ligand_expression_in_BM_sample_heatmap.pdf',width =6,height = 4)

BM_cho_gene_aggregated_exp_df=t(BM_cho_gene_aggregated_exp[,grep('H|F',colnames(BM_cho_gene_aggregated_exp))])
BM_cho_gene_aggregated_exp_df=data.frame(BM_cho_gene_aggregated_exp_df/BM_cho_gene_aggregated_exp_df[,'IGF1'])
BM_cho_gene_aggregated_exp_df=BM_cho_gene_aggregated_exp_df[,c('IGF1','CXCL12','IFNG','TNF','IGF2','LTA','LTB')]
BM_cho_gene_aggregated_exp_df[BM_cho_gene_aggregated_exp_df==0]=NA
boxplot(BM_cho_gene_aggregated_exp_df[grep('CW',rownames(BM_cho_gene_aggregated_exp_df)),],rm.na=T)
round(colMedians(as.matrix(BM_cho_gene_aggregated_exp_df[grep('CW',rownames(BM_cho_gene_aggregated_exp_df)),]),na.rm = T),digits = 1)
#CXCL12 IGF1 IFNG  TNF                    LTA  LTB 
#1.4  1.0  0.7  1.3  1.3  1.1  1.9 
# FBM_agrregated_expression_ref_IGF1_boxplot.pdf,7.5 x 5.6

boxplot(BM_cho_gene_aggregated_exp_df[grep('y',rownames(BM_cho_gene_aggregated_exp_df)),],rm.na=T)
round(colMedians(as.matrix(BM_cho_gene_aggregated_exp_df[grep('y',rownames(BM_cho_gene_aggregated_exp_df)),]),na.rm = T),digits = 1)
# ABM_agrregated_expression_ref_IGF1_boxplot.pdf,7.5 x 5.6
#IGF1 CXCL12   IFNG    TNF   IGF2    LTA    LTB 
#1.0    1.2    0.6    0.6    0.8    0.4    1.0 

YS_altas_seu= readRDS('../NRBC_YS_altas/raw_ref_data/dealt_YS_altas_seu_20251028.rds' )
YS_altas_seu=NormalizeData(YS_altas_seu)
levels(YS_altas_seu$subcelltype)
YS_cho_celltype=c("EO/BASO/MAST" ,"MACROPHAGE","DEF_HSPC" ,"LMPP", "MOP","MONOCYTE","MOMO_MAC_DC","ELP","ILC","NK",  
                  "B_CELL", "MONOCYTE_MACROPHAGE" , "MESOTHELIUM","SMOOTH_MUSCLE" ,
                  "ENDODERM","ENDOTHELIUM","FIBROBLAST")
p4=VlnPlot(subset(YS_altas_seu,subcelltype %in% YS_cho_celltype),group.by = 'subcelltype',features =cho_feature,stack = T)+NoLegend()+ggtitle('YS ALTAS')
p4
ggsave(p4,filename='res_pic/main_figure4/key_ligand_expression_celltype_YS_vlnplot.pdf',width = 6,height = 6)

p4=VlnPlot(subset(YS_altas_seu,subcelltype %in% c('FIBROBLAST','SMOOTH_MUSCLE','MACROPHAGE')),cols = cols,group.by = 'stage',split.by = 'subcelltype',
           features =c('CXCL12','IGF1','IGF2','TNF'),stack = T)
p4
ggsave(p4,filename='res_pic/main_figure4/CXCL12_IGF2_expression_stagetime_YS_vlnplot.pdf',height = 6,width = 4)

YS_cho_gene_aggregated_exp=AggregateExpression(YS_altas_seu,features = cho_feature,group.by = 'id')$RNA
YS_cho_gene_aggregated_exp=log2(YS_cho_gene_aggregated_exp+1)
p=pheatmap(YS_cho_gene_aggregated_exp,cluster_rows = F,cluster_cols = F)
ggsave(as.ggplot(p),filename='res_pic/main_figure4/key_ligand_expression_in_YS_sample_heatmap.pdf',height = 4,width = 6)

YS_cho_gene_aggregated_exp_df=t(YS_cho_gene_aggregated_exp)
YS_cho_gene_aggregated_exp_df=data.frame(YS_cho_gene_aggregated_exp_df/YS_cho_gene_aggregated_exp_df[,'IGF1'])
YS_cho_gene_aggregated_exp_df=YS_cho_gene_aggregated_exp_df[,c('IGF1','EPO','CXCL12','IFNG','TNF','IGF2','LTA','LTB')]
boxplot(as.matrix(YS_cho_gene_aggregated_exp_df))
round(colMedians(as.matrix(YS_cho_gene_aggregated_exp_df),na.rm = T),digits = 1)
#  IGF1    EPO CXCL12   IFNG    TNF   IGF2    LTA    LTB 
# 1.0    0.6    0.9    0.1    1.1    1.2    0.5    0.8 
# YS_agrregated_expression_ref_IGF1_boxplot.pdf,7.5 x 5.6

p=VlnPlot(subset(YS_altas_seu,subcelltype %in% c("MOP","MONOCYTE",'MONOCYTE_MACROPHAGE',"MOMO_MAC_DC")),cols = cols,group.by = 'stage',split.by = 'subcelltype',features =cho_feature[-1],stack = T)+ggtitle('YS MONOCYTE')
ggsave(as.ggplot(p),filename='res_pic/main_figure4/key_ligand_expression_celltype_YS_Mono_vlnplot.pdf',height = 6,width = 6)


