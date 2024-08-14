plot_sp_box <- function(otu,meta,sp,pre_cut = 0.2){
  plot_r <- list()
  data_r <- list() 
  j = 1
  for(i in 1:length(sp)){
    if(! sp[i] %in% row.names(otu)){
      print(paste(c("sp is not exist:",sp[i]),collapse = " "))
      next
    }
    tmp_data <- as.data.frame(t(otu[which(row.names(otu)==sp[i]),meta$sample]))
    #print(sp[i])
    tmp_data$sample <- row.names(tmp_data)
    colnames(tmp_data) <- c("value","sample")
    tmp_data_p <- inner_join(meta,tmp_data,by = "sample")
    if(length(which(tmp_data_p$value != 0)) < nrow(tmp_data_p) * pre_cut ){
      print(paste(c("sp have lower prevalence than threshold:",sp[i]),collapse = " "))
      next
    }
    data_r[[j]] <- tmp_data_p
    tmp_c=unique(tmp_data_p$groups)
    tmp_a <- as.data.frame(t(combn(tmp_c,2))) %>%
      mutate(tmp = str_c(V1,",",V2)) %>% 
      select(tmp)
    tmp_list <- (str_split(tmp_a$tmp,pattern = ","))
    plot_r[[j]] <- ggplot(tmp_data_p,aes(x=groups,y=value,color = groups))+
      geom_boxplot(outlier.shape = NA,)+
      geom_jitter(width = 0.2,height = 0)+
      theme_cowplot()+
      theme(
        legend.position = "none"
      )+
      labs(x="Group",y="",title = sp[i])+
      stat_compare_means(comparisons = tmp_list,
                         method = "wilcox.test", label = "p.format")+
      scale_color_manual(values = color0)
    j = j + 1
  }
  res_r <-list(plot_r,data_r) 
  return(res_r)
}
pathway_beta <- function(x,y,z){
  tmp_meta_pca <- y %>% dplyr::select(sample,groups,age,sex)
  out1 <- paste("output/cross",z,"pca.pdf",sep = ".")
  out2 <- paste("output/cross",z,"pcoa.pdf",sep = ".")
  data_path_pca=x
  data_path_pca <- noise_removal(data_path_pca,percent = 0.1,low = 0.001)
  data_path_pca <- as.data.frame(t(data_path_pca))
  data_path_pca$sample <- row.names(data_path_pca)
  data_path_pca <- data_path_pca %>% filter(sample %in% tmp_sam)
  data_path_pca <- inner_join(data_path_pca,tmp_meta_pca,by = "sample")
  row.names(data_path_pca) <- data_path_pca$sample
  
  tmp_meta_pca <- data_path_pca %>% dplyr::select(sample,groups,age,sex)
  data_path_pca <- data_path_pca %>% dplyr::select(!all_of(c("sample","groups","age","sex")))
  
  tmp_dist = vegdist(as.matrix(data_path_pca), method = 'bray')
  tmp_site = as.data.frame(data.frame(sample = tmp_meta_pca$sample,
                                      group = tmp_meta_pca$groups,
                                      age = tmp_meta_pca$age,
                                      sex = tmp_meta_pca$sex))
  tmp_adonis_result_dis = adonis2(data_path_pca~group, tmp_site)
  res_pca <- prcomp(data_path_pca, scale = T)
  
  p1 <- fviz_pca_ind(res_pca, label="none",axes = c(1, 2), alpha.ind =1,
                     habillage=tmp_meta_pca$groups,invisible="quali",pointsize = 2.5, addEllipses = TRUE,
                     ellipse.level=0.66,palette = color0)+
    theme_bw()+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text.x = element_text(vjust = 0.5, hjust = 0.1)
    )+
    annotate("text", x=max(res_pca$x[,1])-25, y=max(res_pca$x[,2])-5,
             label=paste("p-value(adonis)=",tmp_adonis_result_dis$`Pr(>F)`[1]))+
    annotate("text", x=max(res_pca$x[,1])-25, y=max(res_pca$x[,2])-12, 
             label=paste("R2(adonis)=",round(tmp_adonis_result_dis$R2[1],3)))+
    labs(title = paste(z,"Indiviuals pathway PCA",sep = " "))
  ggsave(out1,p1,width = 6,height = 4)
  
  row.names(tmp_site) <- tmp_site$sample
  tmp_apcoa_res <- aPCoA(tmp_dist ~ age,tmp_site,maincov = group)
  tmp_apcoa_matrix <- data.frame(tmp_apcoa_res$plotMatrix)
  identical(row.names(tmp_apcoa_matrix),tmp_site$sample)
  tmp_apcoa_matrix$group <- tmp_site$group
  p2 <- ggplot(data=tmp_apcoa_matrix,aes(x=X1,y=X2,
                                         color=group,shape=group))+
    geom_point(size=2.5)+
    theme_bw()+
    theme(panel.grid = element_blank(),legend.title = element_blank())+
    geom_vline(xintercept = 0,lty="dashed")+
    geom_hline(yintercept = 0,lty="dashed")+
    labs(x=c("PCoA1"),
         y=c("PCoA2"),
         title = paste(z,"Indiviuals pathway PCoA",sep = " ")
    )+
    stat_ellipse(data=tmp_apcoa_matrix,
                 geom = "polygon",
                 aes(fill=group),
                 alpha=0.1)+
    scale_fill_manual(values = color0)+
    scale_color_manual(values = color0)
  ggsave(out2,p2,width = 6,height = 4)
}

library(Maaslin2)
library(factoextra)
library(vegan)
library(microeco)
library(file2meco)
library(magrittr)
library(aplot)
library(compositions)
dir.create("output/pathway",showWarnings = F)
out_dir <- "output/pathway"

load("workfile/z.phyloseq.Rdata")

data_pathway <- read_csv("input/humann.pathabundance.metacyc.ab.csv")
tmp_data_pathway <- data_pathway[grep("\\|",data_pathway$Pathway,invert = T),]
tmp_data_pathway <- as.data.frame(tmp_data_pathway)
row.names(tmp_data_pathway) <- tmp_data_pathway$Pathway
tmp_data_pathway <- tmp_data_pathway %>% select (! Pathway)

tmp_data_pathway_s <- data_pathway[grep("s__",data_pathway$Pathway),]
tmp_data_pathway_s <- as.data.frame(tmp_data_pathway_s)
row.names(tmp_data_pathway_s) <- tmp_data_pathway_s$Pathway
tmp_data_pathway_s <- tmp_data_pathway_s %>% select (! Pathway)

tmp_data_run1 <- as.data.frame(t(tmp_data_pathway))
tmp_data_run2 <- as.data.frame(t(tmp_data_pathway_s))

tmp_meta <- metadatadf %>%
  #filter(group %in% c("health","infected"))
  filter(groups %in% group[c(1,2)])
tmp_sam <- tmp_meta$sample
tmp_groups <- factor(tmp_meta$groups,level=c(group[1],group[2]))
tmp_meta$groups <- factor(tmp_meta$groups,level=c(group[1],group[2]))
tmp_meta_b <- tmp_meta
#tmp_meta1 <- tmp_meta %>% filter(enterotype == "ET1")
#tmp_meta2 <- tmp_meta %>% filter(enterotype == "ET2")

Maaslin2(
  input_data = tmp_data_run1, 
  input_metadata = tmp_meta, 
  min_prevalence = 0.1,
  normalization = "NONE",
  transform = "LOG",
  analysis_method = "LM",
  correction = "BH",
  output = (paste(out_dir,"cross.pathway_all_metacyc",sep = '/')), 
  fixed_effects = c("groups","age","sex","bmi"),
  reference = c("groups,HC"),
)

Maaslin2(
  input_data = tmp_data_run2, 
  input_metadata = tmp_meta, 
  min_prevalence = 0.1,
  normalization = "NONE",
  transform = "LOG",
  analysis_method = "LM",
  correction = "BH",
  output = paste(out_dir,"cross.pathway_all_metacyc_s",sep = '/'), 
  fixed_effects = c("groups"),
  reference = c("groups,HC"),
)
tmp_maaslin <- read_table(paste(out_dir,"cross.pathway_all_metacyc/significant_results.tsv",sep = '/'))
tmp_maaslin <- tmp_maaslin %>% filter(metadata == 'groups') %>% arrange( desc(coef))
tmp_maaslin$feature <- gsub("\\.\\."," ",tmp_maaslin$feature)
tmp_maaslin$feature <- gsub("\\."," ",tmp_maaslin$feature)
tmp_maaslin$feature <- factor(tmp_maaslin$feature,levels = tmp_maaslin$feature)
tmp_maaslin$log_qval <- -log10(tmp_maaslin$qval)

ggplot(tmp_maaslin)+
  geom_bar(aes(x = feature, y = coef, fill = log_qval),stat = "identity")+
  geom_errorbar(aes(x = feature, ymin = coef - stderr,ymax = coef + stderr),width = 0.2,stat = "identity")+
  theme_bw()+
  theme(
    panel.grid = element_blank()
  )+
  scale_fill_gradient2(low = "grey70",high = "grey50")+
  coord_flip()
ggsave(paste(out_dir,"cross.pathway_diff.pdf",sep = '/'),width = 8,height = 10)

tmp_masslin_up_top10 <- tmp_maaslin %>% arrange(desc(log_qval)) %>% filter(coef > 0) %>% head(10)
tmp_masslin_down_top10 <- tmp_maaslin %>% arrange(desc(log_qval)) %>% filter(coef < 0) %>% head(10)
tmp_masslin_top20 <- bind_rows(tmp_masslin_down_top10,tmp_masslin_up_top10)
tmp_masslin_top20 <- tmp_masslin_top20 %>% arrange(desc(coef))
tmp_masslin_top20$feature <- factor(tmp_masslin_top20$feature,levels = tmp_masslin_top20$feature)
ggplot(tmp_masslin_top20)+
  geom_bar(aes(x = feature, y = coef, fill = log_qval),stat = "identity")+
  geom_errorbar(aes(x = feature, ymin = coef - stderr,ymax = coef + stderr),width = 0.2,stat = "identity")+
  theme_bw()+
  theme(
    panel.grid = element_blank()
  )+
  scale_fill_gradient2(low = "grey70",high = "grey50")+
  coord_flip()
ggsave(paste(out_dir,"cross.pathway_diff_top20.pdf",sep = '/'),width = 8,height = 6)

meco_pathway_metacyc <- humann2meco(in_pathway_metacyc, 
                                    db = "MetaCyc", 
                                    sample_table = tmp_meta 
)
meco_pathway_metacyc$tidy_dataset()
meco_pathway_metacyc$cal_abund(rel = T)
meco_pathway_metacyc$taxa_abund$Superclass1 %<>% .[!grepl("unclass", rownames(.)), ]
meco_pathway_metacyc$taxa_abund$Superclass2 %<>% .[!grepl("unclass", rownames(.)), ]
meco_pathway_metacyc$taxa_abund$pathway %<>% .[!grepl("unclass", rownames(.)), ]
meco_pathway_metacyc$taxa_abund$Species %<>% .[!grepl("unclass", rownames(.)), ]
meco_pathway_metacyc$taxa_abund$Genus %<>% .[!grepl("unclass", rownames(.)), ]

meco_pathway_metacyc_used <- trans_abund$new(meco_pathway_metacyc, taxrank = "Superclass2", 
                                             ntaxa = 10, use_percentage = FALSE)
meco_pathway_metacyc_used$ylabname <- "Relative Abundace"
meco_pathway_metacyc_used$plot_bar(facet = "groups", bar_type = "notfull")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))

ggsave(paste(out_dir,"cross.pathway_metacyc_bar.all.pdf",sep = '/'),width = 20,height = 6)
tmp_plot_data <- meco_pathway_metacyc_used$data_abund %>% arrange(desc(all_mean_abund))
tmp_plot_data$Taxonomy <- factor(tmp_plot_data$Taxonomy,levels = rev(unique(tmp_plot_data$Taxonomy)))
tmp_plot_pathway <- unique(tmp_plot_data$Taxonomy)[1:15]

ggplot(tmp_plot_data %>% filter(Taxonomy %in% tmp_plot_pathway),
       aes(x=groups,y=Abundance,fill = Taxonomy)) +
  stat_summary(fun=median, geom="bar" ,width = 0.6,position = "stack")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(
    values = rev(c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02",
                   "#A6761D","#666666","#8DCEBB","#ECAF80","#BAB7D9","#F394C4",
                   "#B2D28E","#F2D480","#D2BA8E"))
  )+
  scale_y_continuous(breaks = seq(0,1,0.1),limits = c(0,1))
ggsave((paste(out_dir,"cross.pathway_metacyc_bar.mean.pdf",sep = '/')),width = 8,height = 6)

tmp_data_stat <- data.frame()
for(i in 1:length(tmp_plot_pathway)){
  tmp_data <- tmp_plot_data %>% filter(Taxonomy == tmp_plot_pathway[i])
  tmp <- summary_stat(tmp_data,3,11) %>% mutate(Taxaonomy = tmp_plot_pathway[i]) 
  tmp_data_stat <- rbind(tmp_data_stat,tmp)
}
tmp_data_stat_short <- tmp_data_stat %>% 
  select(Group,Mean,Taxaonomy) %>% 
  spread(Group,Mean) %>% 
  mutate(Var = TA - HC)
tmp_data_stat_short$Taxaonomy <-  factor(tmp_data_stat_short$Taxaonomy,
                                         levels = rev((tmp_data_stat_short$Taxaonomy)))

ggplot(tmp_data_stat_short)+
  geom_bar(aes(x=Taxaonomy,y = Var,fill = Taxaonomy),
           stat = "identity",position = "identity")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(
    values = c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02",
               "#A6761D","#666666","#8DCEBB","#ECAF80","#BAB7D9","#F394C4",
               "#B2D28E","#F2D480","#D2BA8E")
  )+
  coord_flip()+
  labs(y="",x="CPM variation")
ggsave(paste(out_dir,"cross.pathway_metacyc_bar.var.pdf",sep = '/'),width = 11,height = 6)

meco_pathway_metacyc_diff <- trans_diff$new(meco_pathway_metacyc,
                                            method =  "lefse",group = "groups",
                                           taxa_level = "pathway",
                                           alpha = 0.05,p_adjust_method = "fdr"
                                           #reference = c("Type,Health"),fixed_effects = c("Type","Sex")
                                           )
meco_pathway_metacyc_diff$plot_diff_bar(use_number = 1:60)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_color_manual(
    values = color0[c(1,2)]
  )+
  scale_fill_manual(
    values = color0[c(1,2)]
  )
  
ggsave(paste(out_dir,"cross.pathway_metacyc_bar.lefse.pdf",sep = '/'),width = 10,height = 15)
meco_pathway_metacyc$cal_betadiv(method = "bray")
meco_pathway_metacyc_distance <- as.dist(meco_pathway_metacyc$beta_diversity$bray)
meco_pathway_metacyc_beta <- trans_beta$new(dataset = meco_pathway_metacyc, 
                                            group = "groups", measure = "bray")

meco_pathway_metacyc_beta$cal_ordination(method = "PCoA")
meco_pathway_metacyc_beta_adonis <- adonis2(meco_pathway_metacyc_distance~groups+age+bmi+sex,tmp_meta)

meco_pathway_metacyc_beta_tmp <- meco_pathway_metacyc_beta$res_ordination$scores
meco_pathway_metacyc_beta_tmp2 <- trans_env$new(dataset = meco_pathway_metacyc, 
                                                add_data = meco_pathway_metacyc_beta_tmp[, 1:2])
# 'KW_dunn' for non-parametric test
tmp_list = combine_list(meco_pathway_metacyc_beta_tmp$groups)
a <- ggplot()+
  geom_point(data = meco_pathway_metacyc_beta_tmp,
             mapping = aes(x=PCo1,y=PCo2,color=groups,shape=groups),
             size=2.5,alpha = 0.7)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.88,0.88),
        axis.text.y = element_text(angle = 90,hjust = 0.5,vjust = 1)
  )+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  labs(x=c("PCoA1"),y=c("PCoA2"))+
  stat_ellipse(data=meco_pathway_metacyc_beta_tmp,geom = "polygon",
               aes(x=PCo1,y=PCo2,fill=groups),
               color =NA,alpha=0.1)+
  scale_fill_manual(values = color0)+
  scale_color_manual(values = color0)+
  scale_x_continuous(limits = c(-0.7,0.4))+
  scale_y_continuous(limits = c(-0.3,0.3))

b <- ggplot(data=meco_pathway_metacyc_beta_tmp,aes(x=groups,y=PCo1,fill=groups))+
  #geom_violin(aes(x=tmp,y=X1,fill=group),width = 1,trim = F,color = "white")+
  geom_boxplot(width = 0.8,alpha = 0.7)+
  theme_bw()+
  theme(
    legend.position = "NONE",
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )+
  labs(x=element_blank(),y=element_blank())+
  coord_flip()+
  scale_fill_manual(values = color0)+
  scale_y_continuous(limits = c(-0.7,0.4))+
  stat_compare_means(
    comparisons = tmp_list,
    method = "wilcox.test", label = "p.format")

c <- ggplot(data=meco_pathway_metacyc_beta_tmp,aes(x=groups,y=PCo2,fill=groups))+
  #geom_violin(aes(x=tmp,y=X1,fill=group),width = 1,trim = F,color = "white")+
  geom_boxplot(width = 0.8,alpha = 0.7)+
  theme_bw()+
  theme(
    legend.position = "NONE",
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )+
  labs(x=element_blank(),y=element_blank())+
  scale_fill_manual(values = color0)+
  scale_y_continuous(limits = c(-0.3,0.3))+
  stat_compare_means(
    comparisons = tmp_list,method = "wilcox.test", label = "p.format"
  )

meco_pathway_metacyc_beta_adonis$item <- row.names(meco_pathway_metacyc_beta_adonis)
meco_pathway_metacyc_beta_adonis_p <- meco_pathway_metacyc_beta_adonis[!is.na(meco_pathway_metacyc_beta_adonis$F),]
meco_pathway_metacyc_beta_adonis_p <- as.data.frame(t(meco_pathway_metacyc_beta_adonis_p))
meco_pathway_metacyc_beta_adonis_p <- meco_pathway_metacyc_beta_adonis_p[3:5,]
meco_pathway_metacyc_beta_adonis_p$posy <- c(0.25,0.5,0.75)
meco_pathway_metacyc_beta_adonis_p$posx <- 1
meco_pathway_metacyc_beta_adonis_p$labels <- paste(row.names(meco_pathway_metacyc_beta_adonis_p),round(as.numeric(meco_pathway_metacyc_beta_adonis_p[,1]),3),sep = ':')

d <- ggplot()+
  geom_bar(data = as.data.frame(matrix(c(1),nrow = 1)),aes(x = V1,y = V1),stat ="identity",fill="white")+
  geom_text(data = meco_pathway_metacyc_beta_adonis_p,mapping = aes(x = posx,y = posy,label = labels))+
  theme_bw()+theme(
    legend.position = "NONE",
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
#plot_grid(b,d,a,c,align = "hv",rel_widths = c(1,0.225),rel_heights = c(0.225,1))
aplot::plot_list(b,d,a,c,widths = c(1,0.225),heights = c(0.225,1))
ggsave(paste(out_dir,"cross.pathway_metacyc.pcoa.pdf",sep = '/'),width = 5,height = 5)

meco_pathway_metacyc_beta$plot_ordination(plot_color = "groups", 
                                                plot_shape = "groups", 
                                                plot_type = c("point", "ellipse"))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_color_manual(
    values = color0[c(1,2)]
  )+
  scale_fill_manual(
    values = color0[c(1,2)]
  )

tmp_meta <- metadatadf %>% filter(groups %in% group[c(1,2)]) %>% select(sample,groups)
tmp_otu <- as.data.frame(otu_table(phyloseqin))
tmp_otu_g <- as.data.frame(otu_table(phyloseqin_all))
tmp_otu_g <- tmp_otu_g[grep("g__",row.names(tmp_otu_g),perl = TRUE),]
tmp_otu_g <- tmp_otu_g[,colnames(tmp_otu_g) %in% tmp_meta$sample,]
tmp_data_pathway <- tmp_data_pathway[,tmp_meta$sample]

tmp_res <- data.frame(pathway="",cor = 0, p = 0, prevelance = "")
n = 1
for(i in 1:nrow(tmp_data_pathway)){
  tmp_d1 <- as.data.frame(t(tmp_data_pathway[i,]))
  tmp_d2 <- as.data.frame(t(tmp_otu_g)[,"g__Bifidobacterium"])
  tmp_p <- colnames(tmp_d1)
  tmp_d1$sample <- row.names(tmp_d1)
  tmp_d2$sample <- row.names(tmp_d2)
  tmp_d1 <- inner_join(tmp_d1,tmp_d2,by = "sample")
  colnames(tmp_d1) <- c("pathway","sample","sp")

  #tmp_d1$sample <- row.names(tmp_d1)
  #tmp_d <- inner_join(tmp_d1,tmp_meta,by = "sample")
  tmp_rate <- length((which(tmp_d1$pathway != 0)))/nrow(tmp_d1)
  if(tmp_rate > 0.2){
    tmp_cor <- cor.test(tmp_d1$pathway,tmp_d1$sp,method = "spearman")
    tmp_res[n,] <- c(tmp_p,tmp_cor$estimate,tmp_cor$p.value,tmp_rate)
    n = n + 1
  }
}
tmp_res[,2:4] <- apply(tmp_res[,2:4],2,as.numeric)
tmp_res$padj <- p.adjust(tmp_res$p,method = "fdr")
tmp_res_sig <- tmp_res %>% filter(padj < 0.05)
write_csv(tmp_res,paste(out_dir,"pathway2Bifido.cov.csv",sep = '/'))

tmp_d1 <- as.data.frame( t(tmp_data_pathway["PWY-6527: stachyose degradation",]))
colnames(tmp_d1) <- "pathway"
tmp_d1$sample <- row.names(tmp_d1)
tmp_bfido <- c("s__Bifidobacterium_adolescentis","s__Bifidobacterium_bifidum","s__Bifidobacterium_dentium",
               "s__Bifidobacterium_longum","s__Bifidobacterium_pseudocatenulatum")
tmp_db <- as.data.frame(t(tmp_otu))[,tmp_bfido]
tmp_db$sample <- row.names(tmp_db)
tmp_db <- inner_join(tmp_db,tmp_d1,by ="sample")
tmp_db <- inner_join(tmp_db,tmp_meta,by = "sample")
plot_l <- list()
for(i in 1:length(tmp_bfido)){
  tmp_d <- tmp_db[,c(tmp_bfido[i],"pathway","groups")]
  colnames(tmp_d) <- c("sp","pathway","groups")
  tmp_d[tmp_d == 0] <- 0.000001
  tmp_d$pathway <- as.numeric(clr(tmp_d$pathway))
  tmp_d$sp <- as.numeric(clr(tmp_d$sp))
  p <- ggplot(tmp_d,aes(y=pathway,x=sp)) + 
    stat_smooth(method='lm',formula = y~x,color = "grey60")+
    geom_point(aes(color = groups),size=3,alpha = 0.7)+
    stat_cor(cor.coef.name = "rho",method = "spearman",
             label.y.npc = 1,label.x.npc = 0.35)+
    #stat_poly_eq(aes(label = ..eq.label..),formula = tmp_formula,label.y.npc = 0.95,label.x.npc=0.87, parse = TRUE,)+
    theme_bw()+
    theme(panel.grid = element_blank())+
    scale_color_manual(
      values = color0
    )+
    labs(title = tmp_bfido[i], x ="CLR count",y = "stachyose degradation")
  plot_l[[i]] <- p
}
plot_grid(plot_l[[1]],plot_l[[2]],plot_l[[3]],plot_l[[4]],plot_l[[5]],nrow = 2)
ggsave(paste(out_dir,"stachy_pathway2Bifdo.cov.pdf",sep = '/'),width = 11,height= 6)

tmp_maaslin_s <- read_table(paste(out_dir,"cross.pathway_all_metacyc_s/significant_results.tsv",sep = '/'))
tmp_maaslin_s <- tmp_maaslin_s[grepl("PWY.6527",tmp_maaslin_s$feature),]
tmp_maaslin_s <- tmp_maaslin_s %>% filter(metadata == 'groups') %>% arrange( desc(coef))
tmp_maaslin_s$feature <- gsub("\\.\\."," ",tmp_maaslin_s$feature)
tmp_maaslin_s$feature <- gsub("\\."," ",tmp_maaslin_s$feature)
tmp_maaslin_s$feature <- gsub(".*s__","",tmp_maaslin_s$feature)
tmp_maaslin_s$feature <- factor(tmp_maaslin_s$feature,levels = tmp_maaslin_s$feature)
tmp_maaslin_s$log_qval <- -log10(tmp_maaslin_s$qval)

ggplot(tmp_maaslin_s)+
  geom_bar(aes(x = feature, y = coef, fill = log_qval),stat = "identity")+
  geom_errorbar(aes(x = feature, ymin = coef - stderr,ymax = coef + stderr),width = 0.2,stat = "identity")+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1)
  )+
  scale_fill_gradientn(colors = colorRampPalette(brewer.pal(5,"BuPu"))(50))+
  labs(title = "PWY-6527: stachyose degradation")
ggsave(paste(out_dir,"cross.pathway_diff_stachy.pdf",sep = '/'),width = 6,height = 5)

tmp_d1$sample <- row.names(tmp_d1)
tmp_d1 <- inner_join(tmp_d1,tmp_d2,by = "sample")
colnames(tmp_d1) <- c("pathway","sample","sp")
tmp_d1 <- inner_join(tmp_d1,tmp_meta,by = "sample")
tmp_d1 <- tmp_d1 %>% filter(sp < 40)
tmp_d1[tmp_d1 == 0] <- 0.000001
tmp_d1$pathway <- as.numeric(clr(tmp_d1$pathway))
tmp_d1$sp <- as.numeric(clr(tmp_d1$sp))

ggplot(tmp_d1,aes(y=pathway,x=sp)) + 
  stat_smooth(method='lm',formula = y~x,color = "grey60")+
  geom_point(aes(color = groups),size=3,alpha = 0.7)+
  stat_cor(cor.coef.name = "rho",method = "spearman",
           label.y.npc = 1,label.x.npc = 0.35)+
  #stat_poly_eq(aes(label = ..eq.label..),formula = tmp_formula,label.y.npc = 0.95,label.x.npc=0.87, parse = TRUE,)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_color_manual(
    values = color0
  )+
  labs(x = "Bifidobacterium", y = "PWY-6527: stachyose degradation")
ggsave(paste(out_dir,"stachy_pathway2Bifdo.genus.cov.pdf",sep = '/'),width = 4.5,height= 4)

tmp_sp <- read.table("input/stachyose_degradation.txt")
tmp_sp_b <- tmp_sp
tmp_sp <- tmp_sp[,1]
tmp_sp_stat <- data.frame()
tmp_pre_cut = 0.2
i = 1
for(i in 1:length(tmp_sp)){
  if(! tmp_sp[i] %in% row.names(tmp_otu)){
    print(paste(c("sp is not exist:",tmp_sp[i]),collapse = " "))
    next
  }
  tmp_d <- as.data.frame(t(tmp_otu[which(row.names(tmp_otu)==tmp_sp[i]),tmp_meta$sample]))
  
  tmp_d$sample <- row.names(tmp_d)
  colnames(tmp_d) <- c("value","sample")
  if(length(which(tmp_d$value != 0)) < nrow(tmp_d) * tmp_pre_cut ){
    print(paste(c("sp have lower prevalence than threshold:",tmp_sp[i]),collapse = " "))
    next
  }
  tmp_d_s <- inner_join(tmp_meta,tmp_d,by = "sample")
  tmp_d_s$tmp <- tmp_sp[i]
  tmp_stat <- summary_stat(tmp_d_s,3,4)
  tmp_sp_stat <- bind_rows(tmp_sp_stat,tmp_stat)
}
tmp_sp_stat %<>% arrange(desc(Median))
tmp_sp <- tmp_sp_stat$Group

tmp_res <- plot_sp_box(tmp_otu,tmp_meta,tmp_sp)
tmp_res_p <- tmp_res[[1]]
plot_grid(tmp_res_p[[1]],tmp_res_p[[2]],tmp_res_p[[3]],tmp_res_p[[4]],
          tmp_res_p[[5]],tmp_res_p[[6]],tmp_res_p[[7]],tmp_res_p[[8]],
          tmp_res_p[[9]],tmp_res_p[[10]],tmp_res_p[[11]],tmp_res_p[[12]],
          tmp_res_p[[13]],tmp_res_p[[14]],tmp_res_p[[15]],tmp_res_p[[16]],
          tmp_res_p[[17]],tmp_res_p[[18]],tmp_res_p[[19]],tmp_res_p[[20]],
          tmp_res_p[[21]],tmp_res_p[[22]],tmp_res_p[[23]],tmp_res_p[[24]],
          tmp_res_p[[25]],tmp_res_p[[26]],tmp_res_p[[27]],nrow = 4)
ggsave(paste(out_dir,"stachy_pathway.sp.pdf",sep = '/'),width = 25,height = 16)

#tmp_sp_metaphlan <- read.table("input/stachyose.degradation_sp_metaphlan.txt")[,1]
tmp_sp_s <- tmp_sp_b %>% filter (V2 != "non")
tmp_sp_s <- tmp_sp_s[,1]
tmp_com <- t(combn(tmp_sp_s,2))
tmp_res <- data.frame(sp1="",sp2="",cor_pearson = 0, p_pearson = 0, cor_spearman = 0, p_spearman = 0,prevelance = "")
n = 1
for(i in 1:nrow(tmp_com)){
  if(! tmp_com[i,1] %in% row.names(tmp_otu)){
    print(paste(c("sp is not exist:",tmp_com[i,1]),collapse = " "))
    next
  }
  if(! tmp_com[i,2] %in% row.names(tmp_otu)){
    print(paste(c("sp is not exist:",tmp_com[i,2]),collapse = " "))
    next
  }
  tmp_data1 <- t(tmp_otu[tmp_com[i,1],tmp_meta$sample])
  tmp_data2 <- t(tmp_otu[tmp_com[i,2],tmp_meta$sample])
  
  tmp_n1 <- colnames(tmp_data1)
  tmp_n2 <- colnames(tmp_data2)
  tmp_cor1 <- cor.test(tmp_data1[,1],tmp_data2[,1])
  tmp_cor2 <- cor.test(tmp_data1[,1],tmp_data2[,1],method = "spearman")
  tmp_zero2 <- length(which(tmp_data2[,1]!=0))
  tmp_tn <- length(tmp_data1[,1])
  tmp_zeroall <- paste(tmp_zero2,tmp_tn,sep = "|")
  tmp_res[n,] <- c(tmp_n1,tmp_n2,tmp_cor1$estimate,tmp_cor1$p.value,tmp_cor2$estimate,tmp_cor2$p.value,tmp_zeroall)
  n = n + 1
  print(i)
}
tmp_res[,3:6] <- apply(tmp_res[,3:6],2,as.numeric)
write_csv(tmp_res,paste(out_dir,"stachyose_pathway.sigsp_cov.csv",sep = '/'))

tmp_res_sig <- tmp_res %>% filter(sp1 != sp2,p_spearman < 0.05)
tmp_res_sig <- tmp_res_sig %>% filter(sp1 == "s__Bifidobacterium_longum"|sp2 == "s__Bifidobacterium_longum")

tmp_sp_sig <- unique(c(tmp_res_sig$sp1,tmp_res_sig$sp2))
tmp_sp_sig <- tmp_sp_sig[tmp_sp_sig != "s__Bifidobacterium_longum"]

plot_l <- list()
for(i in 1:length(tmp_sp)){
  tmp_d <- tmp_otu[c("s__Bifidobacterium_longum",tmp_sp[i]),tmp_meta$sample]
  tmp_d[tmp_d == 0] <- 0.000001
  tmp_d <-as.data.frame(clr(as.matrix((tmp_d))))
  tmp_d <- as.data.frame(t(tmp_d))
  colnames(tmp_d) <- c("sp1","sp2")
  tmp_d$sample <- row.names(tmp_d)
  tmp_d <- inner_join(tmp_d,tmp_meta,by = "sample")
  tmp_formula = y ~ x
  p <- ggplot(tmp_d,aes(x=sp1,y=sp2)) + 
    stat_smooth(method='lm',formula = tmp_formula,color = "grey60")+
    geom_point(aes(color = groups),size=3,alpha = 0.7)+
    stat_cor(cor.coef.name = "rho",method = "spearman",
             label.y.npc = 1,label.x.npc = 0.35)+
    #stat_poly_eq(aes(label = ..eq.label..),formula = tmp_formula,label.y.npc = 0.95,label.x.npc=0.87, parse = TRUE,)+
    theme_bw()+
    theme(panel.grid = element_blank())+
    scale_color_manual(
      values = color0
    )+
    labs(x = "s__Bifidobacterium_longum", y = tmp_sp[i])
  
  plot_l[[i]] <- p
}
plot_grid(plot_l[[1]],plot_l[[2]],plot_l[[3]],plot_l[[4]],plot_l[[5]],plot_l[[6]],plot_l[[7]],plot_l[[8]],plot_l[[9]],
          plot_l[[10]],plot_l[[11]],plot_l[[12]],plot_l[[13]],plot_l[[14]],plot_l[[15]],plot_l[[16]],plot_l[[17]],plot_l[[18]],
          plot_l[[19]],plot_l[[20]],plot_l[[21]],plot_l[[22]],plot_l[[23]],plot_l[[24]],plot_l[[25]],plot_l[[26]],plot_l[[27]],
          nrow = 5)
ggsave(paste(out_dir,"stachyose_degradation.sp_cov.pdf",sep = '/'),width = 25,height = 16)

tmp_meta2 <- metadatadf %>%
  #filter(group %in% c("health","infected"))
  filter(group1 %in% groupl[c(2,3)])
tmp_sam <- tmp_meta2$sample
tmp_group1 <- factor(tmp_meta2$group1,level=c(groupl[2],groupl[3]))
tmp_meta2$group1 <- factor(tmp_meta2$group1,level=c(groupl[2],groupl[3]))

fit_type = Maaslin2(
  input_data = tmp_data_run1, 
  input_metadata = tmp_meta2, 
  min_prevalence = 0.1,
  max_significance = 0.5,
  normalization = "NONE",
  transform = "LOG",
  analysis_method = "LM",
  correction = "BH",
  output = paste(out_dir,"long.pathway_all_metacyc",sep = '/'), 
  fixed_effects = c("group1","interval"),
  reference = c("group1,BE"),
  random_effects = c("pair"),
)

fit_type = Maaslin2(
  input_data = tmp_data_run2, 
  input_metadata = tmp_meta2, 
  min_prevalence = 0.1,
  normalization = "NONE",
  transform = "LOG",
  analysis_method = "LM",
  correction = "BH",
  output = paste(out_dir,"long.pathway_all_metacyc_s",sep = '/'), 
  fixed_effects = c("group1","interval"),
  reference = c("group1,BE"),
  random_effects = c("pair"),
)

meco_pathway_metacyc <- humann2meco(in_pathway_metacyc, 
                                    db = "MetaCyc", 
                                    sample_table = tmp_meta2 
)
