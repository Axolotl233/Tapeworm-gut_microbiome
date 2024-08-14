library(Maaslin2)
library(patchwork)
load("workfile/z.phyloseq.Rdata")

dir.create("output/grid",showWarnings = F)
out_dir <- "output/grid"
data_grid <- read.table("input/z.merge.result.GRiD.txt",header = T,row.names = 1)
data_grid_bak <- data_grid

phyloseqin_used <- subset_samples(phyloseqin, groups %in% c(group[1],group[2]))
tmp_otu <- as.data.frame(t(otu_table(phyloseqin_used)))
identical(row.names(tmp_otu),data_grid$sample)
tmp_meta <- metadatadf %>% dplyr::filter(groups %in% group[c(1,2)])
tmp_meta_cross <- metadatadf %>% dplyr::filter(groups %in% group[c(1,2)]) %>%  dplyr::select(sample,groups)
tmp_list <- combine_list(tmp_meta_cross$groups)


#data_grid_cov <- read.table("input/z.merge.result.GRiD_cov.txt",header = T,row.names = 1)

data_grid <- as.data.frame(t(data_grid)) %>% filter(row.names(.) %in% tmp_meta_cross$sample)
data_grid$sample <- row.names(data_grid)
data_grid <- left_join(tmp_meta_cross,data_grid,by = "sample")
res_cross = data.frame(sp = "",gr1_median="",gr2_median="",p = 0)
n = 1
for(i in 3:ncol(data_grid)){
  tmp_data <- data_grid[,c(1,2,i)]
  tmp_n <- colnames(data_grid)[i]
  colnames(tmp_data) <- c("sample","groups","sp")
  if(length(which(is.na(tmp_data[,3]))) > 0.5*(length(tmp_data[,3]))){
    next
  }else{
    tmp_data <- tmp_data[!is.na(tmp_data[,3]),]
    tmp1 <- tmp_data %>% group_by(groups) %>% summarise(med = median(sp))
    tmp2 <- wilcox.test(sp ~ groups,data = tmp_data)
    res_cross[n,] <- c(tmp_n,tmp1$med[1],tmp1$med[2],tmp2$p.value)
    n = n + 1
  }
}
res_cross[,c(2,3,4)] <- apply(res_cross[,c(2,3,4)],2,as.numeric)
res_cross$padj <- p.adjust(res_cross$p,method = "BH")
res_cross <- res_cross[!is.na(res_cross$padj),]
res_cross <- res_cross[order(res_cross$padj),]
write_csv(x = res_cross,paste(out_dir,"cross.grid.csv",sep = '/'))

data_grid_masslin2 <- as.data.frame(t(data_grid_bak))
data_grid_masslin2 <- apply(data_grid_masslin2, c(1,2), function(x){ifelse(is.na(x),0,x)})
Maaslin2(
  input_data = data_grid_masslin2, 
  input_metadata = tmp_meta, 
  min_prevalence = 0.1,
  normalization = "NONE",
  transform = "LOG",
  analysis_method = "LM" ,
  output = paste(out_dir,"cross.grid",sep = '/'),
  correction = "BH", 
  fixed_effects = c("groups","sex","age","bmi"),
  reference = c("groups,HC"),
)

sp_specific_g <- c("Bifidobacterium_pseudocatenulatum",
                 "Bifidobacterium_longum",
                 "Bifidobacterium_adolescentis")
tmp_data <- data_grid[,c(1,2,which(colnames(data_grid) %in% sp_specific_g))]
tmp_data_l <- tmp_data %>% gather(key = "species",value = "value",sp_specific_g)
tmp_data_l$value <- ifelse(is.na(tmp_data_l$value),0,tmp_data_l$value)
ggplot(tmp_data_l,aes(x=species,y=value))+
  geom_boxplot(aes(color=groups),outlier.shape = NA)+
  geom_jitter(aes(color=groups,shape=groups),position = position_jitterdodge(dodge.width = 0.75,jitter.width = 0.15, jitter.height = 0),
              size = 2,alpha = 0.7)+
  theme_cowplot()+
  theme(axis.title.x=element_blank(), 
#        axis.text.x = element_text(angle = 90,vjust = -0.5),
        panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
  scale_color_manual(
    values = color0
  )+
  labs(y="Observed value")
ggsave(paste(out_dir,"grid.bifido.pdf",sep = '/'),width = 4,height = 5)

for(i in 1:length(sp_specific)){
  tmp_data <- data_grid[,c(1,2,which(colnames(data_grid) == sp_specific_g[i]))]
  colnames(tmp_data) <- c("sample","groups","value")
  p <- ggplot(tmp_data,aes(x=groups,y=value))+
    geom_boxplot(aes(color=groups),outlier.shape = NA)+
    geom_jitter(aes(color=groups,shape=groups),size = 2,height=0,width = 0.3)+
    theme_cowplot()+
    theme(axis.title.x=element_blank(), 
          legend.position="none",
          panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
    scale_color_manual(
      values = color0
    )+
    labs(y="Observed value",title = tmp_n)+
    stat_compare_means(
      comparisons = tmp_list,
      method = "wilcox.test", label = "p.format")
  plot_l[[i]] <- p
}

plot_grid(plot_l[[1]],plot_l[[2]],plot_l[[3]],nrow = 1)
ggsave(paste(out_dir,"cross.stachyose_degradation.grid.ab.pdf",sep = '/'),width = 25,height = 25)

sp_specific <- read.csv("input/stachyose_degradation_sp.csv",head = T)
sp_specific <- sp_specific[!is.na(sp_specific$grid),]
tmp_otu <- as.data.frame(t(otu_table(phyloseqin_used)))
identical(row.names(tmp_otu),data_grid$sample)

plot_l1 <- list()
plot_l2 <- list()
i = 1
for(i in 1:nrow(sp_specific)){
  tmp_n = sp_specific[i,2]
  tmp_data <- data_grid[,c(1,2,which(colnames(data_grid) == sp_specific[i,1]))]
  colnames(tmp_data) <- c("sample","groups","grid")
  tmp_data$ab <- tmp_otu[,tmp_n]
  tmp_data <- tmp_data[!is.na(tmp_data[,3]),]
  tmp_data_l <- tmp_data %>% gather(key = "class",value = "value",grid,ab)
  p1 <- ggplot(tmp_data_l,aes(x=groups,y=value))+
    geom_boxplot(aes(color=groups),outlier.shape = NA)+
    geom_jitter(aes(color=groups,shape=groups),size = 2,height=0,width = 0.3)+
    theme_cowplot()+
    theme(axis.title.x=element_blank(), 
          legend.position="none",
          panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
    scale_color_manual(
      values = color0
    )+
    labs(y="Observed value",title = tmp_n)+
    stat_compare_means(
      comparisons = tmp_list,
      method = "wilcox.test", label = "p.format")+
    facet_wrap(.~class,scales = "free")
  plot_l1[[i]] <- p1
  tmp_formula <- y ~ x
  p2 <- ggplot(tmp_data,aes(x=ab,y=grid)) + 
    stat_smooth(method='lm',formula = tmp_formula,color = "grey70")+
    geom_point(aes(color=groups),size=4)+
    scale_color_manual(
      values = c(color0)
    )+
    stat_cor(cor.coef.name = "rho",method = "spearman",
             label.y.npc = 1,label.x.npc = 0.35)+
    #stat_poly_eq(aes(label = ..eq.label..),formula = tmp_formula,label.y.npc = 0.95,label.x.npc=0.87, parse = TRUE,)+
    theme_light()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position="none"
    )+
    labs(title = tmp_n)
  plot_l2[[i]] <- p2
}

plot_grid(plot_l1[[1]],plot_l1[[2]],plot_l1[[3]],plot_l1[[4]],plot_l1[[5]],plot_l1[[6]],
          plot_l1[[7]],plot_l1[[8]],plot_l1[[9]],plot_l1[[10]],plot_l1[[11]],plot_l1[[12]],
          plot_l1[[13]],plot_l1[[14]],plot_l1[[15]],plot_l1[[16]],plot_l1[[17]],plot_l1[[18]],
          plot_l1[[19]],plot_l1[[20]],plot_l1[[21]],plot_l1[[22]],plot_l1[[23]],plot_l1[[24]],
          nrow = 6)
ggsave(paste(out_dir,"cross.stachyose_degradation.grid.ab.pdf",sep = '/'),width = 25,height = 25)
plot_grid(plot_l2[[1]],plot_l2[[2]],plot_l2[[3]],plot_l2[[4]],plot_l2[[5]],plot_l2[[6]],
          plot_l2[[7]],plot_l2[[8]],plot_l2[[9]],plot_l2[[10]],plot_l2[[11]],plot_l2[[12]],
          plot_l2[[13]],plot_l2[[14]],plot_l2[[15]],plot_l2[[16]],plot_l2[[17]],plot_l2[[18]],
          plot_l2[[19]],plot_l2[[20]],plot_l2[[21]],plot_l2[[22]],plot_l2[[23]],plot_l2[[24]],
          nrow = 4)
ggsave(paste(out_dir,"cross.stachyose_degradation.grid.ab.cor.pdf",sep = '/'),width = 25,height = 25)

tmp_meta_long <- metadatadf %>% 
  dplyr::filter(group1 %in% groupl[c(2,3)]) %>% dplyr::select(sample,pair,group1)
data_grid <- as.data.frame(t(data_grid_bak)) %>% filter(row.names(.) %in% tmp_meta_long$sample)
data_grid$sample <- row.names(data_grid)
data_grid <- left_join(tmp_meta_long,data_grid,by = "sample")

tmp_data <- data_grid[,c(1,2,3,which(colnames(data_grid) %in% sp_specific_g))]
tmp_data_l <- tmp_data %>% gather(key = "species",value = "value",sp_specific_g)
tmp_data_l$group1 <- factor(tmp_data_l$group1,levels = groupl[c(2,3)])
tmp_data_l$value <- ifelse(is.na(tmp_data_l$value),0,tmp_data_l$value)
tmp_list = combine_list(tmp_data_l$group1)
ggplot(tmp_data_l,aes(x=group1,y=value,color=group1))+
  geom_boxplot(aes(color=group1),outlier.shape = NA)+
  geom_point(size = 2,alpha = 0.7)+
  geom_line(aes(group = pair),color = "grey80",linetype = "dashed")+
  theme_cowplot()+
  facet_wrap(.~species)+
  theme(axis.title.x=element_blank(), 
        #        axis.text.x = element_text(angle = 90,vjust = -0.5),
        panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
  scale_color_manual(
    values = color0[c(2,3)]
  )+
  labs(y="Observed value")+
  stat_compare_means(method="wilcox.test",
                     comparisons = tmp_list,
                     label = "p.format",paired = T)

ggsave(paste(out_dir,"grid.bifido.long.pdf",sep = '/'),width = 5,height = 5)

data_grid_masslin2_l <- read.table("input/z.merge.result.GRiD.txt",header = T,row.names = 1) %>% 
  dplyr::select(tmp_meta_long$sample)

data_grid_masslin2_l <- noise_removal(data_grid_masslin2_l,percent = 0.8,method = "na_cut")
data_grid_masslin2_l <- noise_removal(data_grid_masslin2_l,low = 1,method = "mean_cut")
data_grid_masslin2_l <- apply(data_grid_masslin2_l, c(1,2), function(x){ifelse(is.na(x),0,x)})
Maaslin2(
  input_data = data_grid_masslin2_l, 
  input_metadata = tmp_meta_long, 
  min_prevalence = 0,
  normalization = "NONE",
  transform = "LOG",
  analysis_method = "LM" ,
  output = paste(out_dir,"long.grid",sep = '/'),
  correction = "BH", 
  fixed_effects = c("group1"),
  reference = c("group1,BE"),
  random_effects = c("pair"),
  max_significance = 1
)

data_grid <- as.data.frame(t(data_grid)) %>% filter(row.names(.) %in% tmp_meta_cross$sample)
data_grid$sample <- row.names(data_grid)
data_grid <- left_join(tmp_meta_cross,data_grid,by = "sample")

