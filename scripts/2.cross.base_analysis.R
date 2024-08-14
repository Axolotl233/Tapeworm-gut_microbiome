merge_samples_mean <- function(physeq, group){
  group_sums <- as.matrix(table(sample_data(physeq)[ ,group]))[,1]
  merged <- merge_samples(physeq, group)
  x <- as.matrix(otu_table(merged))
  if(taxa_are_rows(merged)){ x<-t(x) }
  out <- t(x/group_sums)
  out <- otu_table(out, taxa_are_rows = TRUE)
  otu_table(merged) <- out
  return(merged)
}

library(vegan)
#library(aPCoA)
library(factoextra)
library(microbiome)
library(microbiomeSeq)
library(pheatmap)
library(Maaslin2)
dir.create("output/cross_base_analysis",showWarnings = F)
out_dir <- "output/cross_base_analysis"

load("workfile/z.phyloseq.Rdata")

phyloseqin_used <- subset_samples(phyloseqin, groups %in% group[c(1,2)])
tmp_meta <- metadatadf %>% filter( groups %in% group[c(1,2)])
phyloseqin_used <- filter_taxa(phyloseqin_used, function(x) sum(x > 0) > 0, TRUE)
tmp_data_bar <- merge_samples_mean(phyloseqin_used, "groups")
tmp_data_bar_count <- 
  transform_sample_counts(tmp_data_bar, function(x) x/sum(x) * 100)
tmp_data_bar_count <-tax_glom(tmp_data_bar_count, 'Phylum', NArm = T)
tmp_data_bar_count@sam_data$groups <- factor(tmp_data_bar_count@sam_data$groups,levels = group[c(1,2)])

plot_bar(tmp_data_bar_count, fill='Phylum')+
  geom_bar(stat='identity',position='stack')+
  #scale_fill_manual(values = color19)+
  scale_fill_manual(values = color_bar_st)+
  ylim(c(0,101)) + theme_cowplot()+ theme(legend.key=element_blank())
ggsave(paste(out_dir,"cross.bar_phylum.pdf",sep = '/'),width = 6,height = 6)

tmp <- as.data.frame(tmp_data_bar_count@otu_table)
tmp_taxa <- as.data.frame(tmp_data_bar_count@tax_table)
identical(row.names(tmp),row.names(tmp_taxa))
tmp$delat <- tmp[,2]-tmp[,1]
tmp$taxa <- tmp_taxa$Phylum
row.names(tmp) <- 1:nrow(tmp)
write_csv(tmp,paste(out_dir,"cross.bar_phylum.csv",sep = '/'))

phyloseqin_phylum <- tax_glom(phyloseqin_used,taxrank = "Phylum")
tmp_otu_p <- as.data.frame(otu_table(phyloseqin_phylum))
tmp_otu_p <- sweep(tmp_otu_p,2,colSums(tmp_otu_p),'/')*100
tmp_tax <- as.data.frame(tax_table(phyloseqin_phylum))

identical(row.names(tmp_otu_p),row.names(tmp_tax))
row.names(tmp_otu_p) <- tmp_tax$Phylum
tmp_tax_p <- c("Actinobacteria","Bacteroidetes","Firmicutes","Proteobacteria")
tmp_otu_p <- tmp_otu_p %>% mutate(taxon = row.names(tmp_otu_p)) %>%
  filter(taxon %in% tmp_tax_p) %>%
  select(-taxon)
tmp_otu_p <- as.data.frame(t(tmp_otu_p))
tmp_otu_p$sample <- row.names(tmp_otu_p)
tmp_meta_p <- tmp_meta %>% select(sample,groups)
tmp_otu_pl <- inner_join(tmp_otu_p,tmp_meta,by = "sample") %>% 
  gather(key = "taxon",value = "value", tmp_tax_p)

plot_l = list() 
for(i in 1:length(tmp_tax_p)){
  tmp_data <- tmp_otu_pl %>% filter(taxon == tmp_tax_p[i])
  tmp_list <- combine_list(tmp_data$groups)
  p <- ggplot(tmp_data,aes(x=groups,y=value,color = groups,shape =groups))+
    geom_boxplot(outlier.shape = NA,)+
    geom_jitter(size = 2,width = 0.2,height = 0,alpha = 0.7)+
    theme_cowplot()+
    theme(
      legend.position = "none"
    )+
    labs(x="Group",y="",title = tmp_tax_p[i])+
    stat_compare_means(comparisons = tmp_list,
                       method = "wilcox.test", label = "p.format")+
    scale_color_manual(values = color0)
  plot_l[[i]] <- p
}
plot_grid(plot_l[[1]],plot_l[[2]],
          plot_l[[3]],plot_l[[4]],nrow = 2)
ggsave(paste(out_dir,"cross.phylum.pdf",sep = '/'),width = 7,height = 9)
tmp_otu_p$BF_ratio <- tmp_otu_p$Bacteroidetes/tmp_otu_p$Firmicutes
tmp_otu_p2 <- tmp_otu_p %>% select(sample,BF_ratio)
tmp_otu_p2 <- inner_join(tmp_otu_p2,tmp_meta,by = "sample")

ggplot(tmp_otu_p2,aes(x=groups,y=BF_ratio,shape=groups,color = groups))+
  geom_boxplot(outlier.shape = NA,)+
  geom_jitter(size = 2,width = 0.2,height = 0,alpha = 0.7)+
  theme_cowplot()+
  labs(x="Group",y = "Values",title="B/F_ratio")+
  stat_compare_means(comparisons = tmp_list,
                     method = "wilcox.test", label = "p.format")+
  scale_color_manual(values = color0)
ggsave(paste(out_dir,"cross.phylum_bf_ratio.pdf",sep = '/'),width = 5,height = 5)

#=====> diversity analysis
## alpha
phyloseqin_used <- subset_samples(phyloseqin, groups %in% group[c(1,2)])
phyloseqin_used <- filter_taxa(phyloseqin_used, function(x) sum(x > 0) > 0, TRUE)
tmp_otu <- as.data.frame(phyloseqin_used@otu_table)
tmp_sample <- phyloseqin_used@sam_data$sample
tmp_group <- factor(phyloseqin_used@sam_data$groups,levels = group[c(1,2)])
tmp_sex <- phyloseqin_used@sam_data$sex
tmp_age <- phyloseqin_used@sam_data$age
tmp_bmi <- phyloseqin_used@sam_data$bmi

tmp_shannon <- vegan::diversity(t(tmp_otu),index = "shannon")
tmp_simpson <-  vegan::diversity(t(tmp_otu),index = "simpson")
tmp_invsimpson <- 1/tmp_simpson

tmp_alpha_raw <- data.frame(sample=tmp_sample
                             ,type=tmp_group
                             ,sex=tmp_sex
                            ,age=tmp_age
                            ,bmi=tmp_bmi
                            ,shannon=tmp_shannon
                            ,simpson=tmp_simpson
                            ,invsimpson=tmp_invsimpson
)

tmp_observe <- data.frame(sample = "",observed ="")
for(i in 1:ncol(tmp_otu)){
  tmp_sam = colnames(tmp_otu)[i]
  tmp_ob = length(which(tmp_otu[,i] != 0))
  tmp_observe[i,] = c(tmp_sam,tmp_ob)
}
tmp_alpha_raw <- inner_join(tmp_alpha_raw,tmp_observe,by = "sample")
tmp_alpha_raw$observed <- as.numeric(tmp_alpha_raw$observed)
write_csv(x = tmp_alpha_raw,paste(out_dir,"cross.alpha_diversity.csv",sep = '/'))

phyloseqin_used2 <- subset_samples(phyloseqin_count, groups %in% group[c(1,2)])
richness_index <- c("Observed","Shannon","Simpson","InvSimpson")
tmp_other <- estimate_richness(phyloseqin_used2,measures = richness_index)
tmp_other$sample <- row.names(tmp_other)
tmp_alpha_raw2 <- data.frame(sample=tmp_sample
                            ,type=tmp_group
                            ,sex=tmp_sex
                            ,age=tmp_age
                            ,bmi=tmp_bmi
)
tmp_alpha_raw2 <- left_join(tmp_alpha_raw2,tmp_other,by = "sample")
write_csv(x = tmp_alpha_raw2,paste(out_dir,"cross.alpha_diversity_phyloseq.csv",sep = '/'))

tmp_alpha_raw_plot <- tmp_alpha_raw %>% dplyr::select(sample,type,observed,shannon,simpson) %>%
  gather(key = "key",value = "value",observed,shannon,simpson)

tmp_list <- combine_list(tmp_alpha_raw_plot$type)
ggplot(tmp_alpha_raw_plot,aes(x=type,y=value)) + 
  geom_boxplot(aes(color=type),outlier.shape = NA)+
  geom_jitter(aes(color=type,shape=type),size = 2.5,height=0,width = 0.2,alpha = 0.7)+
  theme_bw()+
  theme(axis.title.x=element_blank(), 
        legend.position="none",
        panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
  scale_color_manual(values = color0)+
  labs(y="Observed value")+
  stat_compare_means(
    comparisons = tmp_list,
    method = "wilcox.test", label = "p.format")+
  facet_wrap(.~key,scales = "free")
ggsave(paste(out_dir,"cross.alpha_diversity_raw.pdf",sep = '/'),width = 6,height = 5)

row.names(tmp_alpha_raw) <- tmp_alpha_raw$sample
tmp_data <- as.data.frame(t(tmp_alpha_raw[,c(6:9)]))
Maaslin2(
  input_data = tmp_data, 
  input_metadata = tmp_alpha_raw, 
  min_prevalence = 0,
  normalization = "NONE",
  transform = "LOG",
  analysis_method = "LM" ,
  output = paste(out_dir,"cross.alpha_cov_adjust",sep = '/'), 
  fixed_effects = c("type","sex","age","bmi"),
  reference = c("type,HC")
)

## beta
#tmp_otu2 <- noise_removal(tmp_otu,percent = 0.1,low = 0.001)
phyloseqin_pca <- subset_samples(phyloseqin, groups %in% c(group[1],group[2]))
phyloseqin_pca_f <- filter_taxa(phyloseqin_pca, function(x) sum(x > 0.005) > (0.2*length(x)), TRUE)
tmp_data_pca <- as.data.frame(t(phyloseqin_pca_f@otu_table))
all(row.names(tmp_data_pca) == phyloseqin_pca@sam_data$sample)
tmp_group <- factor(phyloseqin_pca@sam_data$groups,levels = group[c(1,2)])
tmp_res_pca <- prcomp(tmp_data_pca, scale = T)
tmp_dist = vegdist(tmp_data_pca, method = 'bray')
tmp_sex <- factor(phyloseqin_pca@sam_data$sex,levels = c("female","male"))
tmp_age <- phyloseqin_pca@sam_data$age
tmp_bmi <- phyloseqin_pca@sam_data$bmi
tmp_site = data.frame(sample = rownames(tmp_data_pca),group = tmp_group,sex =tmp_sex,age = tmp_age,bmi = tmp_bmi)
tmp_adonis_result_dis = adonis2(tmp_dist~group+age+sex+bmi, tmp_site)
fviz_pca_ind(tmp_res_pca, label="none",axes = c(1, 2), alpha.ind =1,
             habillage=tmp_group,invisible="quali",pointsize = 2.5, addEllipses = TRUE,
             ellipse.level=0.95,palette = color0)+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(vjust = 0.5, hjust = 0.1)
  )+
  annotate("text", x=min(tmp_res_pca$x[,1])+3, y=min(tmp_res_pca$x[,2])+17,
           label=paste("p-value(adonis)=",tmp_adonis_result_dis$`Pr(>F)`[1]))+
  annotate("text", x=min(tmp_res_pca$x[,1])+3, y=min(tmp_res_pca$x[,2])+18, 
           label=paste("R2(adonis)=",round(tmp_adonis_result_dis$R2[1],3)))

ggsave(paste(out_dir,"cross.beta_pca.pdf",sep = '/'),width = 6,height = 5)
tmp_anosim_result<-anosim(tmp_data_pca,tmp_group)
pdf(paste(out_dir,"cross.beta_anosim.pdf",sep = '/'),width = 6,height = 4)
plot(tmp_anosim_result, col = c('#FFD700',"#00B5E2","#d95f02","#7570b3"))
dev.off()

#tmp_site = data.frame(sample = rownames(tmp_data_pca),group = tmp_group,sex=tmp_sex)
tmp_pcoa_res <- ordinate(phyloseqin_pca_f, method = "PCoA", distance = "bray")
tmp_pcoa_res_table <- as.data.frame(tmp_pcoa_res$vectors[,1:3])
colnames(tmp_pcoa_res_table) <- c("dim1","dim2","dim3")
all(row.names(tmp_pcoa_res_table) == tmp_meta$sample)
tmp_pcoa_res_table$group <-  tmp_meta$groups

#row.names(tmp_site) <- tmp_site$sample
#tmp_apcoa_res <- aPCoA(tmp_dist ~ age+sex+bmi,tmp_site,maincov = group)
#tmp_apcoa_matrix <- data.frame(tmp_apcoa_res$plotMatrix)
#identical(row.names(tmp_apcoa_matrix),tmp_site$sample)
#tmp_apcoa_matrix$group <- tmp_site$group
#tmp_apcoa_matrix$tmp = "tmp"

tmp_list <- combine_list(tmp_pcoa_res_table$group)
a <- ggplot()+
  geom_point(data = tmp_pcoa_res_table,
             mapping = aes(x=dim1,y=dim2,color=group,shape=group),
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
  stat_ellipse(data=tmp_pcoa_res_table,geom = "polygon",
               aes(x=dim1,y=dim2,fill=group),
               color =NA,alpha=0.1)+
  scale_fill_manual(values = color0)+
  scale_color_manual(values = color0)+
  scale_x_continuous(limits = c(-0.5,0.5))+
  scale_y_continuous(limits = c(-0.5,0.5))

b <- ggplot(data=tmp_pcoa_res_table,aes(x=group,y=dim1,fill=group))+
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
  scale_y_continuous(limits = c(-0.5,0.5))+
  stat_compare_means(
    comparisons = tmp_list,
    method = "wilcox.test", label = "p.format")

c <- ggplot(data=tmp_pcoa_res_table,aes(x=group,y=dim2,fill=group))+
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
  scale_y_continuous(limits = c(-0.5,0.5))+
  stat_compare_means(
    comparisons = tmp_list,method = "wilcox.test", label = "p.format"
    )

tmp_adonis_result_dis$item <- row.names(tmp_adonis_result_dis)
tmp_adonis_result_dis_p <- tmp_adonis_result_dis[!is.na(tmp_adonis_result_dis$F),]
p1 <- ggplot(tmp_adonis_result_dis_p)+
  geom_bar(aes(x = item, y = `F`,fill = item),stat = "identity",width = 0.6)+
  geom_text(aes(x = item, y = `F`+ 0.2,label= paste("P = ",`Pr(>F)`,sep = "")))+
  theme_cowplot()+
  scale_fill_manual(values = c("grey80","grey80","#5698c4","grey80"))+
  theme(legend.position = "none")
p2 <- ggplot(tmp_adonis_result_dis_p)+
  geom_bar(aes(x = item, y = R2,fill = item),stat = "identity",width = 0.6)+
  geom_text(aes(x = item, y = R2 + 0.002,label= paste("P = ",`Pr(>F)`,sep = "")))+
  theme_cowplot()+
  scale_fill_manual(values = c("grey80","grey80","#5698c4","grey80"))+
  theme(legend.position = "none")
p1
ggsave(plot = p1,filename = paste(out_dir,"cross.beta_adonis_F.pdf",sep = '/'),width = 3,height = 4)
tmp_adonis_result_dis_p <- as.data.frame(t(tmp_adonis_result_dis_p))
tmp_adonis_result_dis_p <- tmp_adonis_result_dis_p[3:5,]
tmp_adonis_result_dis_p$posy <- c(0.25,0.5,0.75)
tmp_adonis_result_dis_p$posx <- 1
tmp_adonis_result_dis_p$labels <- paste(row.names(tmp_adonis_result_dis_p),round(as.numeric(tmp_adonis_result_dis_p[,1]),3),sep = ':')

d <- ggplot()+
  geom_bar(data = as.data.frame(matrix(c(1),nrow = 1)),aes(x = V1,y = V1),stat ="identity",fill="white")+
  geom_text(data = tmp_adonis_result_dis_p,mapping = aes(x = posx,y = posy,label = labels))+
  theme_bw()+theme(
    legend.position = "NONE",
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
#plot_grid(b,d,a,c,align = "hv",rel_widths = c(1,0.225),rel_heights = c(0.225,1))
aplot::plot_list(b,d,a,c,widths = c(1,0.225),heights = c(0.225,1))
ggsave(paste(out_dir,"cross.beta_pcoa.pdf",sep = '/'),width = 5,height = 5)

## distance
ann_sam <-data.frame(phyloseqin_used@sam_data$groups)
row.names(ann_sam) <- phyloseqin_used@sam_data$sample
colnames(ann_sam) <- c("type")

ann_colors <- list(type=c(HC=color0[1],TA=color0[2]))
dist <- distance(phyloseqin_used,method="bray")
dist_matrix <- as.matrix(dist)
tmp_hclust <- hclust(dist,method="complete") 
plot(tmp_hclust)   
rect.hclust(tmp_hclust,k=2)    
pheatmap(
  dist_matrix,
  cluster_rows = T,
  cluster_cols = T,
  clustering_method ="complete",
  #scale = "row",
  annotation_col = ann_sam,
  annotation_colors=ann_colors,
  color = colorRampPalette(colors = c("blue","white","red"))(200)
)
ann_sam2 <- as.data.frame(ann_sam)
ann_sam2$sample <- rownames(ann_sam2)

group1 <- ann_sam2[ann_sam2$type == group[1],2]
group2 <- ann_sam2[ann_sam2$type == group[2],2]

dist_group1 <- as.data.frame(dist_matrix[group1,group1])
dist_group2 <- as.data.frame(dist_matrix[group2,group2])
dist_inter <- as.data.frame(dist_matrix[group1,group2])

dist_group1$sample <- rownames(dist_group1)
dist_group2$sample <- rownames(dist_group2)
dist_inter$sample <- rownames(dist_inter)

dist_group1_long <- gather(dist_group1,"sample2","value",-sample) %>% 
  mutate(class = group[1])
dist_group2_long <- gather(dist_group2,"sample2","value",-sample) %>% 
  mutate(class = group[2])
dist_inter_long <- gather(dist_inter,"sample2","value",-sample) %>% 
  mutate(class = "Inter")

dist_plot_type2 <- bind_rows(dist_group1_long,dist_group2_long,dist_inter_long) %>% 
  filter(sample!=sample2)
dist_plot_type2$class <- factor(dist_plot_type2$class,levels = c("Inter",group[1],group[2]))

ggplot(dist_plot_type2,aes(x=class,y=value))+
  geom_violin(aes(fill = class),color = "white")+
  geom_boxplot(width = 0.3)+
  stat_compare_means(method="wilcox.test",
                     comparisons = list(c(group[1],group[2]),
                                        c("Inter",group[1]),
                                        c(group[2],"Inter")),
                     label = "p.format"
  )+
  labs(y = "Bray-Curtis distance")+
  theme_cowplot()+
  theme(legend.position="none", 
        strip.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_fill_manual(
    values = c("gray",color0[c(1,2)])
  )
ggsave(paste(out_dir,"cross.distance_bray.pdf",sep = '/'),width = 4,height = 4)


rm(list=ls(pattern = "tmp"))
rm(list=ls(pattern = "dist"))
rm(list=ls(pattern = "ann"))
rm(list=ls(pattern = "phyloseqin"))
