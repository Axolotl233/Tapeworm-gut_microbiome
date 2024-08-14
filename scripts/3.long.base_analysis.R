library(vegan)
library(microbiome)
library(microbiomeSeq)
library(Maaslin2)
library(car)

dir.create("output/long_base_analysis",showWarnings = F)
out_dir <- "output/long_base_analysis"
load("workfile/z.phyloseq.Rdata")

#=====> diversity analysis
## alpha
phyloseqin_used <- subset_samples(phyloseqin, group1 %in% c(groupl[2],groupl[3]))
tmp_otu <- as.data.frame(phyloseqin_used@otu_table)
tmp_sample <- phyloseqin_used@sam_data$sample
tmp_group <- factor(phyloseqin_used@sam_data$group1,levels = groupl[c(2,3)])
tmp_age <- phyloseqin_used@sam_data$age
tmp_sex <- phyloseqin_used@sam_data$sex
tmp_pair <- phyloseqin_used@sam_data$pair
tmp_bmi <- phyloseqin_used@sam_data$bmi

tmp_shannon <- vegan::diversity(t(tmp_otu),index = "shannon")
tmp_simpson <-  vegan::diversity(t(tmp_otu),index = "simpson")
tmp_invsimpson <- 1/tmp_simpson

tmp_alpha_raw <- data.frame(sample=tmp_sample
                            ,type=tmp_group
                            ,pair=tmp_pair
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
write_csv(tmp_alpha_raw,paste(out_dir,"long.alpha_diversity.csv",sep = '/'))

identical(tmp_alpha_raw$pair[1:(nrow(tmp_alpha_raw)/2)],
          tmp_alpha_raw$pair[((nrow(tmp_alpha_raw)/2)+1):nrow(tmp_alpha_raw)])

wilcox.test(tmp_alpha_raw$shannon[1:(nrow(tmp_alpha_raw)/2)],
            tmp_alpha_raw$shannon[((nrow(tmp_alpha_raw)/2)+1):nrow(tmp_alpha_raw)],
            paired = T)

wilcox.test(tmp_alpha_raw$simpson[1:(nrow(tmp_alpha_raw)/2)],
            tmp_alpha_raw$simpson[((nrow(tmp_alpha_raw)/2)+1):nrow(tmp_alpha_raw)],
            paired = T)

wilcox.test(tmp_alpha_raw$observed[1:(nrow(tmp_alpha_raw)/2)],
            tmp_alpha_raw$observed[((nrow(tmp_alpha_raw)/2)+1):nrow(tmp_alpha_raw)],
            paired = T)
tmp_alpha_raw_plot <- tmp_alpha_raw %>% dplyr::select(sample,type,pair,observed,shannon,simpson) %>%
  gather(key = "key",value = "value",observed,shannon,simpson)
tmp_list <- combine_list(tmp_alpha_raw_plot$type)

ggplot(tmp_alpha_raw_plot,aes(x=type,y=value)) + 
  geom_boxplot(aes(color=type),outlier.shape = NA)+
  geom_point(aes(color=type,shape=type),size = 2.5,alpha = 0.7)+
  geom_line(aes(group=pair),color = "grey80",linetype = "dashed",size = 0.3)+
  theme_bw()+
  theme(axis.title.x=element_blank(), 
        legend.position="none",
        panel.grid.minor = element_blank(),panel.grid.major = element_blank())+
  scale_color_manual(values = color0[c(2,3)])+
  scale_shape_manual(values = c(17,15))+
  labs(y="Observed value")+
  stat_compare_means(
    comparisons = tmp_list,
    method = "wilcox.test", label = "p.format")+
  facet_wrap(.~key,scales = "free")
ggsave(paste(out_dir,"long.alpha_diversity_raw.pdf",sep = '/'),width = 6,height = 5)

row.names(tmp_alpha_raw) <- tmp_alpha_raw$sample
tmp_data <- as.data.frame(t(tmp_alpha_raw[,c(7:10)]))
Maaslin2(
  input_data = tmp_data, 
  input_metadata = tmp_alpha_raw, 
  min_prevalence = 0,
  normalization = "NONE",
  transform = "LOG",
  analysis_method = "LM" ,
  output = paste(out_dir,"long.alpha_cov_adjust",sep = '/'), 
  fixed_effects = c("type","sex","age","bmi"),
  reference = c("type,BE"),
  random_effects = c("pair")
)

group1 <- metadatadf[metadatadf$group1 == "BE",c("sample","pair")]
group2 <- metadatadf[metadatadf$group1 == "AF",c("sample","pair")]
group_all <- inner_join(group1,group2,by = "pair")
tmp_delat_alpha <- data.frame(pair = "",shannon = "",simpson = "",observed = "")

for(i in 1:nrow(group_all)){
  tmp_data <- tmp_alpha_raw %>% filter(sample %in% c(group_all[i,1],group_all[i,3]))
  tmp_delat_alpha[i,] = c(group_all[i,2],tmp_data[2,7]-tmp_data[1,7],
                          tmp_data[2,8]-tmp_data[1,8],tmp_data[2,10]-tmp_data[1,10])
}
tmp_delat_alpha[,2:4] <- apply(tmp_delat_alpha[,2:4],2,as.numeric)

length(which(tmp_delat_alpha$shannon> 0))
length(which(tmp_delat_alpha$simpson> 0))
length(which(tmp_delat_alpha$observed> 0))
tmp_formula = y ~ x
ggplot(tmp_delat_alpha,aes(x=shannon,y=observed)) + 
  stat_smooth(method='lm',formula = tmp_formula,color = "grey90")+
  geom_point(size=4,color = "firebrick3",alpha = 0.7)+
  stat_cor(cor.coef.name = "rho",method = "spearman",
           label.y.npc = 1,label.x.npc = 0.15)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )


## beta
tmp_meta <- metadatadf %>% filter( group1 %in% groupl[c(2,3)])
phyloseqin_used <- subset_samples(phyloseqin, group1 %in% c(groupl[2],groupl[3]))
phyloseqin_used <- filter_taxa(phyloseqin_used, function(x) sum(x > 0.001) > (0.2*length(x)), TRUE)

tmp_pcoa_res <- ordinate(phyloseqin_used, "PCoA", "bray", weighted =FALSE)
tmp_pcoa_res_table <- as.data.frame(tmp_pcoa_res$vectors)
tmp_pcoa_res_gr1 <- tmp_pcoa_res_table[row.names(tmp_pcoa_res_table) %in% tmp_meta$sample[tmp_meta$group1 == groupl[2]],]
tmp_pcoa_res_gr2 <- tmp_pcoa_res_table[row.names(tmp_pcoa_res_table) %in% tmp_meta$sample[tmp_meta$group1 == groupl[3]],]
#summary( procrustes(tmp_pcoa_res_gr1,tmp_pcoa_res_gr2,symmetric = T))

pro_pcoa <- protest(tmp_pcoa_res_gr1,tmp_pcoa_res_gr2, symmetric = T, permutations = 999)
summary(pro_pcoa)

tmp_pcoa_res_table <- tmp_pcoa_res_table[c(1:3)]
colnames(tmp_pcoa_res_table) <- c("dim1","dim2","dim3")
all(row.names(tmp_pcoa_res_table) == tmp_meta$sample)
tmp_pcoa_res_table$group <-  tmp_meta$group1
tmp_pcoa_res_table$pair <-  tmp_meta$pair
tmp_list <- combine_list(tmp_pcoa_res_table$group)

a <- ggplot(data=tmp_pcoa_res_table,aes(x=dim1,y=dim2))+
  #geom_vline(xintercept = 0,lty="dashed",linewidth = 0.3)+
  #geom_hline(yintercept = 0,lty="dashed",linewidth = 0.3)+
  geom_polygon(aes(fill = group),stat = "ellipse",alpha=0.05)+
  geom_line(aes(group=pair),color = "grey80",linetype = "dashed",linewidth = 0.3)+
  geom_point(aes(shape=group,color=group),size=2.5,alpha = 0.7)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.1,0.1),
        axis.text.y = element_text(angle = 90,hjust = 0.5,vjust = 1)
  )+
  labs(x=c("PCoA1"),y=c("PCoA2"))+
  scale_fill_manual(values = color0[c(2,3)])+
  scale_color_manual(values = color0[c(2,3)])+
  scale_shape_manual(values = c(17,15))+
  scale_x_continuous(limits = c(-0.6,0.6))+
  scale_y_continuous(limits = c(-0.6,0.6))

b <- ggplot(data=tmp_pcoa_res_table,aes(x=group,y=dim1,fill=group))+
  #geom_violin(aes(x=tmp,y=X1,fill=group),width = 1,trim = F,color = "white")+
  geom_line(aes(group=pair),color = "grey80",linetype = "dashed",linewidth = 0.3)+
  geom_point(aes(shape=group,color=group),size=1.5,alpha = 0.7)+
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
  scale_color_manual(values = color0[c(2,3)])+
  scale_shape_manual(values = c(17,15))+
  scale_y_continuous(limits = c(-0.6,0.6))+
  stat_compare_means(
    comparisons = tmp_list,
    method = "wilcox.test", label = "p.format")
c <- ggplot(data=tmp_pcoa_res_table,aes(x=group,y=dim2,fill=group))+
  #geom_violin(aes(x=tmp,y=X1,fill=group),width = 1,trim = F,color = "white")+
  geom_line(aes(group=pair),color = "grey80",linetype = "dashed",linewidth = 0.3)+
  geom_point(aes(shape=group,color=group),size=1.5,alpha = 0.7)+
  theme_bw()+
  theme(
    legend.position = "NONE",
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )+
  labs(x=element_blank(),y=element_blank())+
  scale_color_manual(values = color0[c(2,3)])+
  scale_shape_manual(values = c(17,15))+
  scale_y_continuous(limits = c(-0.6,0.6))+
  stat_compare_means(
    comparisons = tmp_list,method = "wilcox.test", label = "p.format"
  )

d <- ggplot()+
  geom_bar(data = as.data.frame(matrix(c(1),nrow = 1)),aes(x = V1,y = V1),stat ="identity",fill="white")+
  annotate("text", x=1, y=0.7,label=paste("M2=",round(pro_pcoa$ss,3)))+
  annotate("text", x=1, y=0.3,label=paste("P=",pro_pcoa$signif))+
  theme_bw()+theme(
    legend.position = "NONE",
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
#plot_grid(b,d,a,c,align = "hv",rel_widths = c(1,0.225),rel_heights = c(0.225,1))
aplot::plot_list(b,d,a,c,widths = c(1,0.225),heights = c(0.225,1))
ggsave(paste(out_dir,"long.beta_pcoa.pdf",sep = '/'),width = 5,height = 5)

wilcox.test(tmp_pcoa_res_table$dim1[1:(nrow(tmp_pcoa_res_table)/2)],
            tmp_pcoa_res_table$dim1[((nrow(tmp_pcoa_res_table)/2)+1):nrow(tmp_pcoa_res_table)],
            paired = T)
wilcox.test(tmp_pcoa_res_table$dim2[1:(nrow(tmp_pcoa_res_table)/2)],
            tmp_pcoa_res_table$dim2[((nrow(tmp_pcoa_res_table)/2)+1):nrow(tmp_pcoa_res_table)],
            paired = T)


#group_dedim <- function(phylo,gr){
#  tmp_data_pca <- as.data.frame(t(phylo@otu_table))
# tmp_res_pca <- prcomp(tmp_data_pca, scale = T)
#  tmp_res_pcoa <- ordinate(phylo, "PCoA", "bray", weighted =FALSE)
#  tmp_res <- list(tmp_res_pca,tmp_res_pcoa,phylo@sam_data$pair)
#  return(tmp_res)
#}
#tmp_dedim_gr1 <- group_dedim(subset_samples(phyloseqin_used, group1 %in% c(groupl[2])))
#tmp_dedim_gr2 <- group_dedim(subset_samples(phyloseqin_used, group1 %in% c(groupl[3])))
#all(tmp_dedim_gr1[[3]] == tmp_dedim_gr2[[3]])

#tmp_pcoa_res_gr1 <- tmp_dedim_gr1[[2]]$vectors
#tmp_pcoa_res_gr2 <- tmp_dedim_gr2[[2]]$vectors

#pro_pcoa <- protest(tmp_pcoa_res_gr1,tmp_pcoa_res_gr2, symmetric = TRUE, permutations = 999)
#plot(pro_pcoa, kind = 2)
#residuals(pro_pcoa)

## distance
#load("workfile/z.phyloseq.Rdata")
dist <- distance(phyloseqin,method="bray")
dist_matrix <- as.matrix(dist)
group1 = group_all$sample.x
group2 = group_all$sample.y
group3 <- metadatadf[metadatadf$groups == "HC",2]

dist_inter1 <- as.data.frame(dist_matrix[group1,group3])
dist_inter2 <- as.data.frame(dist_matrix[group2,group3])

dist_inter1$sample <- rownames(dist_inter1)
dist_inter2$sample <- rownames(dist_inter2)

dist_inter1_long <- gather(dist_inter1,"sample2","value",-sample) %>% 
  mutate(class = "HC_BE")
dist_inter2_long <- gather(dist_inter2,"sample2","value",-sample) %>% 
  mutate(class = "HC_AF")

dist_pair <- as.data.frame(dist_matrix[group1,group2])
identical(row.names(dist_pair),group1)
identical(colnames(dist_pair),group2)
dist_pair_long <- data.frame(sample="",sample2="",value = 0,class ="")
for(i in 1:nrow(dist_pair)){
  dist_pair_long[i,] = c(group1[i],group2[i],dist_pair[i,i],group_all[i,2])
}
dist_pair_long$value <- as.numeric(dist_pair_long$value)
write_csv(x = dist_pair_long,paste(out_dir,"long.distance_BC_paired.csv",sep = '/'))
#dist_pair_long <- read_csv("output/long_base_analysis/long.distance_BC_paired.csv")
dist_pair_long$class <- "pair"
dist_plot_type2 <- bind_rows(dist_inter1_long,dist_inter2_long,dist_pair_long) %>% 
  filter(sample!=sample2)
dist_plot_type2$class <- factor(dist_plot_type2$class,levels = c("HC_BE","HC_AF","pair"))

ggplot(dist_plot_type2,aes(x=class,y=value))+
  geom_violin(aes(fill = class),color = "white")+
  geom_boxplot(width = 0.3)+
  stat_compare_means(method="wilcox.test",
                     comparisons = list(c("HC_BE","HC_AF"),
                                        c("pair","HC_AF"),
                                        c("HC_BE","pair")),
                     label = "p.format"
  )+
  ylab("Bray-Curtis distance")+
  theme_cowplot()+
  theme(legend.position="none", 
        strip.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        )+
  scale_fill_manual(
    values = c(color0[c(2,3)],"gray80")
  )
wilcox.test(dist_inter1_long$value,dist_inter2_long$value,paired = T)
ggsave(paste(out_dir,"long.distance_bray.pdf",sep = '/'),width = 4,height = 4)

dist_inter1_mean <- dist_inter1_long %>% group_by(sample) %>% summarise(value=mean(value))
dist_inter1_mean$class <- "HC_BE"
dist_inter2_mean <- dist_inter2_long %>% group_by(sample) %>% summarise(value=mean(value))
dist_inter2_mean$class <- "HC_AF"
dist_plot_type3 <- bind_rows(dist_inter1_mean,dist_inter2_mean,dist_pair_long[,c(1,3,4)])

ggplot(dist_plot_type3,aes(x=class,y=value))+
  geom_violin(aes(fill = class),color = "white")+
  geom_boxplot(width = 0.3)+
  stat_compare_means(method="wilcox.test",
                     comparisons = list(c("HC_BE","HC_AF"),
                                        c("pair","HC_AF"),
                                        c("HC_BE","pair")),
                     label = "p.format"
  )+
  ylab("Bray-Curtis distance")+
  theme_cowplot()+
  theme(legend.position="none", 
        strip.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
  )+
  scale_fill_manual(
    values = c(color0[c(2,3)],"gray80")
  )

ggplot(dist_pair_long)+
  geom_histogram(aes(x = value),color = "navy",fill = "navy",alpha = 0.7,binwidth = 0.05,boundary = 0)+
  #geom_density(aes(x = value),color = "navy",alpha = 0.7,binwidth = 0.05)+
  theme_cowplot()+
  scale_x_continuous(breaks = seq(0,1,0.1),limits = c(0,1))+
  labs(x = "Bray-Curtis distance" , y = "Density")
ggsave(paste(out_dir,"long.distance_bray_histogram.pdf",sep = '/'),width = 4,height = 4)
#metadatadf2 <- read_table(in_meta_2)
#tmp_res <- data.frame(p1 = "",p2 = "",p3 ="")
#for(i in 1:1000){
#  tmp_data1 <- dist_inter1_long[sample(1:nrow(dist_inter1_long),24),]
#  tmp_sample1 <- tmp_data1$sample
#  tmp_sample2 <- tmp_data1$sample2
#  tmp_sample3 <- metadatadf2[metadatadf2$sample1 %in% tmp_sample1,"sample2"]
#  tmp_res <- wilcox.test(tmp_data$value,data2$value)
#  result[i,] <- c(tmp_res$p.value)
#}

tmp_meta <- metadatadf %>% select(sample,interval,interval_class)
dist_pair_other <- inner_join(dist_pair_long,tmp_meta,by = 'sample')
dist_pair_other <- inner_join(dist_pair_other,tmp_alpha_raw,by = "sample")# %>% filter(sample != "T65be")
dist_pair_other <- dist_pair_other %>% filter(shannon > 2)

tmp_formula <- y ~ x
p1 <- ggplot(dist_pair_other,aes(x=interval,y=value)) + 
  stat_smooth(method='lm',formula = tmp_formula,color = "grey90")+
  geom_point(size=2.5,color = "navy")+
  stat_cor(cor.coef.name = "rho",method = "spearman",
           label.y.npc = 1,label.x.npc = 0.15)+
  #stat_poly_eq(aes(label = ..eq.label..),formula = tmp_formula,label.y.npc = 0.95,label.x.npc=0.87, parse = TRUE,)+
  theme_light()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )+
  scale_y_continuous(#breaks = seq(0,400,100),
    breaks = seq(0,1,0.1),limits = c(0.2,1)
    #labels = scales::scientific
    )+
  #scale_x_continuous(breaks = seq(-2.5,7.5,2.5),limits = c(-3,7.5))+
  labs(y="Bray-Cruit distance",x="time interval")
p2 <- ggplot(dist_pair_other,aes(x=shannon,y=value)) + 
  stat_smooth(method='lm',formula = tmp_formula,color = "grey90")+
  geom_point(size=2.5,color = "navy")+
  stat_cor(cor.coef.name = "rho",method = "spearman",
           label.y.npc = 1,label.x.npc = 0.15)+
  #stat_poly_eq(aes(label = ..eq.label..),formula = tmp_formula,label.y.npc = 0.95,label.x.npc=0.87, parse = TRUE,)+
  theme_light()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )+
  scale_y_continuous(
    breaks = seq(0,1,0.1),limits = c(0.2,1)
    #labels = scales::scientific
    )+
  #scale_x_continuous(breaks = seq(-2.5,7.5,2.5),limits = c(-3,7.5))+
  labs(y="Bray-Cruit distance",x="Shannon index")
p3 <- ggplot(dist_pair_other,aes(x=simpson,y=value)) + 
  stat_smooth(method='lm',formula = tmp_formula,color = "grey90")+
  geom_point(size=2.5,color = "navy")+
  stat_cor(cor.coef.name = "rho",method = "spearman",
           label.y.npc = 1,label.x.npc = 0.15)+
  #stat_poly_eq(aes(label = ..eq.label..),formula = tmp_formula,label.y.npc = 0.95,label.x.npc=0.87, parse = TRUE,)+
  theme_light()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )+
  scale_y_continuous(#breaks = seq(0,400,100),
    breaks = seq(0,1,0.1),limits = c(0.2,1)
  )+
  #scale_x_continuous(breaks = seq(-2.5,7.5,2.5),limits = c(-3,7.5))+
  labs(y="Bray-Cruit distance",x="Simpson index")
p4 <- ggplot(dist_pair_other,aes(x=observed,y=value)) + 
  stat_smooth(method='lm',formula = tmp_formula,color = "grey90")+
  geom_point(size=2.5,color = "navy")+
  stat_cor(cor.coef.name = "rho",method = "spearman",
           label.y.npc = 1,label.x.npc = 0.15)+
  #stat_poly_eq(aes(label = ..eq.label..),formula = tmp_formula,label.y.npc = 0.95,label.x.npc=0.87, parse = TRUE,)+
  theme_light()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )+
  scale_y_continuous(#breaks = seq(0,400,100),
    breaks = seq(0,1,0.1),limits = c(0.2,1)
  )+
  #scale_x_continuous(breaks = seq(-2.5,7.5,2.5),limits = c(-3,7.5))+
  labs(y="Bray-Cruit distance",x="Observed taxa")
plot_grid(p1,p4,p2,p3,nrow = 2)
ggsave(paste(out_dir,"long.distance2alpha.pdf",sep = '/'),width = 6,height = 6)

p5 <- ggplot(dist_pair_other)+
  geom_histogram(aes(x = value),color = "navy",fill = "navy",alpha = 0.7,binwidth = 0.05,boundary = 0)+
  #geom_density(aes(x = value),color = "navy",alpha = 0.7,binwidth = 0.05)+
  theme_bw()+
  scale_x_continuous(breaks = seq(0,1,0.1),limits = c(0.2,1))+
  theme(
    legend.position = "NONE",
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )+
  labs(x=element_blank(),y=element_blank())+
  coord_flip()
p6 <- ggplot(dist_pair_other)+
  geom_histogram(aes(x = shannon),color = "navy",fill = "navy",alpha = 0.7,binwidth = 0.1,boundary = 0)+
  #geom_density(aes(x = value),color = "navy",alpha = 0.7,binwidth = 0.05)+
  theme_bw()+
  theme(
    legend.position = "NONE",
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )+
  labs(x=element_blank(),y=element_blank())
d <- ggplot()+
  theme_nothing()+theme(
    legend.position = "NONE",
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
aplot::plot_list(p6,d,p2,p5,widths = c(1,0.225),heights = c(0.225,1))
ggsave(paste(out_dir,"long.distance2shannon.pdf",sep = '/'),width = 4,height = 4)

#cor.test(dist_pair_other$value,dist_pair_other$shannon,method = "spearman")
#cor.test(dist_pair_other$value,dist_pair_other$simpson,method = "spearman")
all_otu <- as.data.frame(otu_table(phyloseqin_all))
all_otu <- all_otu[grep("t__",row.names(all_otu),invert = T),dist_pair_long$sample]
all_otu <- noise_removal(all_otu,percent = 0.2)
all_otu <- as.data.frame(t(all_otu))
identical(row.names(all_otu),dist_pair_long$sample)
tmp_res <- data.frame(sp = "" , r_pearson = "", p_pearson="",r_spearman = "",p_spearman = "")
for(i in 1:ncol(all_otu)){
  tmp_n <- colnames(all_otu)[i]
  tmp_data <- all_otu[,i]
  tmp_res1 <- cor.test(tmp_data,dist_pair_long$value)
  tmp_res2 <- cor.test(tmp_data,dist_pair_long$value,method = "spearman")
  tmp_res[i,] <- c(tmp_n,
                   tmp_res1$estimate,tmp_res1$p.value,
                   tmp_res2$estimate,tmp_res2$p.value
                   )
}
tmp_res[,2:5] <- apply(tmp_res[,2:5],2,as.numeric)
tmp_res <- tmp_res[grepl("s__",tmp_res$sp),]
tmp_res$p_pearson_adj <- p.adjust(tmp_res$p_pearson,method = "BH")
tmp_res$p_spearman_adj <- p.adjust(tmp_res$p_spearman,method = "BH")
write_csv(tmp_res,paste(out_dir,"long.BC2sp.csv",sep = '/'))

rm(list=ls(pattern = "tmp"))
rm(list=ls(pattern = "dist"))
rm(list=ls(pattern = "ann"))
rm(list=ls(pattern = "phyloseqin"))
