pathway_beta <- function(x,y,z){
  out1 <- paste("output/long",z,"pca.pdf",sep = ".")
  out2 <- paste("output/long",z,"pcoa.pdf",sep = ".")
  data_path_pca=x
  data_path_pca <- noise_removal(data_path_pca,percent = 0.1,low = 0.001)
  data_path_pca <- as.data.frame(t(data_path_pca))
  data_path_pca$sample <- row.names(data_path_pca)
  data_path_pca <- inner_join(data_path_pca,y,by = "sample")
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
         title = paste(z,"Indiviuals pathway PCA",sep = " ")
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
library(aPCoA)
library(microeco)
library(file2meco)
library(magrittr)
library(aplot)

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
  output = "output/long.pathway_all_metacyc", 
  fixed_effects = c("group1","interval"),
  reference = c("group1,before"),
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
  output = "output/5.pathway_long_s", 
  fixed_effects = c("group1","interval"),
  reference = c("group1,before"),
  random_effects = c("pair"),
)

meco_pathway_metacyc <- humann2meco(in_pathway_metacyc, 
                                    db = "MetaCyc", 
                                    sample_table = tmp_meta2 
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
meco_pathway_metacyc_used$plot_bar(facet = "enterotype", bar_type = "notfull")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))
ggsave("output/5.pathway_metacyc_bar.enterotype.pdf",width = 20,height = 6)
tmp_plot_data <- meco_pathway_metacyc_used$data_abund %>% arrange(desc(all_mean_abund))
tmp_plot_data$Taxonomy <- factor(tmp_plot_data$Taxonomy,levels = rev(unique(tmp_plot_data$Taxonomy)))
tmp_plot_pathway <- unique(tmp_plot_data$Taxonomy)[1:15]

ggplot(tmp_plot_data %>% filter(Taxonomy %in% tmp_plot_pathway),
       aes(x=enterotype,y=Abundance,fill = Taxonomy)) +
  stat_summary(fun=median, geom="bar" ,width = 0.6,position = "stack")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(
    values = rev(c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02",
                   "#A6761D","#666666","#8DCEBB","#ECAF80","#BAB7D9","#F394C4",
                   "#B2D28E","#F2D480","#D2BA8E"))
  )
ggsave("output/5.pathway_metacyc_bar.enterotype.pdf",width = 8,height = 6)

tmp_data_stat <- data.frame()
for(i in 1:length(tmp_plot_pathway)){
  tmp_data <- tmp_plot_data %>% filter(Taxonomy == tmp_plot_pathway[i])
  tmp <- summary_stat(tmp_data,3,8) %>% mutate(Taxaonomy = tmp_plot_pathway[i]) 
  tmp_data_stat <- rbind(tmp_data_stat,tmp)
}
tmp_data_stat_short <- tmp_data_stat %>% 
  select(Group,Mean,Taxaonomy) %>% 
  spread(Group,Mean) %>% 
  mutate(Var = ET2 - ET1)
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
ggsave("output/5.pathway_metacyc_bar.enterotype.pdf",width = 11,height = 6)

#meco_pathway_metacyc_diff <- trans_diff$new(meco_pathway_metacyc,
#                                            method =  "lefse",group = "groups",
#                                           taxa_level = "pathway",
#                                           alpha = 0.05,p_adjust_method = "fdr"
#                                            #reference = c("Type,Health"),fixed_effects = c("Type","Sex")
#                                           )
#meco_pathway_metacyc_diff$plot_diff_bar(use_number = 1:20)
meco_pathway_metacyc$cal_betadiv(method = "bray")
meco_pathway_metacyc_beta <- trans_beta$new(dataset = meco_pathway_metacyc, 
                                            group = "group1", measure = "bray")

meco_pathway_metacyc_beta$cal_ordination(ordination = "PCoA")
meco_pathway_metacyc_beta_tmp <- meco_pathway_metacyc_beta$res_ordination$scores
meco_pathway_metacyc_beta_tmp2 <- trans_env$new(dataset = meco_pathway_metacyc, 
                                                add_data = meco_pathway_metacyc_beta_tmp[, 1:2])
# 'KW_dunn' for non-parametric test
meco_pathway_metacyc_beta_tmp2$cal_diff(group = "enterotype", 
                                        method = "anova")

meco_pathway_metacyc_beta$plot_ordination(plot_color = "group1", 
                                                plot_shape = "group1", 
                                                plot_type = c("point", "ellipse"))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_color_manual(
    values = color1[c(1,2)]
  )+
  scale_fill_manual(
    values = color1[c(1,2)]
  )
