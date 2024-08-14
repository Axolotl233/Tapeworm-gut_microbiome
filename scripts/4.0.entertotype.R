dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
  KLD <- function(x,y) sum(x *log(x/y))
  JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x<=0.0000001,pseudocount,x))
  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) { 
      resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                             as.vector(inMatrix[,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix) 
}
pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
  require(cluster)
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
  return(cluster)
}
pam.medoids=function(x,k) {
  require(cluster)
  medoids = as.vector(pam(as.dist(x), k, diss=TRUE)$medoids)
  return(medoids)
}

library(ade4)
#library(aPCoA)
library(cluster)
library(clusterSim)
library(reshape2)
library(vegan)
library(ggstatsplot)
library(class)
library(caret)
library(ape)
library(ggsankey)
library(Maaslin2)
set.seed(10001)
dir.create("output/entertype",showWarnings = F)
out_dir <- "output/entertype"

load("workfile/z.phyloseq.Rdata")
tmp_meta <- metadatadf %>% filter (groups %in% group[c(1:3)])
tmp_otu <- read.table(in_otu_ab,header = T,row.names = 1)
tmp_otu <- tmp_otu[,tmp_meta$sample]
identical(tmp_meta$sample,colnames(tmp_otu))
tmp_otu_genus <- tmp_otu[grep("g__",row.names(tmp_otu),perl = TRUE),]
tmp_otu_genus <- tmp_otu_genus[grep("s__",row.names(tmp_otu_genus),perl = TRUE,invert = T),]
row.names(tmp_otu_genus) <- gsub(".*\\|","",row.names(tmp_otu_genus))
tmp_otu_genus <- noise_removal(tmp_otu_genus,percent = 0.2,low = 0)
tmp_otu_genus <- noise_removal(tmp_otu_genus,low = 0.01,method = "mean_cut")
tmp_jsd = dist.JSD(tmp_otu_genus)
tmp_jsd_data <- as.data.frame(as.matrix(tmp_jsd))

tmp_train_sample <- tmp_meta[tmp_meta$groups %in% group[c(1,2)],2]
tmp_ex_sample <- tmp_meta[tmp_meta$group1 %in% groupl[c(3)],2]
tmp_meta_train <- tmp_meta[tmp_meta$sample %in% tmp_train_sample,]
tmp_jsd_train <- as.dist(as.matrix(tmp_jsd_data[tmp_train_sample,tmp_train_sample]))
tmp_otu_genus_train <- tmp_otu_genus[,tmp_train_sample]
tmp_nclusters=NULL
#k = 2
for (k in 1:20) { 
  if (k==1) {
    tmp_nclusters[k]=NA 
  } else {
    tmp_data_clusters=pam.clustering(tmp_jsd_train, k)
    tmp_nclusters[k]=index.G1(t(tmp_otu_genus_train),tmp_data_clusters,  d = tmp_jsd_train,
                              centrotypes = "medoids")
  }
}
tmp_nclusters[1] = 0
tmp_ncluster_d <- data.frame(n = 1:20,ch = tmp_nclusters)
ggplot(tmp_ncluster_d,aes(x = n,y = ch))+
  geom_point()+
  geom_line()+
  scale_x_continuous(breaks = seq(1,20,1))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x = "k clusters", y = "Calinski-Harabasz index", title = "Optimal number of clusters")
ggsave(paste(out_dir,"cross.opti_clusters_pam.pdf",sep = '/'),width = 6,height = 4)

k_best = which(tmp_nclusters == max(tmp_nclusters), arr.ind = TRUE)
tmp_cluster=pam.clustering(tmp_jsd_train, k = k_best)
tmp_medoids=pam.medoids(tmp_jsd_train, k = k_best)
tmp_silhouette=mean(silhouette(tmp_cluster, tmp_jsd_train)[,3])
cat(tmp_silhouette)

tmp_meta_train$enterotype <- tmp_cluster
tmp_meta_train$enterotype <- ifelse(tmp_meta_train$enterotype == 1,"ET1","ET2")
tmp_jsd_train_data <- as.data.frame(as.matrix(tmp_jsd_train))
identical(row.names(tmp_jsd_train_data),tmp_meta_train$sample)
tmp_jsd_train_data$class <- tmp_meta_train$enterotype

tmp_pca=dudi.pca(as.data.frame(t(tmp_otu_genus_train)), scannf=F, nf=10)
tmp_bet=bca(tmp_pca, fac=as.factor(tmp_cluster), scannf=F, nf=2) 
tmp_bet_res <- as.data.frame(t(tmp_bet$tab))
tmp_bet_res$taxon = row.names(tmp_bet_res)

tmp_et1 <- tmp_bet_res[which(tmp_bet_res[,1] == max(tmp_bet_res[,1])),3]
tmp_et2 <- tmp_bet_res[which(tmp_bet_res[,2] == max(tmp_bet_res[,2])),3]
res_et <- c(tmp_et1,tmp_et2) 

tmp_pcoa_res <- pcoa(tmp_jsd_train)
tmp_pcoa_res_matrix <- tmp_pcoa_res$vectors[,c(1:3)]
colnames(tmp_pcoa_res_matrix) <- c("dim1","dim2","dim3")
identical(tmp_meta_train$sample,row.names(tmp_pcoa_res_matrix))

tmp_train_res <- tmp_meta_train %>% mutate(PC1 = tmp_pcoa_res_matrix[,1],PC2 = tmp_pcoa_res_matrix[,2]) %>% 
  dplyr::select(PC1,PC2,enterotype,groups,sample)
tmp_train_res$enterotype <- as.character(tmp_train_res$enterotype)
write_csv(tmp_meta_train,paste(out_dir,"entertype_train.csv",sep = '/'))
tmp_adonis_res <- adonis2(tmp_jsd_train~enterotype+groups+age+sex+bmi, tmp_meta_train)

ggplot(data=tmp_train_res)+
  stat_ellipse(data=tmp_train_res,geom = "polygon",aes(x=PC1,y=PC2,fill=enterotype),alpha=0.05)+
  geom_point(aes(x=PC1,y=PC2,color=enterotype,shape=groups),size=2.5,alpha = 0.7)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(angle = 90,hjust = 0.5,vjust = 1),
        legend.position = c(0.1,0.1)
  )+
  labs(x=c("PCoA1"),
       y=c("PCoA2")
  )+
  scale_fill_manual(values = color1)+
  scale_color_manual(values = color1)+
  scale_x_continuous(breaks = seq(-4,4,2))+
  scale_y_continuous(breaks = seq(-4,4,2))

ggsave(paste(out_dir,"cross.pcoa.pdf",sep = '/'),height = 4.5,width = 4.5)
plot_l <- list()
for(i in 1:length(res_et)){
  tmp_d <- as.data.frame(t(tmp_otu_genus_train))
  tmp_d$sample = row.names(tmp_d)
  tmp_d <- tmp_d[,c("sample",res_et[i])]
  colnames(tmp_d) <- c("sample","value")
  tmp_plot <- inner_join(tmp_train_res,tmp_d,by = "sample")
  p1 <- ggplot(tmp_plot, aes(x=enterotype, y =value))+
    geom_boxplot(aes(color=enterotype),outlier.shape = NA)+
    geom_jitter(aes(color=enterotype),width = 0.15,height = 0,size = 3,alpha = 0.7)+
    stat_compare_means(method="wilcox.test",
                       comparisons = list(c("ET1","ET2")),
                       label = "p.format"
    )+
    theme_cowplot()+
    theme(legend.position="none", 
          axis.title.x = element_blank()
    )+
    scale_color_manual(
      values = color1[c(1,2)]
    )+
    labs(y = "Relative Abdunance",title = res_et[i])
  plot_l[[i]] <- p1
}
plot_grid(plot_l[[1]],plot_l[[2]],nrow = 1)
ggsave(paste(out_dir,"cross.entertotype.sp.pdf",sep = '/'),width = 5,height = 5)

Maaslin2(
  input_data = t(tmp_otu_genus_train), 
  input_metadata = tmp_meta_train, 
  min_prevalence = 0,
  normalization = "NONE",
  transform = "LOG",
  analysis_method = "LM" ,
  output = paste(out_dir,"cross.adjust",sep = '/'), 
  fixed_effects = c("enterotype","sex","age","bmi"),
  reference = c("enterotype,ET1")
)

ggstatsplot::ggbarstats(
  data = tmp_meta_train,
  x = groups,
  y = enterotype
)
ggsave(paste(out_dir,"cross.entertotype.ratio2.pdf",sep = '/'),width =6,height = 8)
tmp_sankey_p <- tmp_meta_train %>% dplyr::select(groups,enterotype) %>% make_long(enterotype,groups)
ggplot(tmp_sankey_p, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_alluvial(flow.alpha = 0.3,width = 0.2)+
  geom_alluvial_text(size = 3, color = "white")+
  theme_sankey()+
  theme(
    panel.grid = element_blank(),
    legend.title = element_blank())+
  scale_fill_manual(
    values = c(color1[c(1,2)],color0[c(1,2)])
  )+
  labs(x = "")
ggsave(paste(out_dir,"cross.entertotype.ratio.pdf",sep = '/'),width = 6,height = 5)

tmp_alpha <- read_csv("output/cross_base_analysis/cross.alpha_diversity_phyloseq.csv")
#tmp_alpha2 <- read_csv("output/cross_base_analysis/cross.alpha_diversity.csv")
#cor.test(tmp_alpha$Shannon,tmp_alpha2$shannon)
#cor.test(tmp_alpha$InvSimpson,tmp_alpha2$invsimpson)
#cor.test(tmp_alpha$Simpson,tmp_alpha2$simpson)

identical(tmp_alpha$sample,tmp_meta_train$sample)
tmp_alpha$enterotype <- tmp_meta_train$enterotype

tmp_list <- combine_list(tmp_alpha$enterotype)
tmp_alpha_p <- tmp_alpha %>% dplyr::select(sample,enterotype,Observed,Shannon,Simpson) %>%
  gather(key = "key",value = "value",Observed,Shannon,Simpson)
ggplot(tmp_alpha_p,aes(x=enterotype,y=value)) + 
  geom_boxplot(aes(color=enterotype),outlier.shape = NA)+
  geom_jitter(aes(color=enterotype),size = 3,height=0,width = 0.2,alpha = 0.7)+
  theme_cowplot()+
  theme(axis.title.x=element_blank(), 
        legend.position="none", axis.text.x = element_text(vjust = 0.5, hjust = 0.5))+
  scale_color_manual(values = color1)+
  labs(y="Observed value",title="Shannon index ")+
  stat_compare_means(
    comparisons = tmp_list,
    method = "wilcox.test", label = "p.format")+
  facet_wrap(.~key,scales = "free")
ggsave(paste(out_dir,"cross.entertotype2alpha.pdf",sep = '/'),width = 6,height = 5)
#facet_wrap(.~type)

tmp_alpha$class <- paste(tmp_alpha$type,tmp_alpha$enterotype,sep = "-")
tmp_list <- combine_list(tmp_alpha$class)
tmp_alpha_index <- c("Shannon","Simpson","Observed")
plot_l <- list()
for(i in 1:length(tmp_alpha_index)){
  tmp_data <- tmp_alpha[,c("type","class","enterotype",tmp_alpha_index[i])]
  colnames(tmp_data) <- c("type","class","enterotype","value")
  p <- ggplot(tmp_data,aes(x=class,y=value)) + 
    geom_boxplot(aes(color=enterotype),outlier.shape = NA)+
    geom_jitter(aes(color=enterotype,shape = type),
                position = position_jitterdodge(dodge.width = 0.75,jitter.width = 0.2,jitter.height = 0),
                size = 3,alpha = 0.7)+
    theme_cowplot()+
    theme(axis.title.x=element_blank(), 
          axis.text.x = element_text(vjust = 0.5, hjust = 0.5))+
    scale_color_manual(values = color1)+
    labs(y="Observed value",title=tmp_alpha_index[i])+
    stat_compare_means(
      comparisons = tmp_list,
      method = "wilcox.test", label = "p.format")
  plot_l[[i]] <- p
}
plot_grid(plot_l[[1]],plot_l[[2]],plot_l[[3]],nrow = 3)
ggsave(paste(out_dir,"cross.entertotype_class2alpha.pdf",sep = '/'),width = 5,height = 15)
ggsave(paste(out_dir,"cross.entertotype_class2simpson.pdf",sep = '/'),plot_l[[2]],width = 5,height = 5)
ggsave(paste(out_dir,"cross.entertotype_class2shannon.pdf",sep = '/'),plot_l[[1]],width = 5,height = 5)
ggsave(paste(out_dir,"cross.entertotype_class2observed.pdf",sep = '/'),plot_l[[3]],width = 5,height = 5)

rstatix::dunn_test(Shannon~class,data = tmp_alpha,p.adjust.method = "fdr")
rstatix::dunn_test(Simpson~class,data = tmp_alpha,p.adjust.method = "fdr")
rstatix::dunn_test(Observed~class,data = tmp_alpha,p.adjust.method = "fdr")

tmp_grid = expand.grid(.k = seq(2, 20, by = 1))
tmp_control = trainControl(method = "cv")

model <- caret::train(class~.,tmp_jsd_train_data,
                      method = 'knn',
                      trControl = tmp_control,
                      tuneGrid = tmp_grid)

plot.train(model)
tmp_best_k <- model$bestTune
tmp_t <- tmp_jsd_train_data$class
tmp_jsd_train_data %<>% dplyr::select(-class)
tmp_p <- predict(model,newdata = tmp_jsd_train_data)
confusionMatrix(table(tmp_p,tmp_t))

tmp_jsd_pre <- tmp_jsd_data[tmp_ex_sample,tmp_train_sample]
tmp_pp <- predict(model,newdata = tmp_jsd_pre)
tmp_pre_res <- data.frame(sample = row.names(tmp_jsd_pre),enterotype = tmp_pp)
tmp_meta_pre <- tmp_meta[tmp_meta$sample %in% tmp_pre_res$sample,]
tmp_meta_pre <- inner_join(tmp_meta_pre,tmp_pre_res,by= "sample")
tmp_res2 <- bind_rows(tmp_meta_pre,tmp_meta_train) %>% dplyr::select(sample,groups,group1,pair,enterotype)
row.names(tmp_res2) <- NULL
tmp_res2_pre <- tmp_res2 %>% filter(group1 == "after")

tmp_nclusters=NULL
for (k in 1:20) { 
  if (k==1) {
    tmp_nclusters[k]=NA 
  } else {
    tmp_data_clusters=pam.clustering(tmp_jsd_data, k)
    tmp_nclusters[k]=index.G1(t(tmp_otu_genus),tmp_data_clusters,  d = tmp_jsd_data,
                              centrotypes = "medoids")
  }
}
tmp_nclusters[1] = 0
tmp_ncluster_d <- data.frame(n = 1:20,ch = tmp_nclusters)
ggplot(tmp_ncluster_d,aes(x = n,y = ch))+
  geom_point()+
  geom_line()+
  scale_x_continuous(breaks = seq(1,20,1))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x = "k clusters", y = "Calinski-Harabasz index", title = "Optimal number of clusters")

k_best = which(tmp_nclusters == max(tmp_nclusters), arr.ind = TRUE)
tmp_cluster=pam.clustering(tmp_jsd, k = k_best)
tmp_medoids=pam.medoids(tmp_jsd, k = k_best)
tmp_silhouette=mean(silhouette(tmp_cluster, tmp_jsd)[,3])
cat(tmp_silhouette)

tmp_meta$enterotype <- tmp_cluster
tmp_meta$enterotype <- ifelse(tmp_meta$enterotype == 1,"ET1","ET2")
write_csv(tmp_meta,"workfile/z.meta.enterotype.csv")

tmp_jsd_data <- as.data.frame(as.matrix(tmp_jsd_data))
identical(row.names(tmp_jsd_data),tmp_meta$sample)
tmp_jsd_data$class <- tmp_meta$enterotype

tmp_pca=dudi.pca(as.data.frame(t(tmp_otu_genus)), scannf=F, nf=10)
tmp_bet=bca(tmp_pca, fac=as.factor(tmp_cluster), scannf=F, nf=2) 
tmp_bet_res <- as.data.frame(t(tmp_bet$tab))
tmp_bet_res$taxon = row.names(tmp_bet_res)

tmp_et1 <- tmp_bet_res[which(tmp_bet_res[,1] == max(tmp_bet_res[,1])),3]
tmp_et2 <- tmp_bet_res[which(tmp_bet_res[,2] == max(tmp_bet_res[,2])),3]
res_et <- c(tmp_et1,tmp_et2) 

tmp_meta_after <- tmp_meta %>% filter(group1 == "AF") %>% dplyr::select(sample,interval_class,enterotype) 
identical(tmp_meta_pre$sample,tmp_meta_after$sample)
identical(as.character(tmp_meta_pre$enterotype),tmp_meta_after$enterotype)
tmp_meta_pre$enterotype == tmp_meta_after$enterotype
#tmp_meta_after$enterotype2 <- tmp_meta_pre$enterotype
tmp_meta_after %>% group_by(enterotype) %>% summarise(n())

ggstatsplot::ggbarstats(
  data = tmp_res2 %>% filter (groups %in% group[1:2]),
  x = enterotype,
  y = groups
)+coord_flip()+scale_fill_manual(
  values = color1[c(2,1)]
)
ggsave(paste(out_dir,"cross.entertotype.ratio.bar.pdf",sep = '/'),width = 8,height = 3)

ggstatsplot::ggbarstats(
  data = tmp_res2 %>% filter (group1 %in% groupl[2:3]),
  x = enterotype,
  y = group1
)+coord_flip()+scale_fill_manual(
  values = color1[c(2,1)]
)
ggsave(paste(out_dir,"long.entertotype.ratio.bar.pdf",sep = '/'),width = 8,height = 3)

tmp_meta_2 <- read_table(in_meta_2) %>% filter(sp == "T_asiatica")
tmp_meta_2 <- as.data.frame(tmp_meta_2)
tmp_res <- data.frame(sample="",sample1="",sample2="",sample1_et="",sample2_et="",value=0,interval="")
for(i in 1:nrow(tmp_meta_2)){
  tmp_sam_b <- tmp_meta_2[i,7]
  tmp_sam_a <- tmp_meta_2[i,8]
  tmp_id <- tmp_meta_2[i,1]
  tmp_time <- tmp_meta_2[i,11]
  tmp_et_b <- tmp_res2[which(tmp_res2$sample == tmp_sam_b),5]
  tmp_et_a <- tmp_res2[which(tmp_res2$sample == tmp_sam_a),5]
  if(tmp_et_b == tmp_et_a){
    tmp_num = 0
  }else{
    if(tmp_et_b == "ET1"){
      tmp_num = 1
    }
    if(tmp_et_b == "ET2"){
      tmp_num = -1
    }
  }
  tmp_res[i,] = c(tmp_id,tmp_sam_b,tmp_sam_a,tmp_et_b,tmp_et_a,tmp_num,tmp_time)
}

tmp_res$value <- as.numeric(tmp_res$value)
tmp_res$interval <- as.numeric(tmp_res$interval)
tmp_res  %<>%
  mutate(
    class = case_when(
      value == 1 ~ "ET1_to_ET2",
      value == -1 ~ "ET2_to_ET1",
      TRUE ~ "no_change"
    ),
    class2 = case_when(
      interval <= 6 ~ "short",
      TRUE ~ "long"
      #interval <= 4 ~ "short",
      #interval > 4 & interval < 9  ~ "middle",
      #TRUE ~ "long"
    )
  )
tmp_res %>% group_by(class2,class) %>% summarise(n())
write_csv(x = tmp_res,paste(out_dir,"long.entertype_pre.csv",sep = '/'))

tmp_res$pos <- 1:nrow(tmp_res)
tmp_sample = tmp_res$sample
ggplot(tmp_res)+
  geom_point(aes(x = pos,y = value,color = sample1_et),size = 5)+
  scale_x_continuous(breaks = 1:length(tmp_sample),labels = tmp_sample)+
  scale_y_continuous(breaks = c(-1,0,1),labels = c("ET2_to_ET1","no_change","ET1_to_ET2"))+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    axis.title=element_blank(),
    axis.text.x = element_text(angle = 90,vjust = 0.5)
  )+
  scale_color_manual(
    values = color1[c(2,1)]
  )
ggsave(paste(out_dir,"long.entertotype_change.point.pdf",sep = '/'),width = 7,height = 3)

tmp_res %>% dplyr::select(sample,sample1_et,sample2_et) %>% gather(key = "class",value = "enterotype",-sample) %>%
  group_by(class,enterotype) %>% summarise(n())
tmp_sankey_p2 <- tmp_res %>% dplyr::select(sample,sample1_et,sample2_et) %>% 
  #gather(key = "class",value = "enterotype",-sample) %>% 
  make_long(sample1_et,sample2_et)

ggplot(tmp_sankey_p2, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  geom_alluvial(flow.alpha = 0.3,width = 0.2)+
  geom_alluvial_text(size = 3, color = "white")+
  theme_sankey()+
  theme(
    panel.grid = element_blank(),
    legend.title = element_blank())+
  scale_fill_manual(
    values = c(color1[c(1,2)],color0[c(1,2)])
  )+
  labs(x = "")
ggsave(paste(out_dir,"long.entertotype_change.sankey.pdf",sep = '/'),width = 6,height = 5)

tmp_res <- tmp_res %>% 
  mutate(
    class3 = case_when(
      class == "no_change" ~ "no_change",
      TRUE ~ "change"
    ))
ggstatsplot::ggbarstats(
  data = tmp_res,
  y = sample1_et,
  x = class3
)+coord_flip()+scale_fill_manual(
  values = color1[c(2,1)]
)

plot_l <- list()
for(i in 1:length(res_et)){
  tmp_data <- tmp_res %>% filter(class3 == "change") %>% dplyr::select(sample,sample1,sample2,class)# %>% 
  tmp_d <- as.data.frame(t(tmp_otu_genus))
  tmp_d$sample = row.names(tmp_d)
  tmp_d <- tmp_d[c(tmp_data$sample1,tmp_data$sample2),c("sample",res_et[i])]
  colnames(tmp_d) <- c("sample","ab")
  identical(tmp_d$sample,c(tmp_data$sample1,tmp_data$sample2))
  tmp_d$class <- rep(tmp_data$class,2)
  tmp_d <- inner_join(tmp_d,tmp_meta%>%dplyr::select(sample,group1,pair),by = "sample")
  tmp_d$group1 <- factor(tmp_d$group1,levels = (unique(tmp_d$group1)))
  tmp_list <- combine_list(tmp_d$group1)
  p <- ggplot(tmp_d,aes(x=group1,y= ab)) + 
    geom_point(aes(color=group1),size = 3)+
    geom_line(aes(group = pair),color="grey80",linetype="dashed")+
    facet_wrap(.~class)+
    theme_bw()+
    theme(legend.position="none", 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank()
    )+
    scale_color_manual(values = color1)+
    labs(y = "Relative Abdunance",title = res_et[i])#+
    #stat_compare_means(
    #  comparisons = tmp_list,paired = T,
    #  method = "wilcox.test", label = "p.format")
  plot_l[[i]] <- p
}
plot_grid(plot_l[[1]],plot_l[[2]],nrow = 2)
ggsave(paste(out_dir,"long.entertotype.sp_change2.pdf",sep = '/'),width = 5,height = 6)

tmp_data <- tmp_res %>% filter(class3 == "change") %>% dplyr::select(sample,sample1,sample2,class)
tmp_d <- as.data.frame(t(tmp_otu_genus))
tmp_d$sample = row.names(tmp_d)
tmp_d <- tmp_d[c(tmp_data$sample1,tmp_data$sample2),c("sample",res_et)]
identical(tmp_d$sample,c(tmp_data$sample1,tmp_data$sample2))
tmp_d$class <- rep(tmp_data$class,2)
tmp_d <- inner_join(tmp_d,tmp_meta%>%dplyr::select(sample,group1,pair),by = "sample")
tmp_d$group1 <- factor(tmp_d$group1,levels = (unique(tmp_d$group1)))
tmp_d_p <- tmp_d %>% gather(key = "sp",value = "ab",(colnames(tmp_d))[c(2,3)])
tmp_list <- combine_list(tmp_d$group1)
ggplot(tmp_d_p,aes(x=group1,y= ab)) + 
  geom_point(aes(color=group1),size = 3)+
  geom_line(aes(group = pair),color="grey80",linetype="dashed")+
  facet_grid(class~sp,scales = "free_y")+
  theme_cowplot()+
  theme(legend.position="none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank()
  )+
  scale_color_manual(values = color1)+
  labs(y = "Relative Abdunance")
ggsave(paste(out_dir,"long.entertotype.sp_change3.pdf",sep = '/'),width = 5,height = 6)

plot_et_change_sp <- function(x,y=tmp_res,z=color1[c(1,2)]){
  tmp_t <- y[grep(x,y$class),]
  tmp_t$paired <- 1:nrow(tmp_t)
  tmp_t <- tmp_t %>% 
    dplyr::select(sample1,sample2,paired) %>% 
    gather(key = "tmp",value = "sample",sample1,sample2)
  tmp_l <- list()
  for(i in 1:length(res_et)){
    tmp_d <- as.data.frame(t(tmp_otu_genus))
    tmp_d$sample = row.names(tmp_d)
    tmp_d <- tmp_d[,c("sample",res_et[i])]
    colnames(tmp_d) <- c("sample","value")
    tmp_plot <- inner_join(tmp_t,tmp_d,by = "sample")
    tmp_plot <- inner_join(tmp_plot,tmp_meta,by = "sample")
    tmp_plot$group1 <- factor(tmp_plot$group1,levels = c("BE","AF"))
    p1 <- ggplot(tmp_plot, aes(x=group1, y =value))+
      #geom_boxplot(aes(color=group2),outlier.shape = NA)+
      geom_point(aes(color=group1),size = 2.5)+
      geom_line(aes(group = paired),color="grey80",linetype="dashed")+
      theme_bw()+
      theme(legend.position="none", 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.title.x = element_blank()
      )+
      scale_color_manual(
        values = z
      )+
      labs(y = "Relative Abdunance",title = res_et[i])
    tmp_l[[i]] <- p1
  }
  return(tmp_l)
}
p1 <- plot_et_change_sp("ET1_to_ET2")
plot_l1 <- plot_et_change_sp("ET1_to_ET2",z=color1[c(1,2)])
plot_l2 <- plot_et_change_sp("ET2_to_ET1",z=color1[c(2,1)])
plot_grid(plot_l1[[1]],plot_l1[[2]],
          plot_l2[[1]],plot_l2[[2]],
          nrow = 2)
ggsave(paste(out_dir,"long.entertotype.sp_change.pdf",sep = '/'),width = 5,height = 6)

ggstatsplot::ggbarstats(
  data = tmp_res,# %>% filter (group1 == "health"),
  x = class,
  y = class2
)+coord_flip()+
  scale_fill_manual(
    values = c("#bebada","#fb8072","#80b1d3")
  )
ggsave("output/z.entertotype.bar_change3.pdf",width = 8,height = 3)

tmp_alpha_base <- read_csv("output/cross_base_analysis/cross.alpha_diversity_phyloseq.csv")
tmp_res_use <- tmp_res %>% dplyr::select(sample1,sample1_et,sample2_et,interval,class,class2) 
colnames(tmp_alpha_base) <- c("sample1","type","sex","age","bmi","observed","shannon",
                              "simpson","invsimpson")
colnames(tmp_res_use) <- c("sample1","sample1_et","sample2_et","interval","et_type","interval_class")
tmp_res_long_all <- inner_join(tmp_res_use,tmp_alpha_base, by = "sample1")

dist_pair_long <- read_csv("output/long_base_analysis/long.distance_BC_paired.csv")
colnames(dist_pair_long) <- c("sample1","sample2","dist_value","pair")
tmp_res_long_all <- inner_join(tmp_res_long_all,dist_pair_long, by = "sample1")
tmp_res_long_all <- tmp_res_long_all %>% 
  mutate(
    et_all = case_when(
      et_type == "no_change" ~ "no_change",
      TRUE ~ "change"
    ))

a <- ggplot(tmp_res_long_all,aes(x=et_all,y=interval)) + 
  geom_boxplot(aes(color=et_all),outlier.shape = NA)+
  geom_jitter(aes(color=et_all),size = 3,height=0,width = 0.3)+
  theme_cowplot()+
  theme(axis.title.x=element_blank(), 
        legend.position="none", axis.text.x = element_text(vjust = 0.5, hjust = 0.5))+
  scale_color_manual(
    values = c("#65387d","grey70")
  )+
  labs(y="Observed value",title="Interval Time")+
  stat_compare_means(
    comparisons = list(c("no_change","change")),
    method = "wilcox.test", label = "p.format")

b <- ggplot(tmp_res_long_all,aes(x=et_all,y=shannon)) + 
  geom_boxplot(aes(color=et_all),outlier.shape = NA)+
  geom_jitter(aes(color=et_all),size = 3,height=0,width = 0.2)+
  theme_cowplot()+
  theme(axis.title.x=element_blank(), 
        legend.position="none", axis.text.x = element_text(vjust = 0.5, hjust = 0.5))+
  scale_color_manual(
    values = c("#65387d","grey70")
  )+
  labs(y="Observed value",title="Shannon index ")+
  stat_compare_means(
    comparisons = list(c("no_change","change")),
    method = "wilcox.test", label = "p.format")

c <- ggplot(tmp_res_long_all, aes(x=et_all,y=dist_value)) + 
  geom_boxplot(aes(color=et_all),outlier.shape = NA)+
  geom_jitter(aes(color=et_all),size = 3,height=0,width = 0.3)+
  theme_cowplot()+
  theme(axis.title.x=element_blank(), 
        legend.position="none", axis.text.x = element_text(vjust = 0.5, hjust = 0.5))+
  scale_color_manual(
    values = c("#65387d","grey70")
  )+
  labs(y="Observed value",title="BC dist")+
  stat_compare_means(
    comparisons = list(c("no_change","change")),
    method = "wilcox.test", label = "p.format")

d <- ggplot(tmp_res_long_all, aes(x=et_all,y=observed)) + 
  geom_boxplot(aes(color=et_all),outlier.shape = NA)+
  geom_jitter(aes(color=et_all),size = 3,height=0,width = 0.3)+
  theme_cowplot()+
  theme(axis.title.x=element_blank(), 
        legend.position="none", axis.text.x = element_text(vjust = 0.5, hjust = 0.5))+
  scale_color_manual(
    values = c("#65387d","grey70")
  )+
  labs(y="Observed value",title="Richness")+
  stat_compare_means(
    comparisons = list(c("no_change","change")),
    method = "wilcox.test", label = "p.format")
e <- ggplot(tmp_res_long_all, aes(x=et_all,y=simpson)) + 
  geom_boxplot(aes(color=et_all),outlier.shape = NA)+
  geom_jitter(aes(color=et_all),size = 3,height=0,width = 0.3)+
  theme_cowplot()+
  theme(axis.title.x=element_blank(), 
        legend.position="none", axis.text.x = element_text(vjust = 0.5, hjust = 0.5))+
  scale_color_manual(
    values = c("#65387d","grey70")
  )+
  labs(y="Observed value",title="Simpson index")+
  stat_compare_means(
    comparisons = list(c("no_change","change")),
    method = "wilcox.test", label = "p.format")
plot_grid(a,c,d,b,e,nrow = 1)
ggsave(paste(out_dir,"long.entertotype.change.feature.pdf",sep = '/'),width = 12,height = 4)

tmp_data_m <- tmp_res_long_all %>% dplyr::select(interval,observed,shannon,simpson,dist_value)
row.names(tmp_data_m) <- tmp_res_long_all$pair
row.names(tmp_res_long_all) <- tmp_res_long_all$pair
Maaslin2(
  input_data = tmp_data_m, 
  input_metadata = tmp_res_long_all, 
  min_prevalence = 0,
  normalization = "NONE",
  transform = "LOG",
  analysis_method = "LM" ,
  output = paste(out_dir,"long.entertotype_change2feature_type",sep = '/'), 
  fixed_effects = c("et_type"),
  reference = c("et_type,no_change")
)
Maaslin2(
  input_data = tmp_data_m, 
  input_metadata = tmp_res_long_all, 
  min_prevalence = 0,
  normalization = "NONE",
  transform = "LOG",
  analysis_method = "LM" ,
  output = paste(out_dir,"long.entertotype_change2feature_all",sep = '/'), 
  fixed_effects = c("et_all"),
  reference = c("et_all,no_change")
)


#tmp_meta_entero <- tmp_meta_train %>% dplyr::select(sample,enterotype)
#tmp_meta_entero <- bind_rows(tmp_meta_entero,tmp_meta_after[,c(1,3)])

#tmp_meta_df <- inner_join(metadatadf,tmp_meta_entero,by = "sample")
#write_csv(tmp_meta_df,"workfile/z.meta_add_enterotype.csv")
