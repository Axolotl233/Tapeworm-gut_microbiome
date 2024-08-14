phase_ancom_diff_res <- function(res_tab_tmp,count_data,tmp_meta,gr1,gr2,group_col=8,sample_col=2){
  #res_tab_tmp <- ancombc2_res_tab
  ##count_data = tmp_otu
  #tmp_meta = tmp_meta
  #group_col=8
  #sample_col=2
  #gr1 = group[1]
  #gr2 = group[2]
  res_tab <- res_tab_tmp %>%
    mutate(group1_mean_ab=0,group2_mean_ab=0)%>%
    mutate(group1_nozero=0,group2_nozero=0,prevalence_rate=0)
  tmp_meta_g1 <- tmp_meta[tmp_meta[,group_col] == gr1,]
  tmp_meta_g2 <- tmp_meta[tmp_meta[,group_col] == gr2,]
  print(all(res_tab$taxon_id %in% row.names(count_data)))
  res_tab <- as.data.frame(res_tab)
  #i2 = 1
  for(i2 in 1:nrow(res_tab)){
    tmp_count <- as.data.frame(t(count_data[res_tab[i2,"taxon_id"],]))
    data_1 <- as.numeric(tmp_count[row.names(tmp_count) %in% tmp_meta_g1[,sample_col],])
    data_2 <- as.numeric(tmp_count[row.names(tmp_count) %in% tmp_meta_g2[,sample_col],])
    
    l_1 <- length(which(data_1 != 0))
    l_2 <- length(which(data_2 != 0))
    r_1 <- paste(l_1,length(data_1),sep = '|')
    r_2 <- paste(l_2,length(data_2),sep = '|')
    r_3 <- (l_1 + l_2)/(length(data_1)+length(data_2))
    group1_m <- round(mean(data_1),4)
    group2_m <- round(mean(data_2),4)
    
    #group1_m <- round((sum(data_1)/l_1),4)
    #group2_m <- round((sum(data_2)/l_2),4)
    
    res_tab[i2,ncol(res_tab)-4] <- as.numeric(group1_m)
    res_tab[i2,ncol(res_tab)-3] <- as.numeric(group2_m)
    res_tab[i2,ncol(res_tab)-2] <- r_1
    res_tab[i2,ncol(res_tab)-1] <- r_2
    res_tab[i2,ncol(res_tab)] <- r_3
  }
  res_tab$group1_mean_ab <- as.numeric(res_tab$group1_mean_ab)
  res_tab$group2_mean_ab <- as.numeric(res_tab$group2_mean_ab)
  res_tab$lfc_raw <- log2((res_tab$group2_mean_ab)/(res_tab$group1_mean_ab))
  return(res_tab)
}
phase_ancom_add_tax <- function(tab,tax,target_class = c("Phylum","Order"),ref_class = "Species"){
  tax_use = tax[tax$Species %in% tab$taxon_id,]
  tax_use = tax_use[,c(target_class,ref_class)]
  colnames(tax_use) = c(target_class,"taxon_id")
  res_tab <- merge(tab,tax_use,by.x = "taxon_id")
  return(res_tab)
}
phase_tax_file <- function(tax){
  tax[,c(1:6)] <- apply(tax[,c(1:6)], 2, function(x){gsub("Candidatus_","",x)})
  tax[,c(1:6)] <- apply(tax[,c(1:6)], 2, function(x){ifelse(grepl("unclassified",x),"Unclassified",x)})
  tax <- as.data.frame(tax)
  tax$Class <- ifelse(grepl("FGB",tax$Class),"Unclassified",tax$Class)
  tax$Order <- ifelse(grepl("FGB",tax$Order),"Unclassified",tax$Order)
  tax$Family <- ifelse(grepl("FGB",tax$Family),"Unclassified",tax$Family)
  tax$Genus <- ifelse(grepl("GGB",tax$Genus),"Unclassified",tax$Genus)
  return(tax)
}

load("workfile/z.phyloseq.Rdata")
load("workfile/z.diff_res.Rdata")
library(ANCOMBC)
dir.create("output/different_analysis",showWarnings = F)
out_dir <- "output/different_analysis"

phyloseqin_count_used <- subset_samples(phyloseqin_count,groups %in% group[c(1,2)])
tmp_meta <- metadatadf %>%
  filter(groups %in% group[c(1,2)])

tmp_tax <- as.data.frame(phyloseqin_count@tax_table)
tmp_tax <- phase_tax_file(tmp_tax)

### ancombc2 ancombc 2.20 infect vs hc
ancombc2_out = ancombc2(data =  phyloseqin_count_used,  p_adj_method = "fdr",
                        fix_formula = "groups+age+sex+bmi",group = "groups", tax_level = "Species",
                        prv_cut = 0.1, lib_cut = 1000, struc_zero = F, neg_lb = F,
                        iter_control = list(tol = 1e-5, max_iter = 20,verbose=T),
                        em_control = list(tol = 1e-5, max_iter = 100),
                        alpha = 0.05,global = F,verbose = T,n_cl = 8)

ancombc2_out_hc_vs_infect <- ancombc2_out
ancombc2_res_tab = ancombc2_out$res
ancombc2_res_tab <- ancombc2_res_tab[,c(1,seq(3,ncol(ancombc2_res_tab),5))]
colnames(ancombc2_res_tab) <-  c("taxon_id","lfc","se","w","pvalue","qvalue","diff")
#ancombc2_res_tab <- ancombc2_res_tab %>% dplyr::select(!"w-statistic")
ancombc2_res_tab$log_pvalue <- -log10(ancombc2_res_tab$pvalue)
ancombc2_res_tab$class <- "non"
ancombc2_res_tab <- ancombc2_res_tab %>%
  mutate(
    class = case_when(
      lfc > 0 & diff == "TRUE" ~ "up",
      lfc < 0 & diff == "TRUE" ~ "down",
      #lfc > 1 & pvalue < 0.05 ~ "up",
      #lfc < -1 & pvalue < 0.05 ~ "down",
      TRUE ~ "non"
    )
  )
#ancombc2_count_tab <- ancombc2_out$feature_table
#identical(rownames(ancombc2_count_tab),ancombc2_res_tab$taxon_id)
#ancombc2_res_tab$log_count <- log2(rowSums(ancombc2_count_tab))
tmp_otu <- as.data.frame(otu_table(phyloseqin))
rownames(tmp_otu) <- gsub('s__',"",rownames(tmp_otu))
tmp_otu <- tmp_otu[rownames(tmp_otu) %in% ancombc2_res_tab$taxon_id,]
tmp_otu <- tmp_otu[ancombc2_res_tab$taxon_id,tmp_meta$sample]
identical(row.names(tmp_otu),ancombc2_res_tab$taxon_id)

ancombc2_str_zero <- ancombc2_out$zero_ind
ancombc2_res_tab <- phase_ancom_diff_res(ancombc2_res_tab,tmp_otu,tmp_meta,group[1],group[2],group_col = 8)
ancombc2_res_tab <- phase_ancom_add_tax(tab = ancombc2_res_tab,tax = tmp_tax,target_class = c("Phylum","Genus"))

ancombc2_res_tab$taxon_class = ancombc2_res_tab[,ncol(ancombc2_res_tab)-1]
ancombc2_res_tab$taxon_class <- ifelse(ancombc2_res_tab$qvalue > 0.05,
                                       "ZZNone",ancombc2_res_tab$taxon_class)
ancombc2_res_tab %>% group_by(class) %>% summarise(n())
ancombc2_res_tab %>% filter (Phylum == "Actinobacteria") %>% group_by(class) %>% summarise(n())
View(ancombc2_res_tab %>% filter (Phylum == "Bacteroidetes") %>% group_by(class) %>% summarise(n()))
ancombc2_res_tab %>% group_by(class,Phylum) %>% summarise(n())
#tmp_data <- read_csv("output/different_analysis/1.csv") %>% dplyr::select(taxon_id,Gram)
#ancombc2_res_tab <- left_join(ancombc2_res_tab,tmp_data,by = "taxon_id" )
ancombc2_res_tab_hc_vs_infect <- ancombc2_res_tab
write_csv(ancombc2_res_tab_hc_vs_infect,paste(out_dir,"ancombc.hc_vs_infect.csv",sep = '/'))

ancombc2_res_tab_hc_vs_infect <- read_csv("output/different_analysis/ancombc.hc_vs_infect.csv")
ggplot()+
  geom_hline(aes(yintercept = -log10(0.05)) ,linetype = "dashed",color = "grey40")+
  geom_vline(aes(xintercept = 1) ,linetype = "dashed",color = "grey40")+
  geom_vline(aes(xintercept = -1) ,linetype = "dashed",color = "grey40")+
  geom_point(data = ancombc2_res_tab_hc_vs_infect %>% filter(class != "non" & qvalue <= 0.001),
             aes(x=lfc,y=-log10(qvalue),size=prevalence_rate,
                 color=class,shape = Gram),alpha = 0.8)+
  geom_point(data = ancombc2_res_tab_hc_vs_infect %>% filter(class != "non" & qvalue > 0.001),
             aes(x=lfc,y=-log10(qvalue),size=prevalence_rate,
                 color=class,shape = Gram),alpha = 0.4)+
  geom_point(data = ancombc2_res_tab_hc_vs_infect %>% filter(class == "non"),
             aes(x=lfc,y=-log10(qvalue),size=prevalence_rate,
                 color=class),alpha = 0.2)+
  scale_color_manual(values = c(color0[1], "grey", color0[2]))+
  scale_shape_manual(values = c(15,17,18))+
  labs(x="ANCOM LFC",y="-lg(P-value)",color = "Class",size="Prevalence_rate")+
  scale_x_continuous(breaks = seq(-7,7,1),limits = c(-8,7))+
  scale_y_continuous(breaks = seq (0,12,2))+
  theme_cowplot()+
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )+
  geom_text_repel(data = ancombc2_res_tab_hc_vs_infect %>% filter(class != "non" & qvalue <= 0.05),
                  aes(x=lfc,y=-log10(qvalue),label = taxon_id))

ggsave(paste(out_dir,"ancombc.hc_vs_infect.pdf",sep = '/'),width = 7.5,height = 5)

ggstatsplot::ggbarstats(
  data = ancombc2_res_tab_hc_vs_infect_gram <- ancombc2_res_tab_hc_vs_infect %>% filter(class != "non"),
  x = Gram,
  y = class
) + coord_flip()+theme_cowplot()+
  scale_fill_manual(
    values = c("#bcbddc","#fc9272")
  )
ggsave(paste(out_dir,"ancombc.hc_vs_infect_gram.pdf",sep = '/'),width = 10,height = 5)
ancombc2_res_tab_hc_vs_infect %>% filter(class != "non") %>%
  group_by(class,Gram) %>% summarise(n())### ancombc2 ancombc 2.20 long
phyloseqin_count_used <- subset_samples(phyloseqin_count,group1 %in% groupl[c(2,3)])
tmp_meta <- metadatadf %>%
  #filter(group %in% c("health","infected"))
  filter(group1 %in% groupl[c(2,3)])
tmp_otu <- as.data.frame(otu_table(phyloseqin))
tmp_tax <- as.data.frame(phyloseqin_count@tax_table)
tmp_tax <- phase_tax_file(tmp_tax)

ancombc2_out = ancombc2(data =  phyloseqin_count_used,  p_adj_method = "fdr",
                        fix_formula = "group1+interval",group = "group1",
                        tax_level = "Species",
                        rand_formula = "(1|pair)",
                        prv_cut = 0.1, lib_cut = 0, struc_zero = F, neg_lb = F,
                        iter_control = list(tol = 1e-2, max_iter = 20,verbose=T),
                        em_control = list(tol = 1e-5, max_iter = 100), 
                        alpha = 0.05,global = F,verbose = T,n_cl = 8)
ancombc2_out_long <- ancombc2_out
ancombc2_res_tab = ancombc2_out$res
ancombc2_res_tab <- ancombc2_res_tab[,c(1,seq(3,ncol(ancombc2_res_tab),3))]
colnames(ancombc2_res_tab) <-  c("taxon_id","lfc","se","w","pvalue","qvalue","diff")
#ancombc2_res_tab <- ancombc2_res_tab %>% dplyr::select(!"w-statistic")
ancombc2_res_tab$log_pvalue <- -log10(ancombc2_res_tab$pvalue)
ancombc2_res_tab$class <- "non"
ancombc2_res_tab <- ancombc2_res_tab %>%
  mutate(
    class = case_when(
      lfc > 0 & diff == "TRUE" ~ "up",
      lfc < 0 & diff == "TRUE" ~ "down",
      #lfc > 1 & pvalue < 0.05 ~ "up",
      #lfc < -1 & pvalue < 0.05 ~ "down",
      TRUE ~ "non"
    )
  )
tmp_otu <- as.data.frame(otu_table(phyloseqin))
rownames(tmp_otu) <- gsub('s__',"",rownames(tmp_otu))
tmp_otu <- tmp_otu[rownames(tmp_otu) %in% ancombc2_res_tab$taxon_id,]
tmp_otu <- tmp_otu[ancombc2_res_tab$taxon_id,tmp_meta$sample]
identical(row.names(tmp_otu),ancombc2_res_tab$taxon_id)

ancombc2_str_zero <- ancombc2_out$zero_ind
ancombc2_res_tab <- phase_ancom_diff_res(ancombc2_res_tab,tmp_otu,tmp_meta,groupl[2],groupl[3],group_col = 9)
ancombc2_res_tab <- phase_ancom_add_tax(tab = ancombc2_res_tab,tax = tmp_tax,target_class = c("Phylum","Genus"))
ancombc2_res_tab_long <- ancombc2_res_tab
write_csv(ancombc2_res_tab_long,paste(out_dir,"ancombc.long.csv",sep = '/'))
ancombc2_res_tab_long <- read_csv("output/different_analysis/ancombc.long.csv")
ggplot()+
  geom_point(data = ancombc2_res_tab_long %>% filter(class != "non" & qvalue <= 0.05),
             aes(x=lfc,y=-log10(qvalue),size=prevalence_rate,
                 color=class),alpha = 0.8)+
  geom_point(data = ancombc2_res_tab_long %>% filter(class != "non" & qvalue > 0.05),
             aes(x=lfc,y=-log10(qvalue),size=prevalence_rate,
                 color=class),alpha = 0.4)+
  geom_point(data = ancombc2_res_tab_long %>% filter(class == "non"),
             aes(x=lfc,y=-log10(qvalue),size=prevalence_rate,
                 color=class),alpha = 0.2)+
  scale_color_manual(values = c(color0[1], "grey", color0[2]))+
  scale_shape_manual(values = c(15,17,18))+
  labs(x="ANCOM LFC",y="-lg(P-value)",color = "Class",size="Prevalence_rate")+
  scale_x_continuous(breaks = seq(-8,8,2),limits = c(-8,8))+
  scale_y_continuous(breaks = seq (0,12,2))+
  theme_cowplot()+
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )+
  geom_text_repel(data = ancombc2_res_tab_long %>% filter(class != "non" & qvalue <= 0.05),
                  aes(x=lfc,y=-log10(qvalue),label = taxon_id))
ggsave(paste(out_dir,"ancombc.long.pdf",sep = '/'),width = 7.5,height = 5)

phyloseqin_count_used <- subset_samples(phyloseqin_count,enterotype %in% c("ET1"))
phyloseqin_count_used <- subset_samples(phyloseqin_count_used,groups %in% group[c(1,2)])
tmp_meta <- metadatadf %>% filter(sample %in% phyloseqin_count_used@sam_data$sample)
ancombc2_out = ancombc2(data =  phyloseqin_count_used,  p_adj_method = "fdr",
                        fix_formula = "groups+age+sex+bmi",group = "groups", tax_level = "Species",
                        prv_cut = 0.1, lib_cut = 100, struc_zero = T, neg_lb = F,
                        iter_control = list(tol = 1e-5, max_iter = 20,verbose=T),
                        em_control = list(tol = 1e-5, max_iter = 100), 
                        alpha = 0.05,global = F,verbose = T,n_cl = 8)
ancombc2_out_ET1 <- ancombc2_out
ancombc2_res_tab = ancombc2_out$res
ancombc2_res_tab <- ancombc2_res_tab[,c(1,seq(3,ncol(ancombc2_res_tab),5))]
colnames(ancombc2_res_tab) <-  c("taxon_id","lfc","se","w","pvalue","qvalue","diff")
#ancombc2_res_tab <- ancombc2_res_tab %>% dplyr::select(!"w-statistic")
ancombc2_res_tab$log_pvalue <- -log10(ancombc2_res_tab$pvalue)
ancombc2_res_tab$class <- "non"
ancombc2_res_tab <- ancombc2_res_tab %>%
  mutate(
    class = case_when(
      lfc > 0 & diff == "TRUE" ~ "up",
      lfc < 0 & diff == "TRUE" ~ "down",
      #lfc > 1 & pvalue < 0.05 ~ "up",
      #lfc < -1 & pvalue < 0.05 ~ "down",
      TRUE ~ "non"
    )
  )
#ancombc2_count_tab <- ancombc2_out$feature_table
#identical(rownames(ancombc2_count_tab),ancombc2_res_tab$taxon_id)
#ancombc2_res_tab$log_count <- log2(rowSums(ancombc2_count_tab))
tmp_otu <- as.data.frame(otu_table(phyloseqin))
rownames(tmp_otu) <- gsub('s__',"",rownames(tmp_otu))
tmp_otu <- tmp_otu[rownames(tmp_otu) %in% ancombc2_res_tab$taxon_id,]
tmp_otu <- tmp_otu[ancombc2_res_tab$taxon_id,tmp_meta$sample]
identical(row.names(tmp_otu),ancombc2_res_tab$taxon_id)

ancombc2_str_zero <- ancombc2_out$zero_ind
ancombc2_res_tab <- phase_ancom_diff_res(ancombc2_res_tab,tmp_otu,tmp_meta,group[1],group[2],group_col = 8)
ancombc2_res_tab <- phase_ancom_add_tax(tab = ancombc2_res_tab,tax = tmp_tax,target_class = c("Phylum","Genus"))
ancombc2_res_tab_ET1 <- ancombc2_res_tab
write_csv(ancombc2_res_tab_ET1,paste(out_dir,"ancombc.hc_vs_infect.ET1.csv",sep = '/'))

phyloseqin_count_used <- subset_samples(phyloseqin_count,enterotype %in% c("ET2"))
phyloseqin_count_used <- subset_samples(phyloseqin_count_used,groups %in% group[c(1,2)])
tmp_meta <- metadatadf %>% filter(sample %in% phyloseqin_count_used@sam_data$sample)
ancombc2_out = ancombc2(data =  phyloseqin_count_used,  p_adj_method = "fdr",
                        fix_formula = "groups+age+sex+bmi",group = "groups", tax_level = "Species",
                        prv_cut = 0.1, lib_cut = 100, struc_zero = T, neg_lb = F,
                        iter_control = list(tol = 1e-5, max_iter = 20,verbose=T),
                        em_control = list(tol = 1e-5, max_iter = 100), 
                        alpha = 0.05,global = F,verbose = T,n_cl=8)
ancombc2_out_ET2 <- ancombc2_out
ancombc2_res_tab = ancombc2_out$res
ancombc2_res_tab <- ancombc2_res_tab[,c(1,seq(3,ncol(ancombc2_res_tab),5))]
colnames(ancombc2_res_tab) <-  c("taxon_id","lfc","se","w","pvalue","qvalue","diff")
#ancombc2_res_tab <- ancombc2_res_tab %>% dplyr::select(!"w-statistic")
ancombc2_res_tab$log_pvalue <- -log10(ancombc2_res_tab$pvalue)
ancombc2_res_tab$class <- "non"
ancombc2_res_tab <- ancombc2_res_tab %>%
  mutate(
    class = case_when(
      lfc > 0 & diff == "TRUE" ~ "up",
      lfc < 0 & diff == "TRUE" ~ "down",
      #lfc > 1 & pvalue < 0.05 ~ "up",
      #lfc < -1 & pvalue < 0.05 ~ "down",
      TRUE ~ "non"
    )
  )
#ancombc2_count_tab <- ancombc2_out$feature_table
#identical(rownames(ancombc2_count_tab),ancombc2_res_tab$taxon_id)
#ancombc2_res_tab$log_count <- log2(rowSums(ancombc2_count_tab))
tmp_otu <- as.data.frame(otu_table(phyloseqin))
rownames(tmp_otu) <- gsub('s__',"",rownames(tmp_otu))
tmp_otu <- tmp_otu[rownames(tmp_otu) %in% ancombc2_res_tab$taxon_id,]
tmp_otu <- tmp_otu[ancombc2_res_tab$taxon_id,tmp_meta$sample]
identical(row.names(tmp_otu),ancombc2_res_tab$taxon_id)

ancombc2_str_zero <- ancombc2_out$zero_ind
ancombc2_res_tab <- phase_ancom_diff_res(ancombc2_res_tab,tmp_otu,tmp_meta,group[1],group[2],group_col = 8)
ancombc2_res_tab <- phase_ancom_add_tax(tab = ancombc2_res_tab,tax = tmp_tax,target_class = c("Phylum","Genus"))
ancombc2_res_tab_ET2 <- ancombc2_res_tab
write_csv(ancombc2_res_tab_ET2,paste(out_dir,"ancombc.hc_vs_infect.ET2.csv",sep = '/'))

save(ancombc2_res_tab_ET1,ancombc2_res_tab_ET2,
     ancombc2_res_tab_hc_vs_infect,ancombc2_res_tab_long,
     ancombc2_out_ET1,ancombc2_out_ET2,ancombc2_out_hc_vs_infect,ancombc2_out_long,
     file = "workfile/z.diff_res.Rdata")

#load("workfile/z.diff.Rdata")
rm(list = ls(pattern = "ancom"))
rm(list = ls(pattern = "data"))
rm(list = ls(pattern = "phyloseq"))
rm(list = ls(pattern = "tmp"))
