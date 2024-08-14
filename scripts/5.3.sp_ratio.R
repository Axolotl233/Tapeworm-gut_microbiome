load("workfile/z.phyloseq.Rdata")

phyloseqin_used <- subset_samples(phyloseqin,groups %in% group[c(1,2)])
phyloseqin_phylum <- tax_glom(phyloseqin_used,taxrank = "Phylum")
tmp_otu_p <- as.data.frame(otu_table(phyloseqin_phylum))
tmp_tax <- as.data.frame(tax_table(phyloseqin_phylum))
identical(row.names(tmp_otu_p),row.names(tmp_tax))
row.names(tmp_otu_p) <- tmp_tax$Phylum
tmp_tax <- c("Actinobacteria","Bacteroidetes","Firmicutes","Proteobacteria")
tmp_otu_p <- tmp_otu_p %>% mutate(taxon = row.names(tmp_otu_p)) %>%
  filter(taxon %in% tmp_tax) %>%
  select(-taxon)
tmp_otu_p <- as.data.frame(t(tmp_otu_p))
tmp_otu_p$sample <- row.names(tmp_otu_p)
tmp_meta <- metadatadf %>% select(sample,groups)
tmp_otu_pl <- inner_join(tmp_otu_p,tmp_meta,by = "sample") %>% 
  gather(key = "taxon",value = "value", tmp_tax)

plot_l = list() 
for(i in 1:length(tmp_tax)){
  tmp_data <- tmp_otu_pl %>% filter(taxon == tmp_tax[i])
  tmp_c=unique(tmp_data$groups)
  tmp_a <- as.data.frame(t(combn(tmp_c,2))) %>%
    mutate(tmp = str_c(V1,",",V2)) %>% 
    select(tmp)
  tmp_list <- (str_split(tmp_a$tmp,pattern = ","))
  p <- ggplot(tmp_data,aes(x=groups,y=value,color = groups))+
    geom_boxplot(outlier.shape = NA,)+
    geom_jitter(width = 0.2,height = 0)+
    theme_cowplot()+
    theme(
      legend.position = "none"
    )+
    labs(x="Group",y="",title = tmp_tax[i])+
    stat_compare_means(comparisons = tmp_list,
                       method = "wilcox.test", label = "p.format")+
    scale_color_manual(values = color0)
  plot_l[[i]] <- p
}
tmp_otu_p$BF_ratio <- tmp_otu_p$Bacteroidetes/tmp_otu_p$Firmicutes
tmp_otu_p2 <- tmp_otu_p %>% select(sample,BF_ratio)
tmp_otu_p2 <- inner_join(tmp_otu_p2,tmp_meta,by = "sample")
plot_grid(plot_l[[1]],plot_l[[2]],
          plot_l[[3]],plot_l[[4]],nrow = 2)
ggsave("workfile/phylum_all.pdf",width = 10,height = 10)
ggplot(tmp_otu_p2,aes(x=groups,y=BF_ratio,color = groups))+
  geom_boxplot(outlier.shape = NA,)+
  geom_jitter(width = 0.2,height = 0,size = 2)+
  theme_cowplot()+
  labs(x="Group",y="",title = "B/F_ratio")+
  stat_compare_means(comparisons = tmp_list,
                     method = "wilcox.test", label = "p.format")+
  scale_color_manual(values = color0)
ggsave("workfile/BF_ratio_all.pdf",width = 6,height = 6)

phyloseqin_genus <- tax_glom(phyloseqin_used,taxrank = "Genus")
tmp_otu_g <- as.data.frame(otu_table(phyloseqin_genus))
tmp_tax <- as.data.frame(tax_table(phyloseqin_genus))
identical(row.names(tmp_otu_g),row.names(tmp_tax))
row.names(tmp_otu_g) <- tmp_tax$Genus

tmp_tax <- c("Bacteroides","Prevotella")
tmp_otu_g <- tmp_otu_g %>% mutate(taxon = row.names(tmp_otu_g)) %>%
  filter(taxon %in% tmp_tax) %>%
  select(-taxon)
tmp_otu_g <- as.data.frame(t(tmp_otu_g))
tmp_otu_g$sample <- row.names(tmp_otu_g)
tmp_meta <- metadatadf %>% select(sample,groups)
tmp_otu_gl <- inner_join(tmp_otu_g,tmp_meta,by = "sample") %>% 
  gather(key = "taxon",value = "value", tmp_tax)

plot_l = list() 
for(i in 1:length(tmp_tax)){
  tmp_data <- tmp_otu_gl %>% filter(taxon == tmp_tax[i])
  tmp_c=unique(tmp_data$groups)
  tmp_a <- as.data.frame(t(combn(tmp_c,2))) %>%
    mutate(tmp = str_c(V1,",",V2)) %>% 
    select(tmp)
  tmp_list <- (str_split(tmp_a$tmp,pattern = ","))
  p <- ggplot(tmp_data,aes(x=groups,y=value,color = groups))+
    geom_boxplot(outlier.shape = NA,)+
    geom_jitter(width = 0.2,height = 0)+
    theme_cowplot()+
    theme(
      legend.position = "none"
    )+
    labs(x="Group",y="",title = tmp_tax[i])+
    stat_compare_means(comparisons = tmp_list,
                       method = "wilcox.test", label = "p.format")+
    scale_color_manual(values = color0)
  plot_l[[i]] <- p
}
plot_grid(plot_l[[1]],plot_l[[2]],nrow = 1)
ggsave("workfile/genus_cross.pdf",width = 6,height = 4)
tmp_otu_g$BP_ratio <- tmp_otu_g$Prevotella/tmp_otu_g$Bacteroides
tmp_otu_g2 <- tmp_otu_g %>% filter(Bacteroides != 0) %>% select(sample,BP_ratio)
tmp_otu_g2 <- inner_join(tmp_otu_g2,tmp_meta,by = "sample")
tmp_otu_g2$log_BP_ratio <- log10(tmp_otu_g2$BP_ratio+1)

ggplot(tmp_otu_g2,aes(x=groups,y=log_BP_ratio,color = groups))+
  geom_boxplot(outlier.shape = NA,)+
  geom_jitter(aes(shape = groups) ,width = 0.15,height = 0,size = 3,alpha = 0.7)+
  theme_cowplot()+
  labs(x="Group",y="Log10((Prevotella/Bacteroides)+1)",title = "P/B_ratio")+
  stat_compare_means(comparisons = tmp_list,
                     method = "wilcox.test", label = "p.format")+
  scale_color_manual(values = color0)
ggsave("workfile/PB_ratio_cross.pdf",width = 5,height = 5)

phyloseqin_used <- subset_samples(phyloseqin, group1 %in% c(groupl[2],groupl[3]))
tmp_meta <- metadatadf %>% dplyr::select(sample,pair,group1) %>% filter(group1 %in% groupl[c(2,3)]) 

phyloseqin_genus <- tax_glom(phyloseqin_used,taxrank = "Genus")
tmp_otu_g <- as.data.frame(otu_table(phyloseqin_genus))
tmp_tax <- as.data.frame(tax_table(phyloseqin_genus))
identical(row.names(tmp_otu_g),row.names(tmp_tax))
row.names(tmp_otu_g) <- tmp_tax$Genus

tmp_tax <- c("Bacteroides","Prevotella")
tmp_otu_g <- tmp_otu_g %>% mutate(taxon = row.names(tmp_otu_g)) %>%
  filter(taxon %in% tmp_tax) %>%
  select(-taxon)
tmp_otu_g <- as.data.frame(t(tmp_otu_g))
tmp_otu_g$sample <- row.names(tmp_otu_g)

tmp_otu_g <- inner_join(tmp_otu_g,tmp_meta,by = "sample")
zerp_sample <- unique(tmp_otu_g[which(tmp_otu_g$Bacteroides == 0),4])
tmp_otu_g <- tmp_otu_g %>% filter(! pair %in% zerp_sample)
tmp_otu_g$group1 <- factor(tmp_otu_g$group1,levels = groupl[c(2,3)])
tmp_otu_gl <- tmp_otu_g %>% gather(key = "taxon",value = "value", tmp_tax)
tmp_otu_gl$group1 <- factor(tmp_otu_gl$group1,levels = groupl[c(2,3)])
tmp_otu_g$BP_ratio <- tmp_otu_g$Prevotella/tmp_otu_g$Bacteroides
identical(tmp_otu_g$pair[1:(nrow(tmp_otu_g)/2)],tmp_otu_g$pair[((nrow(tmp_otu_g)/2)+1):nrow(tmp_otu_g)])
wilcox.test(tmp_otu_g$BP_ratio[1:(nrow(tmp_otu_g)/2)],tmp_otu_g$BP_ratio[((nrow(tmp_otu_g)/2)+1):nrow(tmp_otu_g)],paired = T)
wilcox.test(tmp_otu_g$Bacteroides[1:(nrow(tmp_otu_g)/2)],tmp_otu_g$Bacteroides[((nrow(tmp_otu_g)/2)+1):nrow(tmp_otu_g)],paired = T)
wilcox.test(tmp_otu_g$Prevotella[1:(nrow(tmp_otu_g)/2)],tmp_otu_g$Prevotella[((nrow(tmp_otu_g)/2)+1):nrow(tmp_otu_g)],paired = T)

plot_l = list()
i = 1 
for(i in 1:length(tmp_tax)){
  tmp_data <- tmp_otu_gl %>% filter(taxon == tmp_tax[i])
  tmp_c=unique(as.character(tmp_data$group1))
  tmp_a <- as.data.frame(t(combn(tmp_c,2))) %>%
    mutate(tmp = str_c(V1,",",V2)) %>% 
    select(tmp)
  tmp_list <- (str_split(tmp_a$tmp,pattern = ","))
  p <- ggplot(tmp_data,aes(x=group1,y=value))+
    geom_boxplot(aes(color = group1),outlier.shape = NA,)+
    geom_line(aes(group = pair),linetype = "dashed",color = "grey90")+
    geom_jitter(aes(color = group1),width = 0,height = 0,size = 2,alpha = 0.7)+
    theme_cowplot()+
    theme(
      legend.position = "none"
    )+
    labs(x="Group",y="",title = tmp_tax[i])+
    stat_compare_means(comparisons = tmp_list,paired = T,
                       method = "wilcox.test", label = "p.format")+
    scale_color_manual(values = color0[c(2,3)])
  plot_l[[i]] <- p
}
plot_grid(plot_l[[1]],plot_l[[2]],nrow = 1)
ggsave("workfile/genus_long.pdf",width = 6,height = 4)
tmp_otu_g$BP_ratio <- tmp_otu_g$Prevotella/tmp_otu_g$Bacteroides
tmp_otu_g$log_BP_ratio <- log10(tmp_otu_g$BP_ratio+1)

ggplot(tmp_otu_g,aes(x=group1,y=log_BP_ratio))+
  geom_line(aes(group = pair),linetype = "dashed",color = "grey90")+
  geom_boxplot(aes(color = group1),outlier.shape = NA,)+
  geom_jitter(aes(color = group1) ,width = 0,height = 0,size = 3,alpha = 0.7)+
  theme_cowplot()+
  labs(x="Group",y="Log10((Prevotella/Bacteroides)+1)",title = "P/B_ratio")+
  stat_compare_means(comparisons = tmp_list,paired = T,
                     method = "wilcox.test", label = "p.format")+
  scale_color_manual(values = color0[c(2,3)])
ggsave("workfile/PB_ratio_long.pdf",width = 6,height = 6)
