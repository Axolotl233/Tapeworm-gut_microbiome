load("workfile/z.phyloseq.Rdata")
diff_sp <- read_table("input/diff_sp_bp.txt")

sp1_s <- diff_sp$species
tmp_meta <- metadatadf %>% filter(groups %in% group[c(1,2)]) %>% select(sample,groups)
tmp_data_otu_all <- as.data.frame(phyloseqin@otu_table)
row.names(tmp_data_otu_all) <- gsub("s__","",row.names(tmp_data_otu_all))
tmp_data_sp <- tmp_data_otu_all[sp1_s,tmp_meta$sample]
tmp_data_sp <- noise_removal(tmp_data_sp,0.1,0.001)

tmp_data_sp <- as.data.frame(t(tmp_data_sp)) %>% mutate(sample=row.names(.))
tmp_data_sp <- inner_join(tmp_meta,tmp_data_sp,by="sample")
row.names(tmp_data_sp) <- tmp_data_sp$sample
tmp_data_sp_l <- tmp_data_sp
tmp_data_sp_l$groups <- as.character(tmp_data_sp_l$groups)
tmp_data_sp_l$groups <- factor(tmp_data_sp$groups,levels = c(group[c(1,2)]))
tmp_data_sp_l$sample <- factor(tmp_data_sp_l$sample,levels = tmp_data_sp$sample)
tmp_data_sp_l %<>%
  gather(key="key",value="value",-colnames(tmp_data_sp_l)[c(1,2)])
tmp_data_sp_l$value <- tmp_data_sp_l$value + 1e-5
tmp_data_sp_l$log_value <- log10(tmp_data_sp_l$value)

ggplot(tmp_data_sp_l)+
  geom_boxplot(aes(x = key, y = log_value, color = groups,fill = groups),outlier.shape = NA,alpha = 0.7)+
  theme_bw()+
  theme(
    axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid = element_blank()
  )+
  scale_color_manual(values = c(color0[c(1,2)]))+
  scale_fill_manual(values = c(color0[c(1,2)]))+
  coord_flip()

ggsave("output/different_analysis/ancombc.hc_vs_infect.bp.ab.pdf",width = 5,height = 5)

tmp_meta <- metadatadf %>% filter(group1 %in% groupl[c(2,3)]) %>% select(sample,group1,pair)
tmp_meta2 <- read_table(in_meta_2) %>% filter(sp != "T_saginata")
tmp_data_sp <- tmp_data_otu_all[sp1_s,tmp_meta$sample]
tmp_data_sp <- noise_removal(tmp_data_sp,0.1,0.001)

tmp_data_sp <- as.data.frame(t(tmp_data_sp)) %>% mutate(sample=row.names(.))
tmp_data_sp <- inner_join(tmp_meta,tmp_data_sp,by="sample")
row.names(tmp_data_sp) <- tmp_data_sp$sample
tmp_data_sp_l <- tmp_data_sp
tmp_data_sp_l$group1 <- as.character(tmp_data_sp_l$group1)
tmp_data_sp_l$group1 <- factor(tmp_data_sp$group1,levels = c(groupl[c(2,3)]))
tmp_data_sp_l$sample <- factor(tmp_data_sp_l$sample,levels = tmp_data_sp$sample)
tmp_data_sp_l %<>%
  gather(key="key",value="value",-colnames(tmp_data_sp_l)[c(1,2,3)])
tmp_data_sp_l$value <- tmp_data_sp_l$value + 1e-5
tmp_data_sp_l$log_value <- log10(tmp_data_sp_l$value)
tmp_list = combine_list(tmp_data_sp_l$group1)

ggplot(tmp_data_sp_l,aes(x = group1, y = log_value, color = group1))+
  geom_boxplot(outlier.shape = NA,alpha = 0.7)+
  geom_line(aes(group = pair),linetype = "dashed",color = "grey80")+
  geom_point(size = 2.5,alpha = 0.7)+
  facet_wrap(.~key,nrow = 1)+
  theme_cowplot()+
  theme(
    axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid = element_blank()
  )+
  scale_color_manual(values = c(color0[c(2,3)]))+
  scale_fill_manual(values = c(color0[c(2,3)]))+
  stat_compare_means(method="wilcox.test",
                     comparisons = tmp_list,
                     label = "p.format",paired = T)
  

ggsave("output/different_analysis/ancombc.hc_vs_infect.bp.ab.long.pdf",width = 15,height = 5)

