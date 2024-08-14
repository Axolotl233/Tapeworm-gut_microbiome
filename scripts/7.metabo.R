library(microbiomeSeq)
library(microbiome)

#tmp_sp <- c("s__GGB1237_SGB1623", "s__Gemmiger_SGB15292")
metaphlanToPhyloseq <- function(
    tax,
    metadat=NULL,
    simplenames=TRUE,
    roundtointeger=FALSE,
    split="|"){
  ## tax is a matrix or data.frame with the table of taxonomic abundances, rows are taxa, columns are samples
  ## metadat is an optional data.frame of specimen metadata, rows are samples, columns are variables
  ## if simplenames=TRUE, use only the most detailed level of taxa names in the final object
  ## if roundtointeger=TRUE, values will be rounded to the nearest integer
  xnames = rownames(tax)
  shortnames = gsub(paste0(".+\\", split), "", xnames)
  if(simplenames){
    rownames(tax) = shortnames
  }
  if(roundtointeger){
    tax = round(tax * 1e4)
  }
  x2 = strsplit(xnames, split=split, fixed=TRUE)
  taxmat = matrix(NA, ncol=max(sapply(x2, length)), nrow=length(x2))
  colnames(taxmat) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")[1:ncol(taxmat)]
  rownames(taxmat) = rownames(tax)
  for (i in 1:nrow(taxmat)){
    taxmat[i, 1:length(x2[[i]])] <- x2[[i]]
  }
  taxmat = gsub("[a-z]__", "", taxmat)
  taxmat = phyloseq::tax_table(taxmat)
  otutab = phyloseq::otu_table(tax, taxa_are_rows=TRUE)
  if(is.null(metadat)){
    res = phyloseq::phyloseq(taxmat, otutab)
  }else{
    res = phyloseq::phyloseq(taxmat, otutab, phyloseq::sample_data(metadat))
  }
  return(res)
}

data_m <- read_table(in_meta_1)
data_mb <- read.csv(in_metabolite,header = T,row.names = 1)
data_mb[is.na(data_mb)] <- 0
metadatadf <- as.data.frame(data_m) %>% filter(sp != "T_saginata") #%>% filter(!(sample %in% outlier_sample))
metadatadf  %<>%
  mutate(
    interval_class = case_when(interval <= 6 ~ "short",
                               interval > 6 ~ "long",
                               TRUE ~ "No_info"
                               #interval <= 4 ~ "short",
                               #interval > 4 & interval < 9  ~ "middle",
                               #TRUE ~ "long"
    )
  ) %>% filter(sample %in% colnames(data_mb)) 

row.names(metadatadf) <- metadatadf$sample
sampledata <- sample_data(metadatadf)
sampledata$groups <- factor(sampledata$groups,levels = group)
sampledata$group1 <- factor(sampledata$group1,levels = groupl)

mphlanin <- read.csv(in_otu_ab, sep = "\t", strip.white = T, stringsAsFactors = F,row.names=1)
mphlanin <- mphlanin[,colnames(mphlanin) %in% metadatadf$sample]
mphlanin <- noise_removal(mphlanin,method = "max_cut",percent = 0,low = 0)
which(rowSums(mphlanin) == 0)
mphlanin_used_tmp <- mphlanin
#mphlanin_used_tmp<-data.frame(apply(mphlanin, 2, function(x) as.numeric(as.character(x))))
#rownames(mphlanin_used_tmp) <- row.names(mphlanin)
# only Bacteria included
mphlanin_used_tmp2 <- mphlanin_used_tmp[(grepl("k__Bacteria", rownames(mphlanin_used_tmp),perl = TRUE)),]
#filtering # keep only the species level
mphlanin_used <- mphlanin_used_tmp2[grep("s__",rownames(mphlanin_used_tmp2),perl =TRUE,invert=FALSE),]
mphlanin_used <- mphlanin_used[grep("t__",rownames(mphlanin_used),perl =TRUE,invert=TRUE),]
#mphlanin_used <- mphlanin_used[-(unlist(lapply(as.list(tmp_sp), function(x) {grep(x,row.names(mphlanin_used))}))),]
mphlanin_used[is.na(mphlanin_used)] <- 0

all(sampledata$sample %in% colnames(mphlanin_used))
mphlanin_used <- mphlanin_used[,sampledata$sample]
identical(sort(rownames(sampledata)), sort(colnames(mphlanin_used)))
phyloseqin <- metaphlanToPhyloseq(mphlanin_used, metadat = sampledata)
phyloseqin_all <- metaphlanToPhyloseq(mphlanin_used_tmp2, metadat = sampledata)

### count
mphlanin_count <- read.csv(in_otu_count, sep = "\t", strip.white = T, stringsAsFactors = F,row.names=1)
mphlanin_count <- mphlanin_count[,-1]
mphlanin_count <- mphlanin_count[,colnames(mphlanin_count) %in% metadatadf$sample]
#mphlanin_count <- noise_removal(mphlanin_count,method = "max_cut",percent = 0,low = 0)
#which(rowSums(mphlanin_count) == 0)
mphlanin_count_used_tmp<-data.frame(apply(mphlanin_count, c(1,2), function(x) as.numeric(as.character(x))))
mphlanin_count[is.na(mphlanin_count)] <- 0
rownames(mphlanin_count_used_tmp) <- row.names(mphlanin_count)
length(which(rowSums(mphlanin_count_used_tmp) == 0))
mphlanin_count_used_tmp <- noise_removal(mphlanin_count_used_tmp,method = "max_cut",low = 0)
length(which(rowSums(mphlanin_count_used_tmp) == 0))
# only Bacteria included
mphlanin_count_used_tmp2 <- mphlanin_count_used_tmp[(grepl("k__Bacteria", rownames(mphlanin_count_used_tmp),perl = TRUE)),]
#filtering # keep only the species level
mphlanin_count_used <- mphlanin_count_used_tmp2[grep("s__",rownames(mphlanin_count_used_tmp2),perl =TRUE,invert=FALSE),]
mphlanin_count_used <- mphlanin_count_used[grep("t__",rownames(mphlanin_count_used),perl =TRUE,invert=TRUE),]
#mphlanin_count_used <- mphlanin_count_used[-(unlist(lapply(as.list(tmp_sp), function(x) {grep(x,row.names(mphlanin_count_used))}))),]
mphlanin_count_used[is.na(mphlanin_count_used)] <- 0

all(sampledata$sample %in% colnames(mphlanin_used))
mphlanin_count_used <- mphlanin_count_used[,sampledata$sample]
identical(sort(rownames(sampledata)), sort(colnames(mphlanin_count_used)))
phyloseqin_count <- metaphlanToPhyloseq(mphlanin_count_used, metadat = sampledata)
phyloseqin_count_all <- metaphlanToPhyloseq(mphlanin_count_used_tmp2, metadat = sampledata)

save(phyloseqin,phyloseqin_all,phyloseqin_count,phyloseqin_count_all,
     sampledata,metadatadf,data_mb,file="workfile/z.phyloseq_mb.Rdata")
rm(list = ls(pattern = "tmp"))
rm(list = ls(pattern = "mphlanin"))
rm(data_m)

library(vegan)
library(factoextra)
library(ropls)
library(MetaNet)
library(pcutils)
library(ggcorrplot)
library(gplots)
library(amap)
library(KEGGREST)
library(pheatmap)
#library(MetaboAnalystR)
library(Maaslin2)
library(circlize)
dir.create("output/metabo",showWarnings = F)
out_dir <- "output/metabo"

load("workfile/z.phyloseq_mb.Rdata")

data_mb <- read.csv(in_metabolite,header = T,row.names = 1)
data_mb[is.na(data_mb)] <- 0
data_mb <- data_mb[,metadatadf$sample]
data_mb_f <- noise_removal(data_mb,percent = 0.2,low = 0)
tmp_mb_class <- read_csv("input/metabolite_class.csv")

tmp_pathway <- as.data.frame(keggList("pathway"))
tmp_pathway$tmp <- row.names(tmp_pathway)
colnames(tmp_pathway) <- c("description","pathway")
tmp_link <- KEGGREST::keggLink("compound","pathway")
tmp_link <- data_frame(compound = as.character(tmp_link),pathway = names(tmp_link))
tmp_link$compound <- gsub(pattern = "cpd:","",tmp_link$compound)
tmp_link$pathway <- gsub(pattern = "path:","",tmp_link$pathway)
write_csv(x = tmp_link,"workfile/c2p_link.csv")
write_csv(x = tmp_pathway,"workfile/p2d_link.csv")
tmp_c <- unique(tmp_link$compound)
tmp_c2p <- data.frame(KEGG = "", pathway = "")
for(i in 1:length(tmp_c)){
  tmp1 <- tmp_link %>% filter(compound == tmp_c[i])
  tmp2 <- paste(tmp1$pathway,collapse = ";")
  tmp_c2p[i,] = c(tmp_c[i],tmp2)
}
tmp_mb_class <- left_join(tmp_mb_class,tmp_c2p,by = "KEGG")
write_csv(tmp_mb_class,paste(out_dir,"metabolite_class_mod.csv",sep = '/'))

tmp_data_pca <- as.data.frame(t(data_mb_f))
all(rownames(tmp_data_pca) == metadatadf$sample)
tmp_group <- factor(metadatadf$groups,levels = group[c(1,2)])
tmp_dist = vegdist(tmp_data_pca, method = "euclidean")
tmp_sex <- factor(metadatadf$sex,levels = c("female","male"))
tmp_age <- metadatadf$age
tmp_bmi <- metadatadf$bmi
tmp_site = data.frame(sample = rownames(tmp_data_pca),group = tmp_group,sex =tmp_sex,age = tmp_age,bmi = tmp_bmi)
tmp_adonis_result_dis = adonis2(tmp_data_pca~group+age+sex+bmi, tmp_site)
metadatadf %>% group_by(groups) %>% summarise(n())

#tmp_res_pca <- prcomp(tmp_data_pca, scale = T)
tmp_res_pca <- opls(tmp_data_pca)
tmp_res_pca_p <- as.data.frame(tmp_res_pca@scoreMN)
all(rownames(tmp_res_pca_p) == metadatadf$sample)
tmp_res_pca_p$group <- factor(metadatadf$groups,levels = group[c(1,2)])
tmp_res_pca_s <- as.data.frame(tmp_res_pca@summaryDF)

ggplot()+
  geom_point(data=tmp_res_pca_p,aes(x=p1,y=p2,color=group,shape=group),size=2.5)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_text(angle = 90,hjust = 0.5,vjust = 1)
  )+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  labs(x=c("PC1"),y=c("PC2"),title = "PCA",subtitle = paste("R2X(cum)=",tmp_res_pca_s$`R2X(cum)`,sep = " "))+
  stat_ellipse(data=tmp_res_pca_p,
               geom = "polygon",
               aes(x=p1,y=p2,color = group, fill=group),
               alpha=0.1)+
  scale_fill_manual(values = color0)+
  scale_color_manual(values = color0)
ggsave(paste(out_dir,"metabo.pca.pdf",sep = '/'),width = 5,height = 4)

tmp_res_pls <- opls(tmp_data_pca,tmp_group,permI = 100)
tmp_res_pls_s <- as.data.frame(tmp_res_pls@summaryDF)
tmp_res_pls_p <- as.data.frame(tmp_res_pls@scoreMN)
tmp_res_pls_perm <- as.data.frame(tmp_res_pls@suppLs$permMN)[,c(2,3,7)] %>% gather(key = "key",value = "value",-sim)
all(rownames(tmp_res_pls_p) == metadatadf$sample)

tmp_formula <- y ~ x
ggplot(tmp_res_pls_perm,aes(x = sim,y = value,color = key))+
  stat_smooth(method='lm',formula = tmp_formula,alpha = 0.15)+
  geom_point(size=3)+
  theme_bw()

tmp_res_pls_p$group <- factor(metadatadf$groups,levels = group[c(1,2)])
ggplot()+
  geom_point(data=tmp_res_pls_p,aes(x=p1,y=p2,color=group,shape=group),size=2.5)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_text(angle = 90,hjust = 0.5,vjust = 1)
  )+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  labs(x=c("PC1"),y=c("PC2")
  )+
  stat_ellipse(data=tmp_res_pls_p,
               geom = "polygon",
               aes(x=p1,y=p2,color = group, fill=group),
               alpha=0.1)+
  scale_fill_manual(values = color0)+
  scale_color_manual(values = color0)

tmp_res_opls <- opls(tmp_data_pca,tmp_group,predI = 1,orthoI = NA,permI = 100)
tmp_res_opls_s <- as.data.frame(tmp_res_opls@summaryDF)
tmp_res_opls_p <- as.data.frame(tmp_res_opls@scoreMN)
tmp_res_opls_p$o1 <- as.data.frame(tmp_res_opls@orthoScoreMN)[,1]
tmp_res_opls_perm <- as.data.frame(tmp_res_opls@suppLs$permMN)[,c(2,3,7)] %>% gather(key = "key",value = "value",-sim)
all(rownames(tmp_res_opls_p) == metadatadf$sample)

tmp_formula <- y ~ x
ggplot(tmp_res_opls_perm,aes(x = sim,y = value,color = key))+
  stat_smooth(method='lm',formula = tmp_formula,alpha = 0.15)+
  geom_point(size=3)+
  theme_bw()
ggsave(paste(out_dir,"metabo.diff.opls_model.pdf",sep = '/'),width = 6,height = 4)

tmp_res_opls_p$group <- factor(metadatadf$groups,levels = group[c(1,2)])
ggplot()+
  geom_point(data=tmp_res_opls_p,aes(x=p1,y=o1,color=group,shape=group),size=3)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_text(angle = 90,hjust = 0.5,vjust = 1)
  )+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  labs(x=c("PC1"),y=c("PCo1"),title = "OPLS-DA",
       subtitle = paste("R2X(cum)=",tmp_res_opls_s$`R2X(cum)`,"Q2(cum)=",tmp_res_opls_s$`Q2(cum)`,sep = " ") 
  )+
  stat_ellipse(data=tmp_res_opls_p,
               geom = "polygon",
               aes(x=p1,y=o1, fill=group),
               alpha=0.05)+
  scale_fill_manual(values = color0)+
  scale_color_manual(values = color0)
ggsave(paste(out_dir,"metabo.opls.pdf",sep = '/'),width = 5,height = 5)

tmp_res_opls_vipscore = as.data.frame(tmp_res_opls@vipVn)
colnames(tmp_res_opls_vipscore) = 'vip'
tmp_res_opls_vipscore$metabolites = rownames(tmp_res_opls_vipscore)
tmp_res_opls_vipscore = tmp_res_opls_vipscore[order(-tmp_res_opls_vipscore$vip),]
tmp_res_opls_vipscore$metabolites = factor(tmp_res_opls_vipscore$metabolites,
                                           levels = tmp_res_opls_vipscore$metabolites)

tmp_res_opls_loadingscore = as.data.frame(tmp_res_opls@loadingMN)
colnames(tmp_res_opls_loadingscore) <- "loadingscore"
tmp_res_opls_loadingscore$metabolites = rownames(tmp_res_opls_loadingscore)
tmp_res_opls_allscore = inner_join(tmp_res_opls_vipscore, 
                                   tmp_res_opls_loadingscore, by = 'metabolites')

tmp_data_mb <- as.data.frame(t(data_mb_f))
all(row.names(tmp_data_mb) == metadatadf$sample)
tmp_data_mb$groups <- metadatadf$groups
tmp_res <- data.frame(metabolites = "",gr1_mean="",gr2_mean="",all_mean="",lfc = 0,p_t.test = 0,p_wilcox.test = 0)
for(i in 1:(ncol(tmp_data_mb)-1)){
  tmp_data <- tmp_data_mb[,c(i,ncol(tmp_data_mb))]
  tmp_n <- colnames(tmp_data)[1]
  colnames(tmp_data) <- c("mb","groups")
  tmp_stat <- summary_stat(tmp_data)
  tmp_lfc <- log2(tmp_stat$Mean[2]/tmp_stat$Mean[1])  
  
  tmp0 <- mean(tmp_data$mb)
  tmp1 <- t.test(mb ~ groups,data = tmp_data)  
  tmp2 <- wilcox.test(mb ~ groups,data = tmp_data)
  
  tmp_res[i,] <- c(tmp_n,tmp_stat$Mean[1],tmp_stat$Mean[2],tmp0,tmp_lfc,tmp1$p.value,tmp2$p.value)
}
tmp_res[,c(2:7)] <- apply(tmp_res[,c(2:7)],2,as.numeric)
tmp_res <- inner_join(tmp_res,tmp_res_opls_allscore,by = "metabolites")
tmp_res$padj_t.test <- p.adjust(tmp_res$p_t.test,method = "BH")
tmp_res$padj_wilcox.test <- p.adjust(tmp_res$p_wilcox.test,method = "BH")

length(which(tmp_res$padj_wilcox.test < 0.05 & tmp_res$vip > 1))
tmp_res <- tmp_res[order(tmp_res$padj_wilcox.test),]

tmp_res <- tmp_res %>%
  mutate(
    diffclass = case_when(
      lfc > 0 & vip > 1 & padj_wilcox.test < 0.05 ~ "up",
      lfc < 0 & vip > 1 & padj_wilcox.test < 0.05 ~ "down",
      #lfc > 1 & pvalue < 0.05 ~ "up",
      #lfc < -1 & pvalue < 0.05 ~ "down",
      TRUE ~ "non"
    )
  )
tmp_mb_class <- read_csv(paste(out_dir,"metabolite_class_mod.csv",sep = '/'))
tmp_res <- inner_join(tmp_res,tmp_mb_class,by = "metabolites")
tmp_res[is.na(tmp_res)] <- "UnClass"
tmp_res %>% group_by(diffclass) %>% summarise(n())
write_csv(tmp_res,(paste(out_dir,"metabo.diff_raw.csv",sep = '/')))

#Maaslin2(
#  input_data = as.data.frame(t(data_mb_f)), 
#  input_metadata = metadatadf, 
#  min_prevalence = 0,
#  normalization = "NONE",
#  transform = "LOG",
#  analysis_method = "LM",
#  max_significance = 0.5,
#  output = "output/metabo.diff_adjust", 
#  fixed_effects = c("groups","age"),
#  reference = c("groups,health")
#)

ggplot()+
  geom_hline(aes(yintercept = -log10(0.05)) ,linetype = "dashed",color = "grey40")+
  geom_vline(aes(xintercept = 1) ,linetype = "dashed",color = "grey40")+
  geom_vline(aes(xintercept = -1) ,linetype = "dashed",color = "grey40")+
  geom_point(data = tmp_res %>% filter(diffclass != "non" & abs(lfc) > 1),
             aes(x=lfc,y=-log10(p_wilcox.test),size=vip,
                 color=diffclass,shape = SuperClass),alpha = 0.8)+
  geom_point(data = tmp_res %>% filter(diffclass != "non" & abs(lfc) <= 1),
             aes(x=lfc,y=-log10(p_wilcox.test),size=vip,
                 color=diffclass,shape = SuperClass),alpha = 0.4)+
  geom_point(data = tmp_res %>% filter(diffclass == "non"),
             aes(x=lfc,y=-log10(p_wilcox.test),size=vip,color=diffclass),shape = 16,alpha = 0.4)+
  scale_color_manual(
    #values=c("#DC0000B2", "grey", "#00A087B2"),
    values = c("#d73027", "grey", "#1a9850")
  )+
  scale_shape_manual(
    values = c(15,2,17,18,0,1,16,5,7,10),
    #values = c(0,1,2,5,6,7,10,11,8)
  )+
  labs(x="Log2FoldChange",y="-lg(P-value)",color = "Class",size="VIP")+
  scale_x_continuous(breaks = seq(-5,5,1))+
  scale_y_continuous(breaks = seq (0,10,1))+
  theme_bw()+
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )+
  geom_text_repel(data = tmp_res %>% filter(abs(lfc) > 1 & vip > 1 & p_wilcox.test < 0.01),
                  aes(x=lfc,y=-log10(p_wilcox.test),label = metabolites))

ggsave(paste(out_dir,"metabo.diff.pdf",sep = '/'),width = 8,height = 6)

#use https://www.metaboanalyst.ca/MetaboAnalyst/ for kegg analysis
diff_kegg_restab <- read_csv(paste(out_dir,"metabo_diff_pathway_kegg_all.csv",sep = '/'))
diff_kegg_restab$Enrichmentratio <-diff_kegg_restab$Hits/diff_kegg_restab$Expected
diff_kegg_restab$Metabolites <- factor(diff_kegg_restab$Metabolites,
                                       levels = rev(diff_kegg_restab$Metabolites))
ggplot(diff_kegg_restab %>% filter(P < 0.05),aes(x = Enrichmentratio, y = Metabolites))+
  geom_point(aes(size = Impact, color=LogP))+
  scale_size_continuous(range=c(3,6))+
  labs(color="Log10(P-value)",x="Enrichment Ratio",y="Metabolites",title="KEGG pathway")+
  theme_bw()+
  theme(axis.line = element_line(colour = "black"), 
        axis.text = element_text(color = "black",size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        axis.title = element_text(color = "black",size = 14))+
  scale_color_gradient(low = "navy", high = "firebrick3")
ggsave(paste(out_dir,"metabo.diff.kegg.pdf",sep = '/'),width = 8,height = 5)

diff_kegg_restab_mod <- read_csv(paste(out_dir,"metabo_diff_kegg.csv",sep = '/'))
diff_kegg_restab_mod$Enrichmentratio <-diff_kegg_restab_mod$Hits/diff_kegg_restab_mod$Expected
diff_kegg_restab_mod$Metabolites <- factor(diff_kegg_restab_mod$Metabolites,
                                           levels = rev(diff_kegg_restab_mod$Metabolites))
diff_kegg_restab_mod$LogP <- log10(diff_kegg_restab_mod$P)
ggplot(diff_kegg_restab_mod %>% filter(P < 0.05),aes(x = Enrichmentratio, y = Metabolites))+
  geom_point(aes(color=LogP),size = 5)+
  scale_size_continuous(range=c(3,6))+
  labs(color="Log10(P-value)",x="Enrichment Ratio",y="Metabolites Set",title="Metabolite Set Enrichment")+
  theme_bw()+
  theme(axis.line = element_line(colour = "black"), 
        axis.text = element_text(color = "black",size = 12),
        legend.text = element_text(size = 12),
        legend.title=element_text(size=12),
        axis.title = element_text(color = "black",size = 14))+
  scale_color_gradient(low = "navy", high = "firebrick3")
ggsave(paste(out_dir,"metabo.diff.MSEA.pdf",sep = '/'),width = 7,height = 5)

#diff_hmdb <- tmp_res_diff$HMDB
#diff_smpdb <- InitDataObjects("conc", "diff_smpdbora", FALSE)
#diff_smpdb <- Setup.MapData(diff_smpdb, diff_hmdb);
#diff_smpdb <- CrossReferencing(diff_smpdb, "hmdb");
#diff_smpdb <- CreateMappingResultTable(diff_smpdb)
#diff_smpdb <- SetMetabolomeFilter(diff_smpdb, F);
#diff_smpdb <- SetCurrentMsetLib(diff_smpdb, "smpdb_pathway", 2);
#diff_smpdb <- CalculateHyperScore(diff_smpdb)
#diff_smpdb_restab <- as.data.frame(diff_smpdb$analSet$ora.mat)
#diff_smpdb_restab$metabolites <- row.names(diff_smpdb_restab)
#diff_smpdb_restab$enrichmentratio <- diff_smpdb_restab$hits/diff_smpdb_restab$expected

#diff_kegg <- InitDataObjects("conc", "msetora", FALSE)
#diff_kegg <- Setup.MapData(diff_kegg, diff_hmdb);
#diff_kegg <- CrossReferencing(diff_kegg, "hmdb");
#diff_kegg <- CreateMappingResultTable(diff_kegg)
#diff_kegg <- SetMetabolomeFilter(diff_kegg, F);
#diff_kegg <- SetCurrentMsetLib(diff_kegg, "kegg_pathway", 2);
#diff_kegg <- CalculateHyperScore(diff_kegg)
#diff_kegg_restab <- as.data.frame(diff_kegg$analSet$ora.mat)
#diff_kegg_restab$metabolites <- row.names(diff_kegg_restab)
#diff_kegg_restab$enrichmentratio <-diff_kegg_restab$hits/diff_kegg_restab$expected

row.names(tmp_res) <- tmp_res$metabolites
tmp_res$log_all_mean <- log2(tmp_res$all_mean+1)
phyloseqin_used <- subset_samples(phyloseqin, groups %in% c(group[1],group[2]))
tmp_otu <- as.data.frame((otu_table(phyloseqin_used)))
row.names(tmp_otu) <- sub("s__","",row.names(tmp_otu))
tmp_tax <- as.data.frame((tax_table(phyloseqin_used)))
row.names(tmp_tax) <- sub("s__","",row.names(tmp_tax))

tmp_tax$Median <- (apply(tmp_otu,1,median))
tmp_tax$Mean <- (apply(tmp_otu,1,mean))
tmp_tax$Max <- (apply(tmp_otu,1,max))
tmp_tax$Log_Mean <- log2(tmp_tax$Mean+1)

tmp_highab_sp <- tmp_tax$Species[tmp_tax$Max > 0.05]
tmp_res_diff <- tmp_res %>% filter(diffclass != "non")

diff_metabo <- tmp_res_diff$metabolites
diff_sp <- read.table(paste(out_dir,"common_diff_sp.txt",sep = '/'))[,1]

data_mb_diff <- as.data.frame(t(data_mb[diff_metabo,]))
data_otu_diff <- as.data.frame(t(tmp_otu[diff_sp,]))
data_otu_diff <- data_otu_diff[,colnames(data_otu_diff)%in%tmp_highab_sp]

diff2diff_cor <- c_net_calculate(data_otu_diff,data_mb_diff,p.adjust.method = "BH")
diff2diff_cor$class <- diff2diff_cor$p.adjust
diff2diff_cor$class <- apply(diff2diff_cor$class,c(1,2),function(x){if(x > 0.05){x = ""}else if(x <= 0.05 & x > 0.01){x = "*"
}else{
  x = "**"
}})
pheatmap(diff2diff_cor$r)

heatmap.2(diff2diff_cor$r, 
          trace="none",
          dendrogram = c("both"),
          col=colorRampPalette(colors = c("navy","white","firebrick3"))(200),
          sepcolor="gray90",
          sepwidth=c(0.001,0.001),
          cellnote = diff2diff_cor$class,
          notecex = 2,
          notecol = "black",
          margins = c(15, 15)
)

tmp_net<- c_net_build(diff2diff_cor, r_thres = 0, p_thres = 0.05, delete_single = T)
#tmp_net <- c_net_annotate(tmp_net, tmp_tax["Phylum"],mode = "v")
tmp_net <- c_net_set(tmp_net,
                     tmp_tax[c("Phylum","Log_Mean","Median")],
                     tmp_res[c("log_all_mean","SuperClass")],
                     vertex_class = c("Phylum","SuperClass"),
                     vertex_size = c("Log_Mean","log_all_mean")
)
tmp_net_modu <- module_detect(tmp_net,n_node_in_module = 10,delete = T)
#tmp_net_modu <- filter_n_module(tmp_net_modu, n_node_in_module = 10,keep_id = 10)
tmp_net_modu <- c_net_set(tmp_net_modu,
                          vertex_class = "module"
                          
)
tmp_coor_mod <- g_layout_circlepack(tmp_net_modu, group = "module")
pdf(paste(out_dir,"metabo.network_cluster.pdf",sep = '/'),width = 12,height = 10)
c_net_plot(tmp_net_modu,coors = tmp_coor_mod,mark_module = T,vertex.label= get_v(tmp_net_modu)$name)
dev.off()

tmp_net_modu <- zp_analyse(tmp_net_modu)
#zp_plot(tmp_net_modu, mode = 1)
tmp_net_modu_v <- (get_v(tmp_net_modu))
tmp_net_modu_v$Class <- ifelse(tmp_net_modu_v$v_group == "v_group1","Microbiome","Metabilite")
tmp_lab <- c("Peripherals", "Network hubs", "Module hubs", "Connectors")
tmp_color <- c("#FCF6EFFC", "#ffffcc", "#e5d8bd", "#fddaec")
names(tmp_color) <- tmp_lab
tmp_back <- data.frame(x1 = c(0, 0.62, 0, 0.62), x2 = c(0.62, 1, 0.62, 1), 
                       y1 = c(-Inf, 2.5, 2.5, -Inf), y2 = c(2.5, Inf, Inf, 2.5), 
                       lab = factor(tmp_lab, levels = tmp_lab))
ggplot() + 
  geom_rect(data = tmp_back, mapping = aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2, fill = lab), alpha = 0.7)+
  guides(fill = guide_legend(title = "Topological roles"))+
  scale_fill_manual(values = tmp_color) + 
  geom_point(data = tmp_net_modu_v, 
             aes(x = Pi, y = Zi, color = factor(v_class),shape = Class),size = 2) + 
  scale_color_manual(values = setNames(unique(tmp_net_modu_v$color), unique(tmp_net_modu_v$v_class)))+
  theme_bw()+ 
  guides(colour = guide_legend(title = "Modules")) + 
  theme(strip.background = element_rect(fill = "white"),
        panel.grid = element_blank()) + 
  labs(x = "Participation Coefficient", y = "Within-module connectivity z-score")+
  ggrepel::geom_text_repel(data = tmp_net_modu_v %>% filter(roles != "Peripherals"),
                           aes(x = Pi, y = Zi, label = name),size = 4)
ggsave(paste(out_dir,"metabo.network_topo.pdf",sep = '/'),width = 7,height = 6)
tmp_net_modu_v %>% group_by(roles,v_group) %>% summarise(n())
#tmp_coor <- g_layout(tmp_net,
#                     group = "v_group", layout1 = as_polycircle() ,layout2 =in_circle(),
#                     )
tmp_coor <- g_layout_polyarc(tmp_net, group = "v_group")
c_net_plot(tmp_net,coors = tmp_coor,
           vertex.color = get_cols(15, "col1"),
           vertex.label= get_v(tmp_net)$name,
           vertex.label.cex = 0.5,
           vertex_size_range = c(3, 10),
           edge.width = (get_e(tmp_net)$weight),
           edge_width_range = c(0.5, 3),
           edge_legend_title = "Correlation", 
           edge_legend_order = c("positive", "negative"),
           group_legend_title= c("Bacteria","Metabolites"),
           size_legend = T
)
tmp_net_e_tab <- get_e(tmp_net)
tmp_net_v_tab <- get_v(tmp_net)
tmp_net_v_tab_sp <- tmp_net_v_tab %>% filter(v_group == "v_group1")
tmp_net_v_tab_mb <- tmp_net_v_tab %>% filter(v_group == "v_group2")

tmp_net_c_data <- tmp_net_e_tab %>% select(from,to,weight)

for(i in 1:nrow(tmp_net_c_data)){
  sp <- tmp_net_v_tab %>% filter(name == tmp_net_c_data[i,1])
  mb <- tmp_net_v_tab %>% filter(name == tmp_net_c_data[i,2])
  tmp_net_c_data[i,1] = sp$Phylum
  tmp_net_c_data[i,2] = mb$SuperClass
}
tmp_net_c_color <- NULL
tmp_net_c_color[unique(tmp_net_c_data$from)] <- color_bar_st[1:length(unique(tmp_net_c_data$from))]
tmp_net_c_color[unique(tmp_net_c_data$to)] <- color_bar_st[(length(unique(tmp_net_c_data$from))+1):(length(unique(tmp_net_c_data$from))+length(unique(tmp_net_c_data$to)))]
chordDiagram(tmp_net_c_data,
             grid.col = tmp_net_c_color,
             annotationTrack = "grid",
             annotationTrackHeight = c(0.04, 0.1),
             #transparency = 0.01,
             #link.lty = 0,    # 线路类型
             #link.border = 0
)
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "downward", niceFacing = T, cex = 0.8)
  circos.axis(h = "top", labels.cex = 0.4,labels.niceFacing = F, labels.pos.adjust =F)
}, bg.border = NA)
legend("right",pch=20,legend=names(tmp_net_c_color),
       col=tmp_net_c_color,bty="n",
       cex=1,pt.cex=3,border="black")

custom_sp <- read.table(paste(out_dir,"stachyose_degradation.metabo.txt",sep = '/'),head = F)
custom_sp$V1 <- gsub("s__","",custom_sp$V1)
row.names(custom_sp) <- custom_sp$V1

all2all_cor <- c_net_calculate(as.data.frame(t(tmp_otu)),as.data.frame(t(data_mb_f)),
                               p.adjust.method = "BH")
tmp_net_all <- c_net_build(all2all_cor, r_thres = 0, p_thres = 0.05, delete_single = T)
tmp_net_all <- c_net_set(tmp_net_all,
                         tmp_tax[c("Kingdom","Species","Log_Mean","Median")],
                         tmp_res[c("log_all_mean","SuperClass")],
                         vertex_class = c("Kingdom","SuperClass"),
                         vertex_size = c("Log_Mean","log_all_mean")
)
View(get_v(tmp_net_all))
#tmp_net_all <- c_net_filter(tmp_net_all, Species %in% c(custom_sp, NA))

all2all_cor_r <- as.data.frame(all2all_cor$r)
all2all_cor_padj <- as.data.frame(all2all_cor$p.adjust)
all2all_cor_p <- as.data.frame(all2all_cor$p.value)

all2all_cor_r_sp <- all2all_cor_r[rownames(all2all_cor_r) %in% custom_sp,]
all2all_cor_r_sp$sp <- row.names(all2all_cor_r_sp)
all2all_cor_r_sp <- gather(all2all_cor_r_sp,key = "mb",value = "cor",-sp)

all2all_cor_p_sp <- all2all_cor_p[rownames(all2all_cor_p) %in% custom_sp,]
all2all_cor_p_sp$sp <- row.names(all2all_cor_p_sp)
all2all_cor_p_sp <-  gather(all2all_cor_p_sp,key = "mb",value = "p",-sp)

all2all_cor_padj_sp <- all2all_cor_padj[rownames(all2all_cor_r) %in% custom_sp,]
all2all_cor_padj_sp$sp <- row.names(all2all_cor_padj_sp)
all2all_cor_padj_sp <-  gather(all2all_cor_padj_sp,key = "mb",value = "padj",-sp)

identical(all2all_cor_r_sp$sp,all2all_cor_padj_sp$sp)
identical(all2all_cor_r_sp$mb,all2all_cor_padj_sp$mb)
identical(all2all_cor_r_sp$sp,all2all_cor_p_sp$sp)

all2all_cor_r_sp$padj <-  all2all_cor_padj_sp$padj
all2all_cor_r_sp$p <-  all2all_cor_p_sp$p
write_csv(all2all_cor_r_sp,file = paste(out_dir,"meta_sp_cor.csv",sep = '/'))
#View(as.data.frame(diff2diff_cor$r)[rownames(diff2diff_cor$r %in% custom_sp),])

sta2all_cor_tmp <- c_net_calculate(as.data.frame(t(tmp_otu[row.names(tmp_otu)%in% custom_sp$V1,])),
                                   as.data.frame(t(data_mb_f)),
                                   p.adjust.method = "BH")
tmp_net_sta_tmp <- c_net_build(sta2all_cor_tmp, r_thres = 0, p_thres = 0.05, use_p_adj = F)
custom_mb <- (get_e(tmp_net_sta_tmp))
custom_mb <- custom_mb[custom_mb$from == "Bifidobacterium_longum","to"]
sta2all_cor_bif<- c_net_calculate(as.data.frame(t(tmp_otu[row.names(tmp_otu)%in% custom_sp$V1,])),
                                  as.data.frame(t(data_mb_f[row.names(data_mb_f)%in% custom_mb,])),
                                  p.adjust.method = "BH")
tmp_net_sta_bif <- c_net_build(sta2all_cor_bif, r_thres = 0, p_thres = 0.05, use_p_adj = F)
V(tmp_net_sta_bif)$new_attri <- c(rep("C",20),rep("B",24))
V(tmp_net_sta_bif)$new_attri[10] <- "A"
tmp_net_sta_bif <- c_net_set(tmp_net_sta_bif,
                             tmp_tax[c("Kingdom","Species","Log_Mean","Median")],
                             custom_sp[c("V2")],
                             tmp_res[c("log_all_mean","SuperClass")],
                             vertex_class = c("V2","SuperClass"),
                             #vertex_size = c("Log_Mean","log_all_mean")
                             vertex_group = "new_attri" 
)
plot(tmp_net_sta_bif)
tmp_net_sta_coor <- g_layout_polycircle(tmp_net_sta_bif, group = "new_attri")
tmp_net_sta_coor <- c_net_layout(tmp_net_sta_bif,method = as_circle_tree())
c_net_plot(tmp_net_sta_bif,
           #coors = tmp_net_sta_coor,
           vertex.color = get_cols(15, "col1"),
           vertex.label= get_v(tmp_net_sta_bif)$name,
           vertex.label.cex = 0.5,
           vertex_size_range = c(3, 10),
           #edge.width = (get_e(tmp_net_sta_bif)$weight),
           #edge_width_range = c(0.5, 3),
           #edge_legend_title = "Correlation", 
           edge_legend_order = c("positive", "negative"),
           group_legend_title= c("Bacteria","Metabolites")
           #size_legend = T
)
c_net_save(tmp_net_sta_bif, filename = paste(out_dir,"bifido_network",sep = '/'), format = "data.frame")
netD3plot(tmp_net_sta_bif)
