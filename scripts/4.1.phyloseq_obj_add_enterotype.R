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

data_m <- read_csv("workfile/z.meta.enterotype.csv")
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
  )

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
mphlanin_count_used_tmp<-data.frame(apply(mphlanin_count, 2, function(x) as.numeric(as.character(x))))
rownames(mphlanin_count_used_tmp) <- row.names(mphlanin_count)
mphlanin_count_used_tmp <- noise_removal(mphlanin_count_used_tmp,method = "max_cut",percent = 0,low = 0)
which(rowSums(mphlanin_count_used_tmp) == 0)
# only Bacteria included
mphlanin_count_used_tmp2 <- mphlanin_count_used_tmp[(grepl("k__Bacteria", rownames(mphlanin_count_used_tmp),perl = TRUE)),]
#filtering # keep only the species level
mphlanin_count_used <- mphlanin_count_used_tmp2[grep("s__",rownames(mphlanin_count_used_tmp2),perl =TRUE,invert=FALSE),]
mphlanin_count_used <- mphlanin_count_used[grep("t__",rownames(mphlanin_count_used),perl =TRUE,invert=TRUE),]
#mphlanin_count_used <- mphlanin_count_used[-(unlist(lapply(as.list(tmp_sp), function(x) {grep(x,row.names(mphlanin_count_used))}))),]
mphlanin_count_used[is.na(mphlanin_count_used)] <- 0

all(sampledata$sample %in% colnames(mphlanin_count_used))
mphlanin_count_used <- mphlanin_count_used[,sampledata$sample]
identical(sort(rownames(sampledata)), sort(colnames(mphlanin_count_used)))
phyloseqin_count <- metaphlanToPhyloseq(mphlanin_count_used, metadat = sampledata)
phyloseqin_count_all <- metaphlanToPhyloseq(mphlanin_count_used_tmp2, metadat = sampledata)
mphlanin_count_used <- noise_removal(mphlanin_count_used,percent = 0.1,low = 0,method = "pre_cut")
mphlanin_count_used <- mphlanin_count_used + 1
phyloseqin_count_ancombc <-  metaphlanToPhyloseq(mphlanin_count_used, metadat = sampledata)

save(phyloseqin,phyloseqin_all,phyloseqin_count,phyloseqin_count_all,sampledata,metadatadf,file="workfile/z.phyloseq.Rdata")
save(phyloseqin_count_ancombc,file="workfile/z.phyloseq.ancombc.Rdata")
rm(list = ls(pattern = "tmp"))
rm(list = ls(pattern = "mphlanin"))
rm(data_m)

