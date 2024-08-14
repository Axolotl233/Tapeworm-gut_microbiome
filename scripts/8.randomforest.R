library(randomForest)
library(rfUtilities)
library(rfPermute)
library(caret)
library(pROC)

set.seed(10001)
dir.create("output/randomforest",showWarnings = F)
out_dir <- "output/randomforest/"
load("workfile/z.phyloseq.Rdata")

phyloseqin_used <- subset_samples(phyloseqin,groups %in% group[c(1,2)])
tmp_otu <- as.data.frame(otu_table(phyloseqin_used))
tmp_meta <- metadatadf %>% dplyr::select(sample,groups)
tmp_meta <- tmp_meta[tmp_meta$sample %in% colnames(tmp_otu), ]
identical(tmp_meta$sample,colnames(tmp_otu))
tmp_otu <- noise_removal(tmp_otu,percent=0.1,low = 0, method = "pre_cut")
tmp_otu <- sweep(tmp_otu,2,colSums(tmp_otu),"/") * 100
tmp_data <- as.data.frame(t(tmp_otu))
tmp_data$sample <- row.names(tmp_data)
tmp_data <- inner_join(tmp_data,tmp_meta,by="sample") 
row.names(tmp_data) <- tmp_data$sample
tmp_data <- tmp_data %>% select(!sample) 

tmp_train_use <- sample(nrow(tmp_data), nrow(tmp_data)*0.7)
tmp_data_train <- tmp_data[tmp_train_use,]
#tmp_data_pre <- tmp_data[! row.names(tmp_data) %in% row.names(tmp_data_train),]

tmp_errrate <- c(1)
for(i in 1:ncol(tmp_data_train)-1){
  tmp_model <- randomForest(x=tmp_data_train[,1:(ncol(tmp_data_train)-1)] , 
                    y=factor(tmp_data_train$groups,levels = group[c(1,2)]),
                    ntree=500, 
                    importance=TRUE, 
                    proximity=TRUE, 
                    mtry=i,
                    na.action=na.omit)
  tmp_err<-mean(tmp_model$err.rate)
  tmp_errrate[i] <- mean(tmp_err)
}
which.min(tmp_errrate)
model_classify <- randomForest(x=tmp_data_train[,1:(ncol(tmp_data_train)-1)] , 
                               y=factor(tmp_data_train$groups,levels = group[c(1,2)]),
                               ntree=1000, 
                               importance=TRUE, 
                               proximity=TRUE, 
                               replace = TRUE,
                               mtry=which.min(tmp_errrate),
                               na.action=na.omit
                               )
plot(model_classify)
model_classify2 <- randomForest(x=tmp_data_train[,1:(ncol(tmp_data_train)-1)] , 
                               y=factor(tmp_data_train$groups,levels = group[c(1,2)]),
                               ntree=1000, 
                               importance=TRUE, 
                               proximity =TRUE, 
                               replace = TRUE,
                               mtry=which.min(tmp_errrate),
                               na.action=na.omit,
                               oob_score=TRUE
)
model_imp <- importance(model_classify2)
model_imp <- data.frame(predictors = rownames(model_imp), model_imp)
model_imp_sort <- arrange(model_imp, desc(MeanDecreaseAccuracy))
#model_imp_sort_20 <- model_imp_sort[1:45,]
#model_imp_sort_20$predictors <- gsub("s__","",model_imp_sort_20$predictors)
model_classify_pm <- rfPermute(x=tmp_data_train[,1:(ncol(tmp_data_train)-1)] , 
                               y=factor(tmp_data_train$groups,levels = group[c(1,2)]),
                               ntree=1000, 
                               importance=TRUE, 
                               proximity =TRUE, 
                               replace = TRUE,
                               mtry=which.min(tmp_errrate),
                               na.action=na.omit,
                               oob_score=TRUE
)
model_imp_pvalue <- as.data.frame(model_classify_pm$pval[,,2]) %>% select(MeanDecreaseAccuracy)
model_imp_pvalue$predictors <- row.names(model_imp_pvalue)
colnames(model_imp_pvalue) <- c("Pvalue","predictors")
model_imp_sort <- inner_join(model_imp_sort,model_imp_pvalue,by = "predictors")
model_imp_sort <- model_imp_sort %>% mutate(
  sig = case_when(
    Pvalue < 0.01 ~ "**",
    Pvalue < 0.05 & Pvalue >= 0.01 ~ "*",
    TRUE ~ "ns"
  )
)
write_csv(model_imp_sort,file = paste(out_dir,"model_imp.csv",sep = '/'))

model_imp_sort_top <- model_imp_sort[model_imp_sort$MeanDecreaseAccuracy>2 & model_imp_sort$Pvalue<0.05,]
model_imp_sort_top$predictors <- factor(model_imp_sort_top$predictors,levels = rev(model_imp_sort_top$predictors))

tmp_data_train_rfcv <- list()
for(i in 1:5){
  tmp_data_train_rfcv[[i]] <- randomForest::rfcv(tmp_data_train[,1:(ncol(tmp_data_train)-1)] , 
                                          factor(tmp_data_train$groups,levels = group[c(1,2)]),
                                          cv.fold = 10, 
                                          step = 1.2,
                                          ntree=1000,
                                          na.action=na.omit, 
                                          oob_score=TRUE
  )
  print(i)
}
tmp_data_train_rfcv_cv <- data.frame(sapply(tmp_data_train_rfcv, '[[', 'error.cv'))
tmp_data_train_rfcv_cv$otu_num <- row.names(tmp_data_train_rfcv_cv)
tmp_data_train_rfcv_cv <- gather(tmp_data_train_rfcv_cv,key = "sp_num",value = "cv",!otu_num)
tmp_data_train_rfcv_cv_stat <- summary_stat(tmp_data_train_rfcv_cv,value_col = 3,class_col = 1) %>% arrange(Mean)
tmp_data_train_rfcv_cv_stat$Group <- as.numeric(tmp_data_train_rfcv_cv_stat$Group)
ggplot(tmp_data_train_rfcv_cv_stat)+
  geom_point(aes(x = Group, y = Mean))+
  geom_line(aes(x = Group, y = Mean),linetype = "dashed")+
  geom_errorbar(aes(x = Group, y = Mean,ymin = Mean-Se,ymax = Mean+Se))+
  scale_x_continuous(limits = c(0,100),
                     breaks=tmp_data_train_rfcv_cv_stat$Group
                     )+
  theme_bw()+
  theme(
    panel.grid = element_blank()
  )+
  labs(x = "Sp num", y = "Mean Cv")
ggsave(paste(out_dir,"randomforest.error_cv.pdf",sep = '/'))

ggplot(model_imp_sort_top %>% head(28)) +
  geom_bar(aes(x = predictors, y = MeanDecreaseAccuracy),stat = "identity", fill = "grey80") +
  geom_text(aes(x = predictors, y = MeanDecreaseAccuracy+0.25,label = sig),angle = 270)+
  coord_flip() +
  labs(title= "The important species",x = "Species")+
  theme(strip.background = element_blank())+
  theme_bw()
ggsave(paste(out_dir,"cross.rf_imp.pdf",sep = '/'),width = 8,height = 8)

plot_l1 <- list()
for(i in 1:nrow(model_imp_sort_top)){
  tmp_t <- as.data.frame(t(tmp_otu))
  tmp_t$sample <- row.names(tmp_t)
  tmp_d <- inner_join(tmp_meta,tmp_t,by ="sample") %>% dplyr::select(model_imp_sort_top[i,1],sample,groups)
  
  tmp_c=unique(tmp_d$groups)
  tmp_a <- as.data.frame(t(combn(tmp_c,2))) %>%
    mutate(tmp = str_c(V1,",",V2)) %>% 
    select(tmp)
  tmp_list <- (str_split(tmp_a$tmp,pattern = ","))
  
  colnames(tmp_d) <- c("value","sample","groups")
  p <- ggplot(tmp_d,aes(x=groups,y=value,color = groups))+
    geom_boxplot(outlier.shape = NA,)+
    geom_jitter(width = 0.2,height = 0)+
    theme_cowplot()+
    theme(
      legend.position = "none"
    )+
    labs(x="Group",y="",title = model_imp_sort_top[i,1])+
    stat_compare_means(comparisons = tmp_list,
                       method = "wilcox.test", label = "p.format")+
    scale_color_manual(values = color0)
  plot_l1[[i]] <- p
}
plot_grid(plot_l1[[1]],plot_l1[[2]],plot_l1[[3]],plot_l1[[4]],
          plot_l1[[5]],plot_l1[[6]],plot_l1[[7]],plot_l1[[8]],
          plot_l1[[9]],plot_l1[[10]],plot_l1[[11]],plot_l1[[12]],
          plot_l1[[13]],plot_l1[[14]],plot_l1[[15]],plot_l1[[16]],
          plot_l1[[17]],plot_l1[[18]],plot_l1[[19]],plot_l1[[20]],
          nrow = 4)
ggsave(paste(out_dir,"cross.rf_imp_sp20.pdf",sep = '/'),width = 20,height = 16)
plot_grid(plot_l1[[21]],plot_l1[[22]],plot_l1[[23]],plot_l1[[24]],
          plot_l1[[25]],plot_l1[[26]],plot_l1[[27]],plot_l1[[28]],
          plot_l1[[29]],plot_l1[[30]],plot_l1[[31]],plot_l1[[32]],
          plot_l1[[33]],plot_l1[[34]],plot_l1[[35]],plot_l1[[36]],
          plot_l1[[37]],plot_l1[[38]],plot_l1[[39]],plot_l1[[40]],
          nrow = 4)
ggsave(paste(out_dir,"cross.rf_imp_sp40.pdf",sep = '/'),width = 20,height = 16)
plot_grid(plot_l1[[41]],plot_l1[[42]],plot_l1[[43]],plot_l1[[44]],
          plot_l1[[45]],plot_l1[[46]],plot_l1[[47]],plot_l1[[48]],
          plot_l1[[49]],plot_l1[[50]],plot_l1[[51]],plot_l1[[52]],
          plot_l1[[53]],plot_l1[[54]],plot_l1[[55]],plot_l1[[56]],
          plot_l1[[57]],
          nrow = 4)
ggsave(paste(out_dir,"cross.rf_imp_sp60.pdf",sep = '/'),width = 20,height = 16)

model_imp_best <- model_imp_sort[model_imp_sort$MeanDecreaseAccuracy>2 & model_imp_sort$Pvalue < 0.05,]
model_imp_best <- model_imp_best[1:28,] 
#tmp_data_train_f <- tmp_data_train[,c(model_imp_best$predictors,"groups")] 

tmp_train_use2 <- sample(nrow(tmp_data), nrow(tmp_data)*0.7)
tmp_data_train_f <- tmp_data[tmp_train_use2,c(model_imp_best$predictors,"groups")]
tmp_data_pre <- tmp_data[! row.names(tmp_data) %in% row.names(tmp_data_train_f),]
tmp_data_pre_f <- tmp_data[! row.names(tmp_data) %in% row.names(tmp_data_train_f),c(model_imp_best$predictors,"groups")]

model_train_1 <- train(x = tmp_data_train_f[,1:(ncol(tmp_data_train_f)-1)],
                     y = factor(tmp_data_train_f$groups,levels = group[c(1,2)]),
                     method = 'rf',# Use the 'random forest' algorithm
                     trControl = trainControl(method='cv', 
                                              number=10, 
                                              search='grid')
                     )

model_train_2 <- randomForest(x = tmp_data_train_f[,1:(ncol(tmp_data_train_f)-1)],
                              y = factor(tmp_data_train_f$groups,levels = group[c(1,2)]),
                              ntree=1000, 
                              importance=TRUE, 
                              proximity=TRUE,
                              mtry = 20
                              )
train_predict <- predict(model_train_2,tmp_data_train_f[,1:(ncol(tmp_data_train_f)-1)],type="response")
test_predict <- predict(model_train_2,tmp_data_pre_f[,1:(ncol(tmp_data_pre_f)-1)],type="response")

confusionMatrix(train_predict, factor(tmp_data_train_f$groups,levels = group[c(1,2)]))
confusionMatrix(test_predict, factor(tmp_data_pre_f$groups,levels = group[c(1,2)]))

tmp_mds_data <- MDSplot(model_train_2, factor(tmp_data_train_f$groups),k=4,pch=30)
tmp_mds_data <- as.data.frame(tmp_mds_data$points) %>% 
  mutate(sample  = row.names(.)) %>% 
  inner_join(.,tmp_meta,by="sample")

ggplot(tmp_mds_data,aes(`Dim 1`,`Dim 2`,color=groups))+
  geom_point(size = 4)+
  theme_bw()+
  theme(
    panel.grid = element_blank()
  )+
  labs(title = "RandomForest MDS plot")+
  scale_color_manual(
    values = color0
  )
ggsave(paste(out_dir,"cross.rf_mds.pdf",sep = '/'),width = 6,height = 4)
train_predict2 <- predict(model_train_2,tmp_data_train_f[,1:(ncol(tmp_data_train_f)-1)],type="vote")
test_predict2 <- predict(model_train_2,tmp_data_pre_f[,1:(ncol(tmp_data_pre_f)-1)],type="vote")

test_roc <- roc(factor(tmp_data_pre_f$groups,levels = group[c(1,2)]), 
                   test_predict2[,1],  
                   plot=T,
                   ci.method="bootstrap",
                   ci =TRUE,
                   conf.level =0.95
)
train_roc <- roc(factor(tmp_data_train_f$groups,levels = group[c(1,2)]), 
                train_predict2[,1],  
                plot=T,
                ci.method="bootstrap",
                ci =TRUE,
                conf.level =0.95
)

roc_with_ci <- function(obj,color_p) {
  ciobj <- ci.se(obj, specificities = seq(0, 1, l = 20))
  dat.ci <- data.frame(x = as.numeric(rownames(ciobj)),
                       lower = ciobj[, 1],
                       upper = ciobj[, 3])
  dat.ci$lower[dat.ci$lower==1] <- 0.99                     
  ggroc(obj, colour =color_p) +theme_minimal() +
    geom_abline(
      slope = 1,
      intercept = 1,
      linetype = "dashed",
      alpha = 0.7,
      color = "grey"
    ) + coord_equal() +
    geom_ribbon(
      data = dat.ci,
      aes(x = x, ymin = lower, ymax = upper),
      fill = color_p,
      alpha = 0.1
    ) + ggtitle(capture.output(obj$auc))
}

roc_with_ci(test_roc,"firebrick3")
ggsave(paste(out_dir,"cross.rf_roc_pre.pdf"))
roc_with_ci(train_roc,"navy")
ggsave(paste(out_dir,"cross.rf_roc_train.pdf"))

#save(model_train_2,model_classify2,file = "workfile/z.rf.Rdata")

phyloseqin_vaild <- subset_samples(phyloseqin,groups %in% group[c(3)])
tmp_data2 <- as.data.frame(t(otu_table(phyloseqin_vaild)))
tmp_data2 <- tmp_data2[,colnames(tmp_data_train_f)[1:ncol(tmp_data_train_f)-1]]
long_predict <- predict(model_train_2,tmp_data2,type="response")

save(model_train_2,model_classify2,
     model_classify,model_classify_pm,
     model_imp_best,model_imp_sort,
     tmp_data_train_f,tmp_data_pre_f,
     tmp_data,tmp_meta,
     file = "workfile/z.rf.Rdata")

data_sp <- read.table("input/stachyose_degradation.txt",header = F) %>% filter(V2 != "non")
which(colnames(data_sp$V1) %in% data_sp$V1)

tmp_train_use2 <- sample(nrow(tmp_data), nrow(tmp_data)*0.7)
tmp_data_train_f <- tmp_data[tmp_train_use2,c(colnames(tmp_data)[which(colnames(tmp_data) %in% data_sp$V1)],"groups")]
tmp_data_pre <- tmp_data[! row.names(tmp_data) %in% row.names(tmp_data_train_f),]
tmp_data_pre_f <- tmp_data[! row.names(tmp_data) %in% row.names(tmp_data_train_f),c(colnames(tmp_data)[which(colnames(tmp_data) %in% data_sp$V1)],"groups")]

model_train_1 <- train(x = tmp_data_train_f[,1:(ncol(tmp_data_train_f)-1)],
                       y = factor(tmp_data_train_f$groups,levels = group[c(1,2)]),
                       method = 'rf',# Use the 'random forest' algorithm
                       trControl = trainControl(method='cv', 
                                                number=10, 
                                                search='grid')
)

model_train_2 <- randomForest(x = tmp_data_train_f[,1:(ncol(tmp_data_train_f)-1)],
                              y = factor(tmp_data_train_f$groups,levels = group[c(1,2)]),
                              ntree=1000, 
                              importance=TRUE, 
                              proximity=TRUE,
)
train_predict <- predict(model_train_2,tmp_data_train_f[,1:(ncol(tmp_data_train_f)-1)],type="response")
test_predict <- predict(model_train_2,tmp_data_pre_f[,1:(ncol(tmp_data_pre_f)-1)],type="response")

confusionMatrix(train_predict, factor(tmp_data_train_f$groups,levels = group[c(1,2)]))
confusionMatrix(test_predict, factor(tmp_data_pre_f$groups,levels = group[c(1,2)]))

tmp_mds_data <- MDSplot(model_train_2, factor(tmp_data_train_f$groups),k=4,pch=30)
tmp_mds_data <- as.data.frame(tmp_mds_data$points) %>% 
  mutate(sample  = row.names(.)) %>% 
  inner_join(.,tmp_meta,by="sample")

ggplot(tmp_mds_data,aes(`Dim 1`,`Dim 2`,color=groups))+
  geom_point(size = 4)+
  theme_bw()+
  theme(
    panel.grid = element_blank()
  )+
  labs(title = "RandomForest MDS plot")+
  scale_color_manual(
    values = color0
  )
ggsave(paste(out_dir,"cross.rf_mds.sta.pdf"),width = 6,height = 4)
train_predict2 <- predict(model_train_2,tmp_data_train_f[,1:(ncol(tmp_data_train_f)-1)],type="vote")
test_predict2 <- predict(model_train_2,tmp_data_pre_f[,1:(ncol(tmp_data_pre_f)-1)],type="vote")

test_roc <- roc(factor(tmp_data_pre_f$groups,levels = group[c(1,2)]), 
                test_predict2[,1],  
                plot=T,
                ci.method="bootstrap",
                ci =TRUE,
                conf.level =0.95
)
train_roc <- roc(factor(tmp_data_train_f$groups,levels = group[c(1,2)]), 
                 train_predict2[,1],  
                 plot=T,
                 ci.method="bootstrap",
                 ci =TRUE,
                 conf.level =0.95
)

roc_with_ci <- function(obj,color_p) {
  ciobj <- ci.se(obj, specificities = seq(0, 1, l = 20))
  dat.ci <- data.frame(x = as.numeric(rownames(ciobj)),
                       lower = ciobj[, 1],
                       upper = ciobj[, 3])
  dat.ci$lower[dat.ci$lower==1] <- 0.99                     
  ggroc(obj, colour =color_p) +theme_minimal() +
    geom_abline(
      slope = 1,
      intercept = 1,
      linetype = "dashed",
      alpha = 0.7,
      color = "grey"
    ) + coord_equal() +
    geom_ribbon(
      data = dat.ci,
      aes(x = x, ymin = lower, ymax = upper),
      fill = color_p,
      alpha = 0.1
    ) + ggtitle(capture.output(obj$auc))
}

roc_with_ci(test_roc,"firebrick3")
ggsave(paste(out_dir,"cross.rf_roc_pre.sta.pdf"))
roc_with_ci(train_roc,"navy")
ggsave(paste(out_dir,"cross.rf_roc_train.sta.pdf"))


phyloseqin_vaild <- subset_samples(phyloseqin,groups %in% group[c(3)])
tmp_data2 <- as.data.frame(t(otu_table(phyloseqin_vaild)))
tmp_data2 <- tmp_data2[,colnames(tmp_data_train_f)[1:ncol(tmp_data_train_f)-1]]
long_predict <- predict(model_train_2,tmp_data2,type="response")

save(model_train_2,
     tmp_data_train_f,tmp_data_pre_f,
     tmp_data,tmp_meta,
     file = "workfile/z.rf.sta.Rdata")

