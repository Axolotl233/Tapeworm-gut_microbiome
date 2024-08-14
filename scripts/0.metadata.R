table_pvalue <- function(x, ...) {
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    p <- wilcox.test(y ~ g)$p.value
  } else {
    p <- fisher.test(table(y, g))$p.value
  }
  c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
}
library(table1)

dir.create("output/metadata")
out_dir <- "output/metadata"
data_m <- read_delim(in_meta_1)
data_m_1 <- data_m %>% filter(groups %in% group[c(1,2)]) %>% filter(sp != "T_saginata")
tmp_list <- combine_list(data_m_1$groups)
data_m_1 <- as.data.frame(data_m_1)

shapiro.test(data_m_1[data_m_1$groups == group[1],"age"])
shapiro.test(data_m_1[data_m_1$groups == group[2],"age"])

table1::table1(data=data_m_1,~sex+age+bmi|groups, 
               extra.col=list(`P-value`=table_pvalue), 
               overall =F)

a <- ggplot(data_m_1,aes(x=groups,y=age)) + 
  geom_boxplot(aes(color=groups),outlier.shape = NA)+
  geom_jitter(aes(color=groups,shape = groups),size = 2.5,width = 0.18,height = 0,alpha = 0.7)+
  theme_cowplot()+
  theme(legend.position = "none")+
  scale_color_manual(values=color0)+
  stat_compare_means(comparisons =  tmp_list,
                     method = "t.test", label = "p.format")
b <- ggplot(data_m_1,aes(x=groups,y=bmi)) + 
  geom_boxplot(aes(color=groups),outlier.shape = NA)+
  geom_jitter(aes(color=groups,shape = groups),size = 2.5,width = 0.18,height = 0,alpha = 0.7)+
  theme_cowplot()+
  theme(legend.position = "none")+
  scale_color_manual(values=color0)+
  stat_compare_means(comparisons =  tmp_list,
                     method = "t.test", label = "p.format")
plot_grid(a,b)
ggsave(paste(out_dir,"age_bmi.cross.pdf",sep = '/'),width=4,height = 4)

data_m <- read_delim(in_meta_2)
data_m  %>% group_by(interval) %>% summarise(count = n())
ggplot(data_m  %>% group_by(interval) %>% summarise(count = n()))+
  geom_bar(aes(x = interval, y = count),stat = "identity",fill = "gray80")+
  geom_text(aes(x = interval, y = count + 0.5,label = count))+
  scale_x_continuous(breaks = seq(0,14,1),limits = c(0,12.5))+
  theme_cowplot()

ggsave(paste(out_dir,"interval.long.pdf",sep = '/'),width=4,height = 4)

ggplot(data_m_1)+
  geom_density(aes(x = age, color = groups))+
  theme_cowplot()
