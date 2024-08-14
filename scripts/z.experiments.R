dir.create("output/experiment",showWarnings = F)
out_dir <- "output/experiment"

tmp_ed1 <- read_csv(in_experment1)
tmp_ed1$class <- factor(tmp_ed1$class,levels = rev(unique(tmp_ed1$class)))
tmp_list=combine_list(tmp_ed1$class)

ggplot(tmp_ed1,aes(x=class,y=CFU,color=class,fill = class)) +
  stat_summary(fun=mean, geom="bar",position = position_dodge(0.75),width=0.75,alpha = 0.4)+
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               position = position_dodge(0.75),width = 0.2)+
  geom_jitter(position = position_jitterdodge(dodge.width = 0.75,
                                                           jitter.width = 0,
                                                           jitter.height = 0),size = 3)+
  theme_cowplot()+
  labs(x="")+
  facet_wrap(species~.)+
  stat_compare_means(method="t.test",
                     comparisons = tmp_list,
                     label = "p.format")+
  scale_color_manual(values = c("#6EBCBC","#37526D"))+
  scale_fill_manual(values = c("#6EBCBC","#37526D"))
ggsave(paste(out_dir,"esp_wyg.pdf",sep = '/'),width = 5,height = 5) 

tmp_ed2 <- read_csv(in_experment2)
tmp_ed2$concentrate <- factor(tmp_ed2$concentrate,levels = rev(unique(tmp_ed2$concentrate)))
tmp_list <- combine_list(tmp_ed2$concentrate,y = "control")

ggplot(tmp_ed2,aes(x=concentrate,y=value,color=concentrate,fill = concentrate))+
  stat_summary(fun=mean, geom="bar",position = position_dodge(0.75),width=0.75,alpha = 0.4)+
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               position = position_dodge(0.75),width = 0.3)+
  geom_jitter(position = position_jitterdodge(dodge.width = 0.75,
                                              jitter.width = 0,
                                              jitter.height = 0),size = 2)+
  theme_cowplot()+
  theme(
    axis.text.x = element_text(angle = 30,hjust = 1),
    legend.position = "none"
  )+
  labs(x="")+
  facet_wrap(species~.,scales = "free_y")+
  stat_compare_means(method="t.test",
                     comparisons = rev(tmp_list),
                     label = "p.format")+
  scale_color_manual(values = colorx)+
  scale_fill_manual(values = colorx)
ggsave(paste(out_dir,"esp_dyy.pdf",sep = '/'),width = 7,height = 6) 

tmp_ed3 <- read_csv(in_experment3,col_names = F)
colnames(tmp_ed3)  <-  c("concentration","OD","species","time")
tmp_ed3$concentration <- paste(tmp_ed3$concentration,"%",sep = "")
tmp_ed3 %<>% filter(time < 21)
topbar <- function(x){      
  return(mean(x)+sd(x)/sqrt(length(x))) #误差采用了mean+-sem
}
bottombar <- function(x){
  return(mean(x)-sd(x)/sqrt(length(x)))
}
tmp_ed3_1 <- tmp_ed3 %>% filter(species == "Bfido",concentration == "1%")
tmp_ed3_2 <- tmp_ed3 %>% filter(species == "Doera")

ggplot()+
  stat_summary(data = tmp_ed3,aes(x=time,y=OD,color=concentration),
               geom = 'line',fun='mean',linewidth = 1)+
  stat_summary(data = tmp_ed3,aes(x=time,y=OD,color=concentration),
               geom = 'errorbar',fun.min = bottombar,fun.max = topbar,
               linewidth = 1,width = 0.5)+
  stat_summary(data = tmp_ed3,aes(x=time,y=OD,fill=concentration),
               geom = 'point',fun='mean',
               size=2,pch=21,color='black')+
  theme_cowplot()+
  scale_x_continuous(breaks = c(0,8,12,16,20,24))+
  scale_y_continuous(breaks = seq(0,1,0.1))+
  #labs(title = "12_Ketolithocholic_Acid" , y = "inhibition rate", x = "Concentration")+
  scale_color_manual(values = colorx)+
  scale_fill_manual(values = colorx)+
  facet_wrap(.~species,scales = "free_y",nrow = 1)
ggsave(paste(out_dir,"sta_lyq.pdf",sep = '/'),width = 7.5,height = 5)



