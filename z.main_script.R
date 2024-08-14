rm(list=ls())
dir.create("output")
dir.create("workfile")

#library(ggstatsplot)
library(tidyverse)
library(magrittr)
library(phyloseq)
library(ggpubr)
library(cowplot)
library(ggrepel)
library(ggpmisc)
library(RColorBrewer)
library(aplot)
library(speedyseq)
#library(calecopal)
#library(MoMAColors)
source("E:/work/z.code/common_custom_function.R")

#=====> global variant

in_meta_1 <- "input/metadata_2309_1.txt"
in_meta_2 <- "input/metadata_2309_3.txt"
in_otu_ab <- "input/metaphlan.ab.txt"
in_otu_count <- "input/metaphlan.count.txt"
in_pathway_metacyc <- "input/humann.pathabundance.metacyc.cpm.tsv"
in_pathway_unipath <- "input/humann.pathabundance.unipath.cpm.tsv"
in_metabolite <- "input/metabolite.csv"
in_experment1 <- "input/Experiment-wyg1.csv"
in_experment2 <- "input/Experiment-dyy_esp1.csv"
in_experment3 <- "input/Experiment-lyq.csv"

group <- c("HC","TA","Anti_TA")
groupl <- c("HC","BE","AF","TA","Anti_TA")
sp <- c("T_asiatica","T_saginata")

color0 <- c("#1B80AD","#DA5956","#Cf9860")
color1 <- c("#4d8f8b","#ad5a6b")

#color1 <- c("#5E6B7B","#286A81")
#color0 <- c("#4d8f8b","#ad5a6b","#Cf9860")
#color1 <- c("#5E6A8B","#b80422")
#color0 <- c("#95caa6","#f19670","#e9d78e")
#color1 <- c("#008D98","#c94a53")
#color0 <- c("#B9C7E2","#ECAB99","#DCC27A")
#color1 <- c("#84A6A2","#AA767C")
#color1 <- c("#286A81","#954150")
colorx <- c("#A2A098","#5E6B7B","#447c69","#286A81","#BE5A47","#375377","#AA767C")

#color3 <- colorRampPalette(c("navy","cyan","springgreen","yellow","coral","firebrick3"))
color_bar<-sample(c("#51574a","#BE5A47","#74c493","#8e8c6d","#e4bf80","#e9d78e","#e2975d","#f19670","#e16552",
                    "#c94a53","#be5168","#a34974","#993767","#65387d","#4e2472","#9163b6","#e279a3","#e0598b","#7c9fb0",
                    "#5698c4","#9abf88"),21)
color_bar_st <- c("#5698c4","#4e2472","#7c9fb0","#a34974","#e16552","#be5168","#9abf88","#8e8c6d","#c94a53","#993767",
                  "#447c69","#f19670","#74c493","#e279a3","#e4bf80","#9163b6","#e2975d","#e0598b","#65387d","#e9d78e",
                  "#51574a")
#outlier_sample = c("T65be","T65AF")
outlier_sp <- c("s__GGB1237_SGB1623","s__Gemmiger_SGB15292")
