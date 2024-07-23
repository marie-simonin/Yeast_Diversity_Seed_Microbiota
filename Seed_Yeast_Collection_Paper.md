---
title: "Seed Yeast Collection Paper"
author: "Marie Simonin"
date: "2023-12-13"
output: html_document
---
# 1.Bar graph taxonomy of isolate distribution among plants - Figure 3
```{r}
isolates <- read.table("yeast-strain-list.txt", header=TRUE, check.names = FALSE, sep = "\t")
dim(isolates)

Plant_colors <-  c('Phaseolus vulgaris'='#8cb369','Brassica napus'='#ffb563' ,'Solanum lycopersicum'='#F6635C', "Triticum aestivum"="#8ecae6","Erophila verna"="#06d6a0", "Capsella bursa-pastoris"="#ffc8dd", "Raphanus sativus"="#023047", "Cardamine hirsuta"="#4361ee")
```
#Figure 3
```{r}
library(ggh4x)
Figure3 <- ggplot(data=isolates, aes(y=Genus, fill=Plant))+ geom_bar(aes(), position="stack")+ theme_classic()+ theme(axis.title = element_text(color="black", size=12, face="bold"))+ theme(axis.text.y = element_text(color="black", size=10, face="bold"))+theme(axis.text.x = element_text(color="black", size=10, face="bold"))	+ theme(legend.text = element_text(color="black", size=10, face="bold"))+ theme(legend.title = element_text(color="black", size=12, face="bold"))+xlab("Number of yeast isolates ")+scale_fill_manual(values = c(Plant_colors))+ facet_nested(Habitat+Class~., scales = "free_y", space = "free")+ theme(strip.text.y = element_text(size=10, face = "bold", angle = 0))
Figure3
#ggsave(Figure3, dpi=300, device = png, width = 13, height = 9, filename = "Figure3.png")
```


# bar graph seed-seedling (not included in article)
```{r}
bargraph1 <- ggplot(data=isolates, aes(x=Plant, fill=Habitat))+ geom_bar(aes(), position="stack")+ theme_classic()+ theme(axis.title = element_text(color="black", size=12, face="bold"))+ theme(axis.text.x = element_text(color="black", size=9, face="bold", angle = 45, hjust = 1))+ theme(axis.text.y = element_text(color="black", size=10, face="bold"))	+ theme(legend.text = element_text(color="black", size=10, face="bold"))+ theme(legend.title = element_text(color="black", size=12, face="bold"))+ylab("Number of yeast isolates")
bargraph1
```
#bar graph composition percentage (not included in article)
```{r}
bargraph <- ggplot(data=isolates, aes(x=Plant, fill=Genus))+facet_grid(.~Habitat, scales="free", space = "free")+ geom_bar(aes(), position="fill")+ theme_classic()+ theme(axis.title = element_text(color="black", size=12, face="bold"))+ theme(axis.text.x = element_text(color="black", size=9, face="bold", angle = 45, hjust = 1))+ theme(axis.text.y = element_text(color="black", size=10, face="bold"))	+ theme(legend.text = element_text(color="black", size=10, face="bold"))+ theme(legend.title = element_text(color="black", size=12, face="bold"))+ylab("Proportion of yeast isolates (%)")+ theme(strip.text.x = element_text(size=12, face = "bold"))+scale_y_continuous(labels =scales::percent)
bargraph
```

#2. Yeast isolate in Seed Microbiota Database - Figure 4


### Import ITS dataset from Seed Microbiota Database - Subset 3
```{r}
meta2 <- read.table("Metadata_ITS1_ITS2_withDivSubset2_Jan2021.txt", header=TRUE, check.names = FALSE, sep = "\t")
SV<-read.table("Subset3-ITS1_table-FINAL-rarefied-transposed.txt", header=TRUE, check.names = FALSE, sep = "\t")
head(SV)
head(meta2)
dim(meta2)
dim(SV)
SV_ITS1_R<-merge(meta2,SV,by="SampleID")
dim(SV_ITS1_R)
head(SV_ITS1_R)
``` 

```{r}
matrix<-SV_ITS1_R[c(34:ncol(SV_ITS1_R))]
dim(matrix)
head(matrix)
```

```{r}
##Keeping ASVs with at least 100 reads in the meta-analysis dataset
dim(matrix)
matrix_use<-matrix[,colSums(matrix)>=100]
dim(matrix_use)

#transposing
matrix_uset=t(matrix_use)
head(matrix_uset)
dim(matrix_uset)
```


```{r}
## Calculate prevalence and relative abundance for each ASV
#Code Shade lab: https://github.com/ShadeLab/PAPER_Shade_CurrOpinMicro/blob/master/script/Core_prioritizing_script.R
#presence-absence data
SV_ITS1_PA <- 1*((matrix_uset>0)==1)                                              
# occupancy calculation
SV_ITS1_Prevalence <- rowSums(SV_ITS1_PA)/ncol(SV_ITS1_PA) 
# relative abundance  
library(vegan)
library(dplyr)
SV_ITS1_relative_abundance <- apply(decostand(matrix_uset, method="total", MARGIN=2),1, mean)     

# combining occupancy and relative abundance of each SV_ITS1 in a table
SV_ITS1prev_rel <- add_rownames(as.data.frame(cbind(SV_ITS1_Prevalence, SV_ITS1_relative_abundance)),'SV') 
head(SV_ITS1prev_rel)
dim(SV_ITS1prev_rel)

```


```{r}
## Merge prevalence/rel abund data with taxonomic info for each SV
taxo_ITS1<-read.table("Subset1-2_All_studies_merged_ITS1_taxonomy.tsv", header=TRUE, check.names = FALSE, sep = "\t")
SV_ITS1prev_rel_taxo<-merge(SV_ITS1prev_rel,taxo_ITS1,by="SV")
dim(SV_ITS1prev_rel_taxo)

head(SV_ITS1prev_rel_taxo)
```



##Figure 4 Abundance-occupancy graph - ITS1 - Meta-analysis
```{r}
library(ggplot2)
SV_ITS1prev_rel_taxo$Isolate<-ordered(SV_ITS1prev_rel_taxo$Isolate, levels=c("Yes", "Others"))	
color= c("#ee8332", "#4a1777")
figure4=ggplot(data=SV_ITS1prev_rel_taxo, aes(x=SV_ITS1_relative_abundance, y=SV_ITS1_Prevalence)) +
    geom_point(aes(shape=Isolate,color=Core, size=Number_isolates)) + xlab("Mean Taxa (ASV) Relative Abundance") + ylab("Taxa (ASV) Prevalence")+scale_x_log10(labels = scales::percent_format(accuracy = 0.1))+ scale_y_log10(labels = scales::percent_format(accuracy = 1))+ theme_classic()+ theme(axis.title = element_text(color="black", size=12, face="bold"))+ theme(axis.text = element_text(color="black", size=11, face="bold"))+ theme(legend.text = element_text(color="black", size=12, face="bold"))+ theme(legend.title = element_text(color="black", size=12, face="bold")) +ggtitle("Seed Microbiota Database - ITS1 region")+theme(plot.title = element_text(hjust = 0.5, face="bold"))+scale_shape_manual(values=c(16, 1))+scale_color_manual(values=color)+ labs(shape = "Isolate", color = "Taxon type", size = "Number of Isolates")
figure4
ggsave(figure4, dpi=300, device = png, width = 8, height = 6, filename = "Figure4.png")
```


```{r}
###Subset only strains to prepare the table in Figure 1
SV_ITS1prev_rel_taxo_strain=subset(SV_ITS1prev_rel_taxo, Isolate!="Others")
SV_ITS1prev_rel_taxo_strain
dim(SV_ITS1prev_rel_taxo_strain) 
```




