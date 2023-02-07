#The output from HUMAnN and QIIME (distance matrices) are then processed in R for plots and statistical analysis.
```{r}
library(tidyverse)
library(vegan)
library(pairwiseAdonis)
library(ggplot2)
library(cowplot)
library(glmmTMB)
library(multcomp)
library(car)
library(fdrtool)
library(pheatmap)

setwd("C:/Users/user/Desktop/Pathways/NMDS/")
```

```{r}

bray_gf = as.dist(read.table("bray-curtis-distance-matrix-gf.tsv", header = T))

jaccard_gf = as.dist(read.table("jaccard-distance-matrix-gf.tsv", header = T))

bray_pa = as.dist(read.table(file.choose(), header = T)) 

jaccard_pa = as.dist(read.table("jaccard-distance-matrix-pa.tsv", header = T))

bray_cog_full = as.dist(read.table("bray-curtis-distance-matrix-cog.tsv", header = T))

jaccard_cog_full = as.dist(read.table("jaccard-distance-matrix-cog.tsv", header = T))

metadata_gf = read.delim("pepper_metadata_gf.txt", header = T)

metadata_pa = read.delim(file.choose(), header = T)  

metadata_cog_full = read.csv("pepper_metadata_gf.txt", header = T)
```



```{r}

library(vegan)

adonis2(jaccard_pa ~ treatment + innoculation + cultivar + time, data = metadata_pa, 
        by = "margin", permutations = 4999)

adonis2(bray_pa ~ treatment + innoculation + cultivar + time, data = metadata_pa, 
        by = "margin", permutations = 4999)


pairwise.adonis(bray_pa, factors = metadata_pa$treatment, perm = 4999, p.adjust.m = "BH")

pairwise.adonis(bray_pa, factors = metadata_pa$innoculation, perm = 4999, p.adjust.m = "holm")

adonis2(bray_gf ~ cultivar + innoculation, data=metadata_gf, 
        by = "margin", permutations = 4999)


mds_bray_pa<-metaMDS(bray_pa, k=2, trymax=499)
mds_bray_pa_points<-mds_bray_pa$points
mds_bray_pa_points2<-merge(x=mds_bray_pa_points, 
                           y = metadata_pa, 
                                  by.x = "row.names", by.y = "name")
braypa <- ggplot(mds_bray_pa_points2, aes(x = MDS1, y = MDS2, 
                                                 color = treatment, shape = innoculation)) + 
    geom_point(size = 4) +
   scale_color_brewer(palette = 'Set1') + ggtitle("NMDS of pathways across treatment")+  annotate("text", x = -0.2, y = 0.03, label = "Stress = 0.003167382")+ theme(legend.position="bottom")

braypa

ggsave("NMDSpathways.pdf",width = 25, height = 15, units="cm", dpi=700)


#For Cultivar
dune.dist <- vegdist(mds_bray_pa, method="bray")
dispersionC <- betadisper(dune.dist, group=NMDS.env$ultivar)
permutest(dispersionC)
plot(dispersionC, hull=FALSE, ellipse=TRUE) ##sd ellipse

ano = anosim(NMDS, NMDS.env$Cultivar, distance = "bray", permutations = 999)
ano
