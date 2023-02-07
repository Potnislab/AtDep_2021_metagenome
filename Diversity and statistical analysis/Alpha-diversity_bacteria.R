---
title: "Untitled"
output: html_document
date: '2022-04-06'
---
```{r}
#Load the library
library(vegan)
library(tidyverse)
library(devtools)
library(grid) 
library(gridExtra)
library(knitr)
library(png)
library(ggpubr) 
library(RColorBrewer)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(phyloseq)
library(ggnewscale)
library(SRS)
library(MicrobiotaProcess)
library(FSA)
```

```{r}
#load the biom file and metadata
data <- import_biom(file.choose(),parseFunction=parse_taxonomy_default)

metadata <- read.csv(file.choose(), header=TRUE, 
                     row.names=1, check.names=F, comment.char="" )

names(data)
sample_names(data)
rownames(metadata)
all(rownames(metadata) %in% sample_names(data))

sample_data <- sample_data(metadata)

physeq2 <- merge_phyloseq(sample_data, data )
physeq2
View(physeq2@tax_table@.Data)
physeq2@tax_table@.Data <- substring(data@tax_table@.Data, 4)
colnames(physeq2@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
colnames(physeq2@tax_table@.Data)

justbacteria <- subset_taxa(physeq2, !Kingdom %in% c("Eukaryota"))
justbacteria2 <- subset_taxa(justbacteria, !Kingdom %in% c("Viruses"))
justbacteria3 <- subset_taxa(justbacteria2, !Phylum %in% c("Chordata"))

#similarly lets remove the mid season with high xanthomonas
allsamples<- subset_samples(justbacteria3, id != "7EM" & id !="8EM")
allsamples

#To know number of reads per sample and sort them 
sample_sums(allsamples)
sort(sample_sums(allsamples))

#SRS
minControl=425449
keep=names(which(sample_sums(allsamples)>=minControl))
pruned_all=prune_samples(keep, allsamples)
any(taxa_sums(pruned_all) == 0)

pruned_all <- prune_taxa(taxa_sums(pruned_all) > 0, pruned_all)
any(taxa_sums(pruned_all) == 0)



pruned_data_df = as.data.frame(otu_table(pruned_all))

SRS_OUTPUT <- SRS (pruned_data_df,minControl)

rownames(SRS_OUTPUT)=rownames(pruned_data_df)

table(sample_data(pruned_all)[,"HET"])

# transform back into phyloseq object
taxa = otu_table(SRS_OUTPUT, taxa_are_rows = TRUE)
otu_table(pruned_all)=taxa
any(taxa_sums(pruned_all) == 0)

pruned_all <- prune_taxa(taxa_sums(pruned_all) > 0, pruned_all)
any(taxa_sums(pruned_all) == 0)

SeqDepthPruned = sample_sums(pruned_all)
sample_data(pruned_all)$SeqDepthPruned = SeqDepthPruned

# barplot of library sizes
library(ggplot2)
ggplot(meta(pruned_all), aes(id, SeqDepthPruned)) + geom_bar(stat = "identity", aes(fill = environment)) +
  rotate_x_text()

pruned_all

```

```{r}
#Alpha diversity plots for Chao1 and Shannon
A1= plot_richness(pruned_all, x="IET", measures= "Chao1") + 
  theme_grey() + 
  geom_boxplot(aes(fill=sample_data(pruned_all)$time),width=0.5,lwd=0.2) +  
  ggtitle("Chao1 richness") + 
  scale_fill_manual(values=c( "#40E0D0", "#F4A460", "#56B4E9"))  + 
  facet_wrap(~cultivar)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "none") +stat_compare_means(aes(group=cultivar), label = "p.signif", method="anova")

#To put the samples in a specific order i first created a new order
newSTorder = c("CAB","CAM","CAE","COM","COE","IAB","IAM","IAE","IOM","IOE")
#Then convert the x axis factor into character vector
A1$data$IET <- as.character(A1$data$IET)

#Then change the character vector 
A1$data$IET <- factor(A1$data$IET, levels=newSTorder)

A1$layers
A1$layers <- A1$layers[-1]

 N1 <-A1
#Now for shannon diversity

A2 = plot_richness(pruned_all, x="IET", measures= "Shannon") + 
  theme_grey() + 
  geom_boxplot(aes(fill=sample_data(pruned_all)$time), width=0.5,lwd=0.2) +  
  ggtitle("Shannon richness") + 
  scale_fill_manual(values=c( "#40E0D0", "#F4A460", "#56B4E9"))  + 
  stat_compare_means(aes(label = "wilcox.test"), label.x = NULL, label.y = 4.3)  + facet_wrap(~cultivar)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "none")

#To put the samples in a specific order i first created a new order
newSTorder = c("CAB","CAM","CAE","COM","COE","IAB","IAM","IAE","IOM","IOE")
#Then convert the x axis factor into character vector
A2$data$IET <- as.character(A2$data$IET)

#Then change the character vector 
A2$data$IET <- factor(A2$data$IET, levels=newSTorder)

#Use gg arrange to combine two figures into one
ggarrange(A1, A2, 
          labels = c("A", "B"),
          ncol = 1, nrow = 2)
ggsave("Alpha.pdf",width = 30, height = 23, units="cm", dpi=700)


ggarrange(N1,N2,F2,F1, 
          labels = c("A","B", "C", "D"),
          ncol = 2, nrow = 2)
ggsave("mix2.pdf",width = 24, height = 14, units="cm", dpi=700)

#To remove the points 
A2$layers

A2$layers <- A2$layers[-1]

N2 <- A2


```



```{r}
#Lets subset X10R samples

X10R<- subset_samples(pruned_all, cultivar != "ECW")
meta(X10R)
richnessX10R<- estimate_richness(X10R)
plot_richness (X10R)
shapiro.test(richnessX10R$Chao1)
shapiro.test(richnessX10R$Shannon)

#For Chao1 
anova.sh = aov(richnessX10R$Chao1 ~ sample_data(X10R)$IET)
summary(anova.sh)

TukeyHSD(anova.sh)

#For non normal shannon index
kruskal.test(richnessX10R$Shannon ~ sample_data(X10R)$IET)

#If kruskal is significant perform the following test for pairwise comparison
pairwise.wilcox.test(richnessX10R$Shannon, sample_data(X10R)$IET,   p.adjust.method = "BH")
#Dunntest for pairwise comparison

dunnTest(richnessX10R$Shannon, sample_data(X10R)$IET, method="bh")

#Now for ECW cultivars and Chao1
ECW<- subset_samples(pruned_all, cultivar != "X10R")

richnessECW<- estimate_richness(ECW)
plot_richness (ECW)
shapiro.test(richnessECW$Shannon)
shapiro.test(richnessECW$Chao1)

#For Chao1 
set.seed(123)
anova.sh = aov(richnessECW$Chao1 ~ sample_data(ECW)$IET)
summary(anova.sh)
TukeyHSD(anova.sh)

#For non normal shannon index we need to break data further down into ozone and ambient env and remove the base samples to reduce more factors
pruned_ECWoz<- subset_samples(ECW, environment !="Ambient" & IET != "CAB" & IET !="IAB")
meta(pruned_ECWoz)
richnessECWoz<- estimate_richness(pruned_ECWoz)
#Shanon
shapiro.test(richnessECWoz$Shannon)
kruskal.test(richnessECWoz$Shannon ~ sample_data(pruned_ECWoz)$IET)
dunnTest(richnessECWoz$Shannon, sample_data(pruned_ECWoz)$IET, method="bh")
#Chao1
shapiro.test(richnessECWoz$Chao1)
anova.sh = aov(richnessECWoz$Chao1 ~ sample_data(pruned_ECWoz)$IET)
summary(anova.sh)
TukeyHSD(anova.sh)


#For Ambient 
pruned_ECWamb<- subset_samples(ECW, environment !="Ozone" & IET != "CAB" & IET !="IAB")
meta(pruned_ECWamb)
richnessECWamb<- estimate_richness(pruned_ECWamb)
shapiro.test(richnessECWamb$Shannon)
#As shannon is normally distributed we will use anoval here
anova.sh = aov(richnessECWamb$Shannon ~ sample_data(pruned_ECWamb)$IET)
summary(anova.sh)
TukeyHSD(anova.sh)


#Now to look at the effect of innoculation and envon both the cultivars. Here i looked at the effect of environment on both innoculated and control samples and found no significance difference in any of the samples so i have code for one here
pruned_ECWinoc<- subset_samples(ECW, innoculation !="Control")
meta(pruned_ECWinoc)
richnessECWinoc<- estimate_richness(pruned_ECWinoc)
#Shapiro
shapiro.test(richnessECWinoc$Shannon)
kruskal.test(richnessECWinoc$Shannon ~ sample_data(pruned_ECWinoc)$IET)
dunnTest(richnessECWinoc$Shannon, sample_data(pruned_ECWinoc)$IET, method="bh")

#Chao1
shapiro.test(richnessECWinoc$Chao1)
anova.sh = aov(richnessECWinoc$Chao1 ~ sample_data(pruned_ECWinoc)$IET)
summary(anova.sh)
TukeyHSD(anova.sh)

```



```{r}
#Beta diversity for all the samples
#To know number of reads per sample and sort them 
Beta_all<- subset_samples(allsamples,  IET != "CAB" & IET !="IAB")
sample_sums(Beta_all)
sort(sample_sums(Beta_all))

#SRS
minControl=2548613
keep=names(which(sample_sums(Beta_all)>=minControl))
prunedBeta_all=prune_samples(keep, Beta_all)
any(taxa_sums(prunedBeta_all) == 0)

prunedBeta_all <- prune_taxa(taxa_sums(prunedBeta_all) > 0, prunedBeta_all)
any(taxa_sums(prunedBeta_all) == 0)


library(SRS)

pruned_data_df = as.data.frame(otu_table(prunedBeta_all))

SRS_OUTPUT <- SRS (pruned_data_df,minControl)

rownames(SRS_OUTPUT)=rownames(pruned_data_df)

table(sample_data(prunedBeta_all)[,"HET"])

# transform back into phyloseq object
taxa = otu_table(SRS_OUTPUT, taxa_are_rows = TRUE)
otu_table(prunedBeta_all)=taxa
any(taxa_sums(prunedBeta_all) == 0)

prunedBeta_all <- prune_taxa(taxa_sums(prunedBeta_all) > 0, prunedBeta_all)
any(taxa_sums(prunedBeta_all) == 0)

SeqDepthPruned = sample_sums(prunedBeta_all)
sample_data(prunedBeta_all)$SeqDepthPruned = SeqDepthPruned

# barplot of library sizes
library(ggplot2)
ggplot(meta(prunedBeta_all), aes(id, SeqDepthPruned)) + geom_bar(stat = "identity", aes(fill = environment)) +
  rotate_x_text()

prunedBeta_all

#beta-diversity Remove the base samples first
nmds.all <- ordinate(prunedBeta_all, method = "NMDS", distance = "bray")
#Find stess

nmds.all$stress
stressplot(nmds.all)

p <- plot_ordination(prunedBeta_all, nmds.all,  type="samples", 
                      color="cultivar",
                      shape="IET")+   scale_shape_manual(values=c(3, 15, 16, 17, 18, 8,9,10))+
  scale_colour_manual(values=c("ECW" = "#FF7F00", "X10R" = "#5D3FD3"))+ geom_point(size=4)+ ggtitle("NMDS between two cultivars")+  annotate("text", x = -0.5, y = 0.5, label = "Stress = 0.1154452")


p1=p + 
  stat_ellipse(type = "norm", linetype = 2) +
  theme_gray()+theme(
    axis.text.x = element_text(size = 12, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 12), legend.position ="bottom")


ggarrange(p,NULL, 
          labels = c("A",""),
          ncol = 2, nrow = 1, widths = c(1, 0.2))
ggsave("hello10.pdf",width = 15, height = 15, units="cm", dpi=700)
ggsave("NMDS_Cultivar&time.pdf",width = 20, height = 15, units="cm", dpi=700)

#beta-diversity on just control plants
Beta_control <- subset_samples(prunedBeta_all, innoculation != "Innoculated")
nmds.control <- ordinate(Beta_control, method = "NMDS", distance = "bray")

nmds.control$stress

p2 <- plot_ordination(Beta_control, nmds.control,  type="samples", 
                      color="cultivar", 
                      shape="time")+  
  scale_fill_brewer(palette="Paired")+ geom_point(size=4, alpha=0.75)+ ggtitle("NMDS between control samples")+  annotate("text", x = -0.5, y = 0.5, label = "Stress = 0.1162348")


p2 + 
  stat_ellipse(type = "norm", linetype = 2) +
  theme_gray()+theme(
    axis.text.x = element_text(size = 10, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 12), legend.position ="bottom")

#Betadiversity on control samples
Beta_innoc <- subset_samples(prunedBeta_all, innoculation != "Control")
meta(Beta_innoc)
nmds.inoc <- ordinate(Beta_innoc, method = "NMDS", distance = "bray")
nmds.inoc$stress

p3 <- plot_ordination(Beta_innoc, nmds.inoc,  type="samples", 
                      color="cultivar", 
                      shape="time")+  
  scale_fill_brewer(palette="Set1")+ geom_point(size=4, alpha=0.75)+ ggtitle("NMDS between innoculated samples")+  annotate("text", x = -0.5, y = 0.5, label = "Stress = 0.06120877")

p3 + 
  stat_ellipse(type = "norm", linetype = 2) +
  theme_gray()+theme(
    axis.text.x = element_text(size = 10, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 12), legend.position ="bottom")
    

#Use gg arrange to combine two figures into one
ggarrange(p2, p3, 
          labels = c("A", "B"),
          ncol = 1, nrow = 2)
ggsave("NMDS_treatment.pdf",width = 20, height = 30, units="cm", dpi=700)

#NMDS for environmental variable for all samples
E1 <- plot_ordination(prunedBeta_all, nmds.all,  type="samples", 
                      color="environment",
                      shape="cultivar")+  
  scale_fill_brewer(palette="Set1")+ geom_point(size=4, alpha=0.75)+ ggtitle("NMDS between two cultivars")+  annotate("text", x = -0.5, y = 0.5, label = "Stress = 0.1154452")

E1+ 
  stat_ellipse(type = "norm", linetype = 2) +
  theme_gray()+theme(
    axis.text.x = element_text(size = 10, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 12), legend.position ="bottom")

ggsave("NMDSall_Environment.pdf",width = 20, height = 20, units="cm", dpi=700)

#Lets subset them based on cultivars
ECWnmds<- subset_samples(prunedBeta_all, cultivar !="X10R")
nmds.ECW <- ordinate(ECWnmds, method = "NMDS", distance = "bray")
nmds.ECW$stress

E2 <- plot_ordination(ECWnmds, nmds.ECW,  type="samples", 
                      color="environment", 
                      shape="environment")+  
  scale_fill_brewer(palette="Set1")+ geom_point(size=4, alpha=0.75)+ ggtitle("NMDS between env in ECW")+  annotate("text", x = -0.5, y = 0.5, label = "Stress = 0.06357417")

E2 + 
  stat_ellipse(type = "norm", linetype = 2) +
  theme_gray()+theme(
    axis.text.x = element_text(size = 10, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 12), legend.position ="bottom")
    
#lets subset again to control ECW samples
ECWcontnmds<- subset_samples(ECWnmds, innoculation !="Innoculated")
nmdscont.ECW <- ordinate(ECWcontnmds, method = "NMDS", distance = "bray")
nmdscont.ECW$stress

E3 <- plot_ordination(ECWcontnmds, nmdscont.ECW,  type="samples", 
                      color="environment", 
                      shape="environment")+  
  scale_fill_brewer(palette="Set1")+ geom_point(size=4, alpha=0.75)+ ggtitle("NMDS between env in ECW control samples")+  annotate("text", x = -0.5, y = 0.5, label = "Stress = 0.04583639")

E3 + 
  stat_ellipse(type = "norm", linetype = 2) +
  theme_gray()+theme(
    axis.text.x = element_text(size = 10, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 12), legend.position ="bottom")

#ECW innoculated samples
#lets subset again to control ECW samples
ECWinnocnmds<- subset_samples(ECWnmds, innoculation !="Control")
nmdsinnoc.ECW <- ordinate(ECWinnocnmds, method = "NMDS", distance = "bray")
nmdsinnoc.ECW$stress

E4 <- plot_ordination(ECWinnocnmds, nmdsinnoc.ECW,  type="samples", 
                      color="environment", 
                      shape="environment")+  
  scale_fill_brewer(palette="Set1")+ geom_point(size=4, alpha=0.75)+ ggtitle("NMDS between env in ECW innoculated samples")
#annotate("text", x = 4e-04, y = 0.5, label = "Stress = 0.04583639")

E4 + 
  stat_ellipse(type = "norm", linetype = 2) +
  theme_gray()+theme(
    axis.text.x = element_text(size = 10, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 12), legend.position ="bottom")

#Now for X10R samples
envPalette= c()
X10Rnmds<- subset_samples(prunedBeta_all, cultivar !="ECW")
nmds.X10R <- ordinate(X10Rnmds, method = "NMDS", distance = "bray")
nmds.X10R$stress

X1 <- plot_ordination(X10Rnmds, nmds.X10R,  type="samples", 
                      color="environment", 
                      shape="environment")+  
  scale_fill_brewer(palette="Set1")+ geom_point(size=4, alpha=0.75)+ ggtitle("NMDS between env in X10R")+  annotate("text", x = -0.5, y = 0.5, label = "Stress = 0.109644")

X1 + 
  stat_ellipse(type = "norm", linetype = 2) +
  theme_gray()+theme(
    axis.text.x = element_text(size = 10, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 12), legend.position ="bottom")

#lets subset again to control X10R samples
X10Rcontnmds<- subset_samples(X10Rnmds, innoculation !="Innoculated")
nmdscont.X10R <- ordinate(X10Rcontnmds, method = "NMDS", distance = "bray")
nmdscont.X10R$stress

X2 <- plot_ordination(X10Rcontnmds, nmdscont.X10R,  type="samples", 
                      color="cultivar",
                      shape="IET")+   scale_shape_manual(values=c(3, 15, 16, 17, 18, 8,9,10))+
  scale_colour_manual(values=c("ECW" = "#FF7F00", "X10R" = "#5D3FD3"))+ geom_point(size=4)+ ggtitle("NMDS between two cultivars")+  annotate("text", x = -0.25, y = 0.5, label = "Stress = 0.06014424")
  
ggsave("X10R.pdf",width = 10, height = 12, units="cm", dpi=700)                    

X2 + 
  stat_ellipse(type = "norm", linetype = 2) +
  theme_bw()+theme(
    axis.text.x = element_text(size = 10, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 12), legend.position ="bottom")

#ECW innoculated samples
#lets subset again to control ECW samples
X10Rinnocnmds<- subset_samples(X10Rnmds, innoculation !="Control")
nmdsinnoc.X10R <- ordinate(X10Rinnocnmds, method = "NMDS", distance = "bray")
nmdsinnoc.X10R$stress

X3 <- plot_ordination(X10Rinnocnmds, nmdsinnoc.X10R,  type="samples", 
                      color="environment", 
                      shape="environment")+  
 scale_colour_manual(values=c("Ambient" = "Darkgreen", "Ozone" = "Red"))+ geom_point(size=4, alpha=0.75)+ ggtitle("NMDS between env in X10R innoculated samples")+ annotate("text", x = 4e-04, y = 0.5, label = "Stress = 0.05688418")+ scale_shape_manual(values=c(15, 18))

X3=X3 + 
  stat_ellipse(type = "norm", linetype = 2) +
  theme_gray()+theme(
    axis.text.x = element_text(size = 12, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 12), legend.position ="bottom")

#Save the plot from cultivar nmds and x10R innoculated
#Use gg arrange to combine two figures into one
ggarrange(p1, X3, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)
ggsave("NMDS_treatment.pdf",width = 30, height = 20, units="cm", dpi=700)
```



```{r}
#Sometimes it will be clear from nMDS that one group tends to vary more (be more spread out) than another group. We can test this statistically with multivariate homogeneity of group dispersion (variances).
#Lets also divide the samples by timepoint
Mid<- subset_samples(prunedBeta_all, time !="End")
meta(Mid)

#Calculate distance and save as a matrix
set.seed(123)
all_bray_matrix <- phyloseq::distance(Mid, method = "bray")

#Betadisperser to see if the dispersion, variance, of two or more groups are significantly different or not. We will do the dispersion test for all the factors i.e innoculation, cultivar and time
#For innoculation
disprinnoc <- vegan::betadisper(all_bray_matrix, phyloseq::sample_data(Mid)$innoculation)
disprinnoc

plot(disprinnoc, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")

boxplot(disprinnoc, main = "", xlab = "")

permutest(disprinnoc)

#Similarly for cultivars
disprculti <- vegan::betadisper(all_bray_matrix, phyloseq::sample_data(Mid)$cultivar)
disprculti

plot(disprculti, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")

boxplot(disprculti, main = "", xlab = "")

permutest(disprculti)

#For time
disprtime <- vegan::betadisper(all_bray_matrix, phyloseq::sample_data(Mid)$environment)
disprtime

plot(disprtime, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")

boxplot(disprtime, main = "", xlab = "")

permutest(disprtime)

#The  ANOVA identified no significant differences between all the groups dispersion. So, we reject the hypothesis that these groups may have different dispersion.  Groups which have similar dispersion may still be significantly different in regards to their centroids, which will be tested using a PERMANOVA.

#Run PERMANOVA on distances with ADONIS test
adonis2(all_bray_matrix ~ phyloseq::sample_data(Mid)$cultivar+phyloseq::sample_data(Mid)$environment+phyloseq::sample_data(Mid)$innoculation,  permutations = 1000,   by = "margin")


#Effect of environment in innoculated samples (X10R)

X10Rinnoc <- subset_samples(Beta_all, cultivar != "ECW" & innoculation != "Control")
meta(X10Rinnoc)
#Generate distance matrix
innoc_bray_X10R <- phyloseq::distance(X10Rinnoc, method = "bray") 
#betadisper
disprinnocX10R <- vegan::betadisper(innoc_bray_X10R, phyloseq::sample_data(X10Rinnoc)$environment)
permutest(disprinnocX10R)
#ADONIS test
adonis2(innoc_bray_X10R ~ phyloseq::sample_data(X10Rinnoc)$environment+phyloseq::sample_data(X10Rinnoc)$time,   by = "margin")

#X10R control plants
X10Rcont <- subset_samples(Beta_all, cultivar != "ECW" & innoculation != "Innoculated")
meta(X10Rcont)
#Generate distance matrix
cont_bray_X10R <- phyloseq::distance(X10Rcont, method = "bray") 
#betadisper
disprcontX10R <- vegan::betadisper(cont_bray_X10R, phyloseq::sample_data(X10Rcont)$environment)
permutest(disprcontX10R)
#ADONIS test
adonis2(cont_bray_X10R ~ phyloseq::sample_data(X10Rcont)$environment+phyloseq::sample_data(X10Rcont)$time,   by = "margin")

#Now for ECW innoc
ECWinnoc <- subset_samples(Beta_all, cultivar != "X10R" & innoculation != "Control")
meta(ECWinnoc)
#Generate distance matrix
innoc_bray_ECW <- phyloseq::distance(ECWinnoc, method = "bray") 
#betadisper
disprinnocECW <- vegan::betadisper(innoc_bray_ECW, phyloseq::sample_data(ECWinnoc)$environment)
permutest(disprinnocECW)
#ADONIS test
adonis2(innoc_bray_ECW ~ phyloseq::sample_data(ECWinnoc)$environment+phyloseq::sample_data(ECWinnoc)$time,   by = "margin")

#Now for ECW control
ECWcont <- subset_samples(Beta_all, cultivar != "X10R" & innoculation != "Innoculated")
meta(ECWcont)
#Generate distance matrix
cont_bray_ECW <- phyloseq::distance(ECWcont, method = "bray") 
#betadisper
disprcontECW <- vegan::betadisper(cont_bray_ECW, phyloseq::sample_data(ECWcont)$environment)
permutest(disprcontECW)
#ADONIS test
adonis2(cont_bray_ECW ~ phyloseq::sample_data(ECWcont)$environment+phyloseq::sample_data(ECWcont)$time,   by = "margin")

```


```{r}
#to look at the effect of different factors in mid and end season
Mid<- subset_samples(Beta_all, time !="Mid")
meta(Mid)
sample_sums(Mid)
sort(sample_sums(Mid))

#SRS
minControl= 2548613
keep=names(which(sample_sums(Mid)>=minControl))
prunedMid_all=prune_samples(keep, Mid)
any(taxa_sums(prunedMid_all) == 0)
meta(prunedMid_all)
prunedBeta_all <- prune_taxa(taxa_sums(prunedMid_all) > 0, prunedMid_all)
any(taxa_sums(prunedBeta_all) == 0)
meta(prunedBeta_all)

library(SRS)

pruned_data_df = as.data.frame(otu_table(prunedBeta_all))

SRS_OUTPUT <- SRS (pruned_data_df,minControl)

rownames(SRS_OUTPUT)=rownames(pruned_data_df)

table(sample_data(prunedBeta_all)[,"HET"])

# transform back into phyloseq object
taxa = otu_table(SRS_OUTPUT, taxa_are_rows = TRUE)
otu_table(prunedBeta_all)=taxa
any(taxa_sums(prunedBeta_all) == 0)

prunedBeta_all <- prune_taxa(taxa_sums(prunedBeta_all) > 0, prunedBeta_all)
any(taxa_sums(prunedBeta_all) == 0)

SeqDepthPruned = sample_sums(prunedBeta_all)
sample_data(prunedBeta_all)$SeqDepthPruned = SeqDepthPruned

# barplot of library sizes
library(ggplot2)
ggplot(meta(prunedBeta_all), aes(id, SeqDepthPruned)) + geom_bar(stat = "identity", aes(fill = environment)) +
  rotate_x_text()

prunedBeta_all

#Calculate distance and save as a matrix
set.seed(123)
all_bray_matrix <- phyloseq::distance(prunedBeta_all, method = "bray")

#Betadisperser to see if the dispersion, variance, of two or more groups are significantly different or not. We will do the dispersion test for all the factors i.e innoculation, cultivar and time
#For innoculation
disprinnoc <- vegan::betadisper(all_bray_matrix, phyloseq::sample_data(prunedBeta_all)$innoculation)
disprinnoc

plot(disprinnoc, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")

boxplot(disprinnoc, main = "", xlab = "")

permutest(disprinnoc)

#Similarly for cultivars
disprculti <- vegan::betadisper(all_bray_matrix, phyloseq::sample_data(prunedBeta_all)$cultivar)
disprculti

plot(disprculti, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")

boxplot(disprculti, main = "", xlab = "")

permutest(disprculti)

#For time
disprtime <- vegan::betadisper(all_bray_matrix, phyloseq::sample_data(prunedBeta_all)$environment)
disprtime

plot(disprtime, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")

boxplot(disprtime, main = "", xlab = "")

permutest(disprtime)

#The  ANOVA identified no significant differences between all the groups dispersion. So, we reject the hypothesis that these groups may have different dispersion.  Groups which have similar dispersion may still be significantly different in regards to their centroids, which will be tested using a PERMANOVA.
meta(prunedBeta_all)
#Run PERMANOVA on distances with ADONIS test
adonis2(all_bray_matrix ~ phyloseq::sample_data(prunedBeta_all)$cultivar+phyloseq::sample_data(prunedBeta_all)$environment+phyloseq::sample_data(prunedBeta_all)$innoculation,  permutations = 1000,   by = "margin")



ano = anosim(all_bray_matrix, prunedBeta_all$cultivar,  permutations = 9999)
ano

```

