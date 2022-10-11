#Network analysis is an RAM intensive job so we submitted the codes in HPC.

#/bin/bash
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load R/4.1.0
R CMD BATCH inference.R2

#The inference.R2 file contains the following R scripts.

library(devtools)
library(SpiecEasi)
library(microbiome)
library(phyloseq)
library(NetCoMi)
# Load  data
net.c1 <- readRDS("/scratch/aubrrb/netcomi/allsamples.rds")
net.c<- subset_samples(net.c1, environment !="GH" & environment !="Ozone")

X10Rinnocamb_prevalent <- core(net.c, detection = 1^2, prevalence = 20/100)

net.c <- prune_taxa(taxa_sums(X10Rinnocamb_prevalent) >100, X10Rinnocamb_prevalent)
# Extract count matrix and phenotypes
env_counts <- t(net.c@otu_table@.Data)
colnames(env_counts) <- net.c@tax_table@.Data[, "Genus"]
dim(env_counts)



env_pheno <- net.c@sam_data@.Data
names(env_pheno) <- net.c@sam_data@names

environment <- env_pheno$innoculation

#construct net with higherst freq more than 80
net_env <- netConstruct(env_counts, group = environment, filtTax = "highestFreq",
                           filtTaxPar = list(highestFreq = 80),
                           zeroMethod = "none", normMethod = "none",
                           measure = "spring",
                           measurePar = list(nlambda = 50, rep.num = 50,
                                             ncores=30),
                           seed = 20190101)

netprops_env <- netAnalyze(net_env, clustMethod = "cluster_fast_greedy",
                              hubPar = "eigenvector", hubQuant = 0.95,
                              lnormFit = FALSE)

summary(netprops_env, groupNames = c("Control", "Innoculated"))



pdf("./biotic.pdf", width = 50, height = 25)
p <- plot(netprops_env,
          sameLayout = TRUE,
          layoutGroup = 1,
          repulsion = 1,
          groupNames = c("Control","Innoculated"),
          shortenLabels = "intelligent",
          labelLength = 14,
          labelPattern = c(14,"'",3),
          charToRm = "g__",
          labelScale = FALSE,
          nodeFilter = "none",
          nodeFilterPar = 20,
          rmSingles = "inboth",
          nodeSize = "eigen",
          nodeSizeSpread = 3,
          nodeColor = "cluster",
          colorVec = rainbow(12),
          nodeTransp = 65,
          hubTransp = 50,
          hubBorderWidth = 2,
          hubBorderCol = "gray40",
          colorNegAsso = TRUE,
          edgeWidth = 1,
          edgeTranspLow = 70,
          edgeTranspHigh = 30,
          cexNodes = 1,
          cexHubs = 1.3,
          cexTitle = 2.0,
          cexLabels = 0.8,
          showTitle = TRUE,
          mar = c(1,1,3,1))
dev.off()

# Run on server with 15 CPU cores:
net_env_comp<- netCompare(netprops_env, permTest = TRUE,
                             lnormFit = FALSE, jaccQuant = 0.75,
                             nPerm = 1000, cores = 30,
                             seed = 20190101, adjust = "none")

# Rerun with multiple testing adjustment
net_env_comp_adaptbh <- netCompare(netprops_env, permTest = TRUE,
                                      lnormFit = FALSE, jaccQuant = 0.75,
                                      nPerm = 1000, cores = 30,
                                      seed = 20190101, adjust = "adaptBH",
                                      assoPerm = net_env_comp$assoPerm)

summary(net_env_comp, pAdjust = TRUE,
        groupNames = c("Control", "Innoculated"), digitsPval = 6)

summary(net_env_comp_adaptbh, pAdjust = TRUE,
        groupNames = c("Control", "Innoculated"), digitsPval = 6)
