# Modified bNTI script for the Shale S3 (and probably other) datasets

# Modified from https://github.com/stegen/Stegen_etal_ISME_2013/blob/master/bNTI_Local_Machine.r; Stegen et al, 2013
# Modified from https://github.com/stegen/Stegen_etal_ISME_2013/blob/master/Raup_Crick_Abundance.r; Stegen et al, 2013
# RED 2019; danczak.6@osu.edu

#install.packages(c("vegan", "picante", "digest", "dplyr", "plyr", "tidyverse", "abind"))

library(vegan)
library(picante)
library(digest)
library(dplyr)
library(plyr)
library(tidyverse)
library(abind)

###################################
#### Data Loading and cleaning ####
###################################

setwd('~/Dropbox/PMI_3_analyses/multi-omics_data/amplicon/16S/ecology_tables/just_soil-hip/bnti') # Working directory for 16S comm on my local computer

day0 = read.csv("../Day0_CMU.csv") # Importing the organismal data
early = read.csv("../Early_CMU.csv") # Importing the organismal data
active = read.csv("../Active_CMU.csv") # Importing the organismal
advanced = read.csv("../Advanced_CMU.csv") # Importing the organismal
temp1 = full_join(day0,early, by='X')
temp2 = full_join(temp1,active, by='X')
data = full_join(temp2,advanced, by='X')
rownames(data) = data[,1]
data[,1] = NULL
data[is.na(data)] = 0
tree = read.tree("../../tree2.nwk") # Importing the tree

####################################
#### Beginning the bNTI Process ####
####################################

# Matching the tree to the newly rarefied OTU dataset
phylo = match.phylo.data(tree, data)

# Calculating bMNTD for my samples
#bMNTD = as.matrix(comdistnt(t(phylo$data), cophenetic(phylo$phy), abundance.weighted = T))
#saveRDS(bMNTD, "CMU_bMNTD.rds")
bMNTD = readRDS(file="CMU_bMNTD.rds")

# Calculating the bMNTD for 999 random distributions
#bMNTD.rand = array(c(-999), dim = c(ncol(phylo$data), ncol(phylo$data), 999)) # Creating 999 'dummy' matrices to put random results into

#for(i in 1:999){
#  bMNTD.rand[,,i] = as.matrix(comdistnt(t(phylo$data), taxaShuffle(cophenetic(phylo$phy)), abundance.weighted = T, exclude.conspecifics = F))
#  print(c(date(),i)) # Measuring how far the loop has gotten
#} # Performing the calculations on using the OTU table but with randomized taxonomic affiliations
#saveRDS(bMNTD.rand, "UTK_bMNTD_rand1.rds")

# load and combine random calculations
CMU_rand1 = readRDS(file="CMU_bMNTD_rand1.rds")
CMU_rand2 = readRDS(file="CMU_bMNTD_rand2.rds")
CMU_rand3 = readRDS(file="CMU_bMNTD_rand3.rds")
CMU_rand4 = readRDS(file="CMU_bMNTD_rand4.rds")
bMNTD.rand = abind(CMU_rand1,CMU_rand2,CMU_rand3,CMU_rand4)

# sanity check that random arrays combined and match sample array
(sum(dim(CMU_rand1)[3],dim(CMU_rand2)[3],dim(CMU_rand3)[3],dim(CMU_rand4)[3]) + sum(dim(bMNTD.rand)[3])) / 2 == 999
dim(CMU_rand1)[1] == dim(CMU_rand2)[1] | dim(CMU_rand2)[1] == dim(CMU_rand3)[1] | dim(CMU_rand3)[1] == dim(CMU_rand4)[1]
dim(bMNTD.rand)[1] == dim(bMNTD)[1]

# Calculating the bNTI between these randomized communities and the empirical dataset
bNTI = matrix(c(NA), nrow = ncol(phylo$data), ncol = ncol(phylo$data))

for(i in 1:(ncol(phylo$data)-1)){ 
  for(j in (i+1):ncol(phylo$data)){
    m = bMNTD.rand[j,i,] # Just setting all the randomizations for a given comparison to a matrix
    bNTI[j,i] = ((bMNTD[j,i]-mean(m))/sd(m)) # The bNTI calculation
  }
}

rownames(bNTI) = colnames(phylo$data)
colnames(bNTI) = colnames(phylo$data)
bNTI
hist(bNTI)

write.csv(bNTI, "CMU_Community_bNTI_TotalCounts.csv", quote = F)

