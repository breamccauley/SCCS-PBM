library(ape)
library(brms)
library(dplyr)
library(tibble)
library(tidyverse)

## TIGHT SAMPLE

#Read in file of SCCS supertree
SCCS_Tree<-read.tree("~/SCCS_Supertree.tre")

#Ladderize SCCS_Tree
SCCS_Tree<-ladderize(SCCS_Tree, right = TRUE)

# Read data into R
SCCS_PBM_PA <- read.csv("~/SCCS_PBM_Data_Tight.csv")

# Replace 'nd' cells with NA
SCCS_PBM_PA<- replace(SCCS_PBM_PA, SCCS_PBM_PA=='nd', NA)

#Make sure no duplicates
SCCS_PBM_PA = SCCS_PBM_PA[!duplicated(SCCS_PBM_PA$Node),]

# Set row names as 'node' column
rownames(SCCS_PBM_PA) = SCCS_PBM_PA$Node

#Avoid zero edge lengths
SCCS_Tree$edge.length[SCCS_Tree$edge.length==0] = 0.000001

#Order the data frame by the order of tree tiplabels
SCCS_PBM_PA <- SCCS_PBM_PA[SCCS_Tree$tip.label,]

#Ladderize sample
SCCS_Tree <- ladderize(SCCS_Tree, right = TRUE)

# Make sure that the data is numeric
SCCS_PBM_PA$Total_Types = as.numeric(SCCS_PBM_PA$Total_Types)
SCCS_PBM_PA$Total_Types_M = as.numeric(SCCS_PBM_PA$Total_Types_M)
SCCS_PBM_PA$Total_Types_F = as.numeric(SCCS_PBM_PA$Total_Types_F)

SCCS_PBM_PA$Total_Subtypes = as.numeric(SCCS_PBM_PA$Total_Subtypes)
SCCS_PBM_PA$Total_M_Subtypes = as.numeric(SCCS_PBM_PA$Total_M_Subtypes)
SCCS_PBM_PA$Total_F_Subtypes = as.numeric(SCCS_PBM_PA$Total_F_Subtypes)

## also scale the predictors
SCCS_PBM_PA <- SCCS_PBM_PA %>%
  mutate(
    Sources_Num_z = as.numeric(scale(Sources_Num)),
    Pages_Num_z = as.numeric(scale(Pages_Num)),
    PBMSources_Num_z = as.numeric(scale(PBMSources_Num)),
    PBMPages_Num_z = as.numeric(scale(PBMPages_Num))
  )

#create covariance matrix
covarMatrix <- vcv.phylo(SCCS_Tree)

# Run a model predicting number of PBM types by documentation
Types_AllSources <- brm(
  Total_Types ~ Sources_Num_z + Pages_Num_z + PBMSources_Num_z + PBMPages_Num_z +
    (1|gr(Node, cov = CX)), # Line to add the phylogenetic control (links to 'CX' below)
  data = SCCS_PBM_PA,  # main trait data frame
  family = poisson(), 
  data2 = list(CX = covarMatrix), # The secondary data - the covariance matrix
  control = list(adapt_delta = 0.95)) 

saveRDS(Types_AllSources, "~/Types_AllSources.rds")

# Run a model predicting number of PBM types M by documentation
Types_M_AllSources <- brm(
  Total_Types_M ~ Sources_Num_z + Pages_Num_z + PBMSources_Num_z + PBMPages_Num_z +
    (1|gr(Node, cov = CX)), # Line to add the phylogenetic control (links to 'CX' below)
  data = SCCS_PBM_PA,  # main trait data frame
  family = poisson(), 
  data2 = list(CX = covarMatrix), # The secondary data - the covariance matrix
  control = list(adapt_delta = 0.95)) 

saveRDS(Types_M_AllSources, "~/Types_M_AllSources.rds")

# Run a model predicting number of PBM types F by documentation
Types_F_AllSources <- brm(
  Total_Types_F ~ Sources_Num_z + Pages_Num_z + PBMSources_Num_z + PBMPages_Num_z +
    (1|gr(Node, cov = CX)), # Line to add the phylogenetic control (links to 'CX' below)
  data = SCCS_PBM_PA,  # main trait data frame
  family = poisson(), 
  data2 = list(CX = covarMatrix), # The secondary data - the covariance matrix
  control = list(adapt_delta = 0.95)) 

saveRDS(Types_F_AllSources, "~/Types_F_AllSources.rds")

# Run a model predicting number of PBM subtypes by documentation
Subtypes_AllSources <- brm(
  Total_Subtypes ~ Sources_Num_z + Pages_Num_z + PBMSources_Num_z + PBMPages_Num_z +
    (1|gr(Node, cov = CX)), # Line to add the phylogenetic control (links to 'CX' below)
  data = SCCS_PBM_PA,  # main trait data frame
  family = poisson(), 
  data2 = list(CX = covarMatrix), # The secondary data - the covariance matrix
  control = list(adapt_delta = 0.95)) 

saveRDS(Subtypes_AllSources, "~/Subtypes_AllSources.rds")

# Run a model predicting number of male PBM subtypes by documentation
Subtypes_M_AllSources <- brm(
  Total_M_Subtypes ~ Sources_Num_z + Pages_Num_z + PBMSources_Num_z + PBMPages_Num_z +
    (1|gr(Node, cov = CX)), # Line to add the phylogenetic control (links to 'CX' below)
  data = SCCS_PBM_PA,  # main trait data frame
  family = poisson(), 
  data2 = list(CX = covarMatrix), # The secondary data - the covariance matrix
  control = list(adapt_delta = 0.95)) 

saveRDS(Subtypes_M_AllSources, "~/Subtypes_M_AllSources.rds")

# Run a model predicting number of female PBM subtypes by documentation
Subtypes_F_AllSources <- brm(
  Total_F_Subtypes ~ Sources_Num_z + Pages_Num_z + PBMSources_Num_z + PBMPages_Num_z +
    (1|gr(Node, cov = CX)), # Line to add the phylogenetic control (links to 'CX' below)
  data = SCCS_PBM_PA,  # main trait data frame
  family = poisson(), 
  data2 = list(CX = covarMatrix), # The secondary data - the covariance matrix
  control = list(adapt_delta = 0.95)) 

saveRDS(Subtypes_F_AllSources, "~/Subtypes_F_AllSources.rds")

## LOOSE SAMPLE

#Read in file of SCCS supertree
SCCS_Tree<-read.tree("~/SCCS_Supertree.tre")

#Ladderize SCCS_Tree
SCCS_Tree<-ladderize(SCCS_Tree, right = TRUE)

# Read data into R (loose SCCS sample)
SCCS_PBM_PA_Loose <- read.csv("~/SCCS_PBM_Data_Loose.csv")

# Replace 'nd' cells with NA
SCCS_PBM_PA_Loose<- replace(SCCS_PBM_PA_Loose, SCCS_PBM_PA_Loose=='nd', NA)

#Make sure no duplicates
SCCS_PBM_PA_Loose = SCCS_PBM_PA_Loose[!duplicated(SCCS_PBM_PA_Loose$Node),]

# Set row names as 'node' column
rownames(SCCS_PBM_PA_Loose) = SCCS_PBM_PA_Loose$Node

#Avoid zero edge lengths
SCCS_Tree$edge.length[SCCS_Tree$edge.length==0] = 0.000001

#Order the data frame by the order of tree tiplabels
SCCS_PBM_PA_Loose <- SCCS_PBM_PA_Loose[SCCS_Tree$tip.label,]

#Ladderize sample
SCCS_Tree <- ladderize(SCCS_Tree, right = TRUE)

# Make sure that the data is numeric
SCCS_PBM_PA_Loose$Total_Types = as.numeric(SCCS_PBM_PA_Loose$Total_Types)
SCCS_PBM_PA_Loose$Total_Types_M = as.numeric(SCCS_PBM_PA_Loose$Total_Types_M)
SCCS_PBM_PA_Loose$Total_Types_F = as.numeric(SCCS_PBM_PA_Loose$Total_Types_F)

SCCS_PBM_PA_Loose$Total_Subtypes = as.numeric(SCCS_PBM_PA_Loose$Total_Subtypes)
SCCS_PBM_PA_Loose$Total_M_Subtypes = as.numeric(SCCS_PBM_PA_Loose$Total_M_Subtypes)
SCCS_PBM_PA_Loose$Total_F_Subtypes = as.numeric(SCCS_PBM_PA_Loose$Total_F_Subtypes)

## also scale the predictors
SCCS_PBM_PA_Loose <- SCCS_PBM_PA_Loose %>%
  mutate(
    Sources_Num_z = as.numeric(scale(Sources_Num)),
    Pages_Num_z = as.numeric(scale(Pages_Num)),
    PBMSources_Num_z = as.numeric(scale(PBMSources_Num)),
    PBMPages_Num_z = as.numeric(scale(PBMPages_Num))
  )

#create covariance matrix
covarMatrix <- vcv.phylo(SCCS_Tree)

# Run a model predicting number of PBM types by documentation
Types_Loose_AllSources <- brm(
  Total_Types ~ Sources_Num_z + Pages_Num_z + PBMSources_Num_z + PBMPages_Num_z +
    (1|gr(Node, cov = CX)), # Line to add the phylogenetic control (links to 'CX' below)
  data = SCCS_PBM_PA_Loose,  # main trait data frame
  family = poisson(), 
  data2 = list(CX = covarMatrix), # The secondary data - the covariance matrix
  control = list(adapt_delta = 0.95)) 

saveRDS(Types_Loose_AllSources, "~/Types_Loose_AllSources.rds")

# Run a model predicting number of PBM types M by documentation
Types_M_Loose_AllSources <- brm(
  Total_Types_M ~ Sources_Num_z + Pages_Num_z + PBMSources_Num_z + PBMPages_Num_z +
    (1|gr(Node, cov = CX)), # Line to add the phylogenetic control (links to 'CX' below)
  data = SCCS_PBM_PA_Loose,  # main trait data frame
  family = poisson(), 
  data2 = list(CX = covarMatrix), # The secondary data - the covariance matrix
  control = list(adapt_delta = 0.95))

saveRDS(Types_M_Loose_AllSources, "~/Types_M_Loose_AllSources.rds")

# Run a model predicting number of PBM types F by documentation
Types_F_Loose_AllSources <- brm(
  Total_Types_F ~ Sources_Num_z + Pages_Num_z + PBMSources_Num_z + PBMPages_Num_z +
    (1|gr(Node, cov = CX)), # Line to add the phylogenetic control (links to 'CX' below)
  data = SCCS_PBM_PA_Loose,  # main trait data frame
  family = poisson(), 
  data2 = list(CX = covarMatrix), # The secondary data - the covariance matrix
  control = list(adapt_delta = 0.95))

saveRDS(Types_F_Loose_AllSources, "~/Types_F_Loose_AllSources.rds")

# Run a model predicting number of PBM subtypes by documentation
Subtypes_Loose_AllSources <- brm(
  Total_Subtypes ~ Sources_Num_z + Pages_Num_z + PBMSources_Num_z + PBMPages_Num_z +
    (1|gr(Node, cov = CX)), # Line to add the phylogenetic control (links to 'CX' below)
  data = SCCS_PBM_PA_Loose,  # main trait data frame
  family = poisson(), 
  data2 = list(CX = covarMatrix), # The secondary data - the covariance matrix
  control = list(adapt_delta = 0.95)) 

saveRDS(Subtypes_Loose_AllSources, "~/Subtypes_Loose_AllSources.rds")

# Run a model predicting number of male PBM subtypes by documentation
Subtypes_M_Loose_AllSources <- brm(
  Total_M_Subtypes ~ Sources_Num_z + Pages_Num_z + PBMSources_Num_z + PBMPages_Num_z +
    (1|gr(Node, cov = CX)), # Line to add the phylogenetic control (links to 'CX' below)
  data = SCCS_PBM_PA_Loose,  # main trait data frame
  family = poisson(), 
  data2 = list(CX = covarMatrix), # The secondary data - the covariance matrix
  control = list(adapt_delta = 0.95)) 

saveRDS(Subtypes_M_Loose_AllSources, "~/Subtypes_M_Loose_AllSources.rds")

# Run a model predicting number of female PBM subtypes by documentation
Subtypes_F_Loose_AllSources <- brm(
  Total_F_Subtypes ~ Sources_Num_z + Pages_Num_z + PBMSources_Num_z + PBMPages_Num_z +
    (1|gr(Node, cov = CX)), # Line to add the phylogenetic control (links to 'CX' below)
  data = SCCS_PBM_PA_Loose,  # main trait data frame
  family = poisson(), 
  data2 = list(CX = covarMatrix), # The secondary data - the covariance matrix
  control = list(adapt_delta = 0.95)) 

saveRDS(Subtypes_F_Loose_AllSources, "~/Subtypes_F_Loose_AllSources.rds")

