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
SCCS_PBM_PA$Total_Types_M = as.numeric(SCCS_PBM_PA$Total_Types_M, na.rm=TRUE)
SCCS_PBM_PA$Total_Types_F = as.numeric(SCCS_PBM_PA$Total_Types_F, na.rm=TRUE)

SCCS_PBM_PA$Total_M_Subtypes = as.numeric(SCCS_PBM_PA$Total_M_Subtypes, na.rm=TRUE)
SCCS_PBM_PA$Total_F_Subtypes = as.numeric(SCCS_PBM_PA$Total_F_Subtypes, na.rm=TRUE)

## scale the predictors
SCCS_PBM_PA <- SCCS_PBM_PA %>%
  mutate(
    PBMPages_Num_z = as.numeric(scale(PBMPages_Num))
  )

#create covariance matrix
covarMatrix <- vcv.phylo(SCCS_Tree)

# Set weak priors
priors <- c(
  prior(normal(0, 1), class = "b"),
  prior(normal(0, 2), class = "Intercept"),
  prior(exponential(2), class = "sd")
)

# Types
## Each female type is associated with how many increased male types?
Types_MF_subset <- SCCS_PBM_PA %>%
  select(Node,
         Total_Types_M, Total_Types_F, PBMPages_Num_z) 

## male dependent, female independent
Types_MF <- brm(
  Total_Types_M ~ Total_Types_F + PBMPages_Num_z +
    (1 | gr(Node, cov = CX)),        # phylogenetic covariance
  data = Types_MF_subset,
  family = poisson(),
  prior = priors,
  chains = 4,
  iter = 8000,
  data2 = list(CX = covarMatrix),
  control = list(adapt_delta = 0.995, max_treedepth = 15)
)

saveRDS(Types_MF, "~/Types_MF_Tight.rds")

## female dependent, male independent
Types_FM <- brm(
  Total_Types_F ~ Total_Types_M + PBMPages_Num_z +
    (1 | gr(Node, cov = CX)),        # phylogenetic covariance
  data = Types_MF_subset,
  family = poisson(),
  prior = priors,
  chains = 4,
  iter = 8000,
  data2 = list(CX = covarMatrix),
  control = list(adapt_delta = 0.995, max_treedepth = 15)
)

saveRDS(Types_FM, "~/Types_FM_Tight.rds")

# Subtypes
Subtypes_MF_subset <- SCCS_PBM_PA %>%
  select(Node,
         Total_M_Subtypes, Total_F_Subtypes, PBMPages_Num_z) 

## male dependent, female independent
Subtypes_MF <- brm(
  Total_M_Subtypes ~ Total_F_Subtypes + PBMPages_Num_z +
    (1 | gr(Node, cov = CX)),        # phylogenetic covariance
  data = Subtypes_MF_subset,
  family = poisson(),
  prior = priors,
  chains = 4,
  iter = 8000,
  data2 = list(CX = covarMatrix),
  control = list(adapt_delta = 0.995, max_treedepth = 15)
)

saveRDS(Subtypes_MF, "~/Subtypes_MF_Tight.rds")

## female dependent, male independent
Subtypes_FM <- brm(
  Total_F_Subtypes ~ Total_M_Subtypes + PBMPages_Num_z +
    (1 | gr(Node, cov = CX)),        # phylogenetic covariance
  data = Subtypes_MF_subset,
  family = poisson(),
  prior = priors,
  chains = 4,
  iter = 10000,
  data2 = list(CX = covarMatrix),
  control = list(adapt_delta = 0.995, max_treedepth = 15)
)

saveRDS(Subtypes_FM, "~/Subtypes_FM_Tight.rds")

## LOOSE SAMPLE
#Read in file of SCCS supertree
SCCS_Tree<-read.tree("~/SCCS_Supertree.tre")

#Ladderize SCCS_Tree
SCCS_Tree<-ladderize(SCCS_Tree, right = TRUE)

# Read data into R
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
SCCS_PBM_PA_Loose$Total_Types_M = as.numeric(SCCS_PBM_PA_Loose$Total_Types_M, na.rm=TRUE)
SCCS_PBM_PA_Loose$Total_Types_F = as.numeric(SCCS_PBM_PA_Loose$Total_Types_F, na.rm=TRUE)

SCCS_PBM_PA_Loose$Total_M_Subtypes = as.numeric(SCCS_PBM_PA_Loose$Total_M_Subtypes, na.rm=TRUE)
SCCS_PBM_PA_Loose$Total_F_Subtypes = as.numeric(SCCS_PBM_PA_Loose$Total_F_Subtypes, na.rm=TRUE)

## scale the predictors
SCCS_PBM_PA_Loose <- SCCS_PBM_PA_Loose %>%
  mutate(
    PBMPages_Num_z = as.numeric(scale(PBMPages_Num))
  )

#create covariance matrix
covarMatrix <- vcv.phylo(SCCS_Tree)

# Set weak priors
priors <- c(
  prior(normal(0, 1), class = "b"),
  prior(normal(0, 2), class = "Intercept"),
  prior(exponential(2), class = "sd")
)

# Types
Types_MF_Loose_subset <- SCCS_PBM_PA_Loose %>%
  select(Node,
         Total_Types_M, Total_Types_F, PBMPages_Num_z) 

## male dependent, female independent
Types_MF_Loose <- brm(
  Total_Types_M ~ Total_Types_F + PBMPages_Num_z +
    (1 | gr(Node, cov = CX)),        # phylogenetic covariance
  data = Types_MF_Loose_subset,
  family = poisson(),
  prior = priors,
  chains = 4,
  iter = 8000,
  data2 = list(CX = covarMatrix),
  control = list(adapt_delta = 0.995, max_treedepth = 15)
)

saveRDS(Types_MF_Loose, "~/Types_MF_Loose.rds")

## female dependent, male independent
Types_FM_Loose <- brm(
  Total_Types_F ~ Total_Types_M + PBMPages_Num_z +
    (1 | gr(Node, cov = CX)),        # phylogenetic covariance
  data = Types_MF_Loose_subset,
  family = poisson(),
  prior = priors,
  chains = 4,
  iter = 8000,
  backend = "cmdstanr",
  data2 = list(CX = covarMatrix),
  control = list(adapt_delta = 0.995, max_treedepth = 15)
)

saveRDS(Types_FM_Loose, "~/Types_FM_Loose.rds")

# Subtypes
Subtypes_MF_Loose_subset <- SCCS_PBM_PA_Loose %>%
  select(Node,
         Total_M_Subtypes, Total_F_Subtypes, PBMPages_Num_z) 

## male dependent, female independent
Subtypes_MF_Loose <- brm(
  Total_M_Subtypes ~ Total_F_Subtypes + PBMPages_Num_z +
    (1 | gr(Node, cov = CX)),        # phylogenetic covariance
  data = Subtypes_MF_Loose_subset,
  family = poisson(),
  prior = priors,
  chains = 4,
  iter = 8000,
  backend = "cmdstanr",
  data2 = list(CX = covarMatrix),
  control = list(adapt_delta = 0.995, max_treedepth = 15)
)

saveRDS(Subtypes_MF_Loose, "~/Subtypes_MF_Loose.rds")

## female dependent, male independent
Subtypes_FM_Loose <- brm(
  Total_F_Subtypes ~ Total_M_Subtypes + PBMPages_Num_z +
    (1 | gr(Node, cov = CX)),        # phylogenetic covariance
  data = Subtypes_MF_Loose_subset,
  family = poisson(),
  prior = priors,
  chains = 4,
  iter = 8000,
  backend = "cmdstanr",
  data2 = list(CX = covarMatrix),
  control = list(adapt_delta = 0.995, max_treedepth = 15)
)

saveRDS(Subtypes_FM_Loose, "~/Subtypes_FM_Loose.rds")

