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
SCCS_PBM_PA$Total_Types = as.numeric(SCCS_PBM_PA$Total_Types, na.rm=TRUE)
SCCS_PBM_PA$Total_Types_M = as.numeric(SCCS_PBM_PA$Total_Types_M, na.rm=TRUE)
SCCS_PBM_PA$Total_Types_F = as.numeric(SCCS_PBM_PA$Total_Types_F, na.rm=TRUE)

SCCS_PBM_PA$Total_Subtypes = as.numeric(SCCS_PBM_PA$Total_Subtypes, na.rm=TRUE)
SCCS_PBM_PA$Total_M_Subtypes = as.numeric(SCCS_PBM_PA$Total_M_Subtypes, na.rm=TRUE)
SCCS_PBM_PA$Total_F_Subtypes = as.numeric(SCCS_PBM_PA$Total_F_Subtypes, na.rm=TRUE)

## also scale the predictors
SCCS_PBM_PA <- SCCS_PBM_PA %>%
  mutate(
    PBMPages_Num_z = as.numeric(scale(PBMPages_Num))
  )

#create covariance matrix
covarMatrix <- vcv.phylo(SCCS_Tree)

# Sum-code geographic region
SCCS_PBM_PA <- SCCS_PBM_PA %>%
  mutate(
    Geographic_Region = factor(Geographic_Region)
  )

contrasts(SCCS_PBM_PA$Geographic_Region) <- contr.sum(
  nlevels(SCCS_PBM_PA$Geographic_Region)
)

levels(SCCS_PBM_PA$Geographic_Region)

# Set weak priors
priors <- c(
  prior(normal(0, 1), class = "b"),
  prior(normal(0, 2), class = "Intercept"),
  prior(exponential(3), class = "sd")
)

# run PBM types and geographic region
Types_GeoRegion_Pages <- brm(
  Total_Types ~ Geographic_Region + PBMPages_Num_z +
   (1 | gr(Node, cov = CX)),        # phylogenetic covariance
  data = SCCS_PBM_PA,
  family = poisson(),
  chains = 4,
  iter = 15000,
  init = 0,
  data2 = list(CX = covarMatrix),
  control = list(adapt_delta = 0.99, max_treedepth = 20)
)

saveRDS(Types_GeoRegion_Pages, "~/Types_GeoRegion_Pages.rds")

## Males and females
# Males
Types_M_GeoRegion_Pages <- brm(
  Total_Types_M ~ Geographic_Region + PBMPages_Num_z +
    (1 | gr(Node, cov = CX)),        # phylogenetic covariance
  data = SCCS_PBM_PA,
  family = poisson(),
  prior = priors,
  chains = 4,
  iter = 15000,
  init = 0,
  data2 = list(CX = covarMatrix),
  control = list(adapt_delta = 0.995, max_treedepth = 20)
)

saveRDS(Types_M_GeoRegion_Pages, "~/Types_M_GeoRegion_Pages.rds")

# Females
Types_F_GeoRegion_Pages <- brm(
  Total_Types_F ~ Geographic_Region + PBMPages_Num_z +
    (1 | gr(Node, cov = CX)),        # phylogenetic covariance
  data = SCCS_PBM_PA,
  family = poisson(),
  prior = priors,
  chains = 4,
  iter = 15000,
  init = 0,
  data2 = list(CX = covarMatrix),
  control = list(adapt_delta = 0.995, max_treedepth = 20)
)

saveRDS(Types_F_GeoRegion_Pages, "~/Types_F_GeoRegion_Pages.rds")

# run PBM subtypes and geographic region
Subtypes_GeoRegion_Pages <- brm(
 Total_Subtypes ~ Geographic_Region + PBMPages_Num_z +
   (1 | gr(Node, cov = CX)),        # phylogenetic covariance
 data = SCCS_PBM_PA,
 family = poisson(),
 data2 = list(CX = covarMatrix),
 chains = 4,
 init = 0,
 iter = 15000,
 control = list(adapt_delta = 0.995, max_treedepth = 20)
)

saveRDS(Subtypes_GeoRegion_Pages, "~/Subtypes_GeoRegion_Pages.rds")

## Males and Females
# Males
Subtypes_M_GeoRegion_Pages <- brm(
  Total_M_Subtypes ~ Geographic_Region + PBMPages_Num_z +
    (1 | gr(Node, cov = CX)),        # phylogenetic covariance
  data = SCCS_PBM_PA,
  family = poisson(),
  prior = priors,
  data2 = list(CX = covarMatrix),
  chains = 4,
  iter = 15000,
  init = 0,
  control = list(adapt_delta = 0.995, max_treedepth = 20)
)

saveRDS(Subtypes_M_GeoRegion_Pages, "~/Subtypes_M_GeoRegion_Pages.rds")

# Females
Subtypes_F_GeoRegion_Pages <- brm(
  Total_F_Subtypes ~ Geographic_Region + PBMPages_Num_z +
    (1 | gr(Node, cov = CX)),        # phylogenetic covariance
  data = SCCS_PBM_PA,
  family = poisson(),
  prior = priors,
  data2 = list(CX = covarMatrix),
  chains = 4,
  iter = 15000,
  init = 0,
  control = list(adapt_delta = 0.995, max_treedepth = 20)
)

saveRDS(Subtypes_F_GeoRegion_Pages, "~/Subtypes_F_GeoRegion_Pages.rds")

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
SCCS_PBM_PA_Loose$Total_Types = as.numeric(SCCS_PBM_PA_Loose$Total_Types, na.rm=TRUE)
SCCS_PBM_PA_Loose$Total_Types_M = as.numeric(SCCS_PBM_PA_Loose$Total_Types_M, na.rm=TRUE)
SCCS_PBM_PA_Loose$Total_Types_F = as.numeric(SCCS_PBM_PA_Loose$Total_Types_F, na.rm=TRUE)

SCCS_PBM_PA_Loose$Total_Subtypes = as.numeric(SCCS_PBM_PA_Loose$Total_Subtypes, na.rm=TRUE)
SCCS_PBM_PA_Loose$Total_M_Subtypes = as.numeric(SCCS_PBM_PA_Loose$Total_M_Subtypes, na.rm=TRUE)
SCCS_PBM_PA_Loose$Total_F_Subtypes = as.numeric(SCCS_PBM_PA_Loose$Total_F_Subtypes, na.rm=TRUE)

## also scale the predictors
SCCS_PBM_PA_Loose <- SCCS_PBM_PA_Loose %>%
  mutate(
    PBMPages_Num_z = as.numeric(scale(PBMPages_Num))
  )

#create covariance matrix
covarMatrix <- vcv.phylo(SCCS_Tree)

# Sum-code geographic region
SCCS_PBM_PA_Loose <- SCCS_PBM_PA_Loose %>%
  mutate(
    Geographic_Region = factor(Geographic_Region)
  )

contrasts(SCCS_PBM_PA_Loose$Geographic_Region) <- contr.sum(
  nlevels(SCCS_PBM_PA_Loose$Geographic_Region)
)

# Set weak priors
priors <- c(
  prior(normal(0, 1), class = "b"),
  prior(normal(0, 2), class = "Intercept"),
  prior(exponential(3), class = "sd")
)

# run PBM types and geographic region (with PBM pages control)
Types_GeoRegionPages_Loose <- brm(
  Total_Types ~ Geographic_Region + PBMPages_Num_z +
    (1 | gr(Node, cov = CX)),        # phylogenetic covariance
  data = SCCS_PBM_PA_Loose,
  family = poisson(),
  prior = priors,
  chains = 4,
  iter = 15000,
  init = 0,
  data2 = list(CX = covarMatrix),
  control = list(adapt_delta = 0.999, max_treedepth = 20)
)

saveRDS(Types_GeoRegionPages_Loose, "~/Types_GeoRegion_Pages_Loose.rds")

# run PBM subtypes and geographic region (with PBM pages control)
Subtypes_GeoRegionPages_Loose <- brm(
  Total_Subtypes ~ Geographic_Region + PBMPages_Num_z +
    (1 | gr(Node, cov = CX)),        # phylogenetic covariance
  data = SCCS_PBM_PA_Loose,
  family = poisson(),
  prior = priors,
  data2 = list(CX = covarMatrix),
  chains = 4,
  iter = 15000,
  init = 0,
  control = list(adapt_delta = 0.999, max_treedepth = 20)
)

saveRDS(Subtypes_GeoRegionPages_Loose, "~/Subtypes_GeoRegion_Pages_Loose.rds")

### Males and Females
## Types
# Males
Types_M_GeoRegionPages_Loose <- brm(
  Total_Types_M ~ Geographic_Region + PBMPages_Num_z +
    (1 | gr(Node, cov = CX)),        # phylogenetic covariance
  data = SCCS_PBM_PA_Loose,
  family = poisson(),
  prior = priors,
  chains = 4,
  iter = 15000,
  init = 0,
  data2 = list(CX = covarMatrix),
  control = list(adapt_delta = 0.999, max_treedepth = 20)
)

saveRDS(Types_M_GeoRegionPages_Loose, "~/Types_M_GeoRegion_Pages_Loose.rds")

# Females
Types_F_GeoRegionPages_Loose <- brm(
  Total_Types_F ~ Geographic_Region + PBMPages_Num_z +
    (1 | gr(Node, cov = CX)),        # phylogenetic covariance
  data = SCCS_PBM_PA_Loose,
  family = poisson(),
  prior = priors,
  chains = 4,
  iter = 15000,
  init = 0,
  data2 = list(CX = covarMatrix),
  control = list(adapt_delta = 0.999, max_treedepth = 20)
)

saveRDS(Types_F_GeoRegionPages_Loose, "~/Types_F_GeoRegion_Pages_Loose.rds")

## Subtypes
# Males
Subtypes_M_GeoRegionPages_Loose <- brm(
  Total_M_Subtypes ~ Geographic_Region + PBMPages_Num_z +
    (1 | gr(Node, cov = CX)),        # phylogenetic covariance
  data = SCCS_PBM_PA_Loose,
  family = poisson(),
  prior = priors,
  data2 = list(CX = covarMatrix),
  chains = 4,
  iter = 15000,
  init = 0,
  control = list(adapt_delta = 0.999, max_treedepth = 20)
)

saveRDS(Subtypes_M_GeoRegionPages_Loose, "~/Subtypes_M_GeoRegion_Pages_Loose.rds")

# Females
Subtypes_F_GeoRegionPages_Loose <- brm(
  Total_F_Subtypes ~ Geographic_Region + PBMPages_Num_z +
    (1 | gr(Node, cov = CX)),        # phylogenetic covariance
  data = SCCS_PBM_PA_Loose,
  family = poisson(),
  prior = priors,
  data2 = list(CX = covarMatrix),
  chains = 4,
  iter = 15000,
  init = 0,
  control = list(adapt_delta = 0.999, max_treedepth = 20)
)

saveRDS(Subtypes_F_GeoRegionPages_Loose, "~/Subtypes_F_GeoRegion_Pages_Loose.rds")
