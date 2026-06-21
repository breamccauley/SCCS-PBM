## Load the required packages
library(tidyverse)
library(dplyr)
library(ape)
library(phytools)

## TIGHT SAMPLE

#Read in NEXUS file of SCCS_Tree phylogenetic tree
SCCS_Tree<-read.tree("~/SCCS_Supertree.tre")

#Ladderize SCCS_Tree
SCCS_Tree<-ladderize(SCCS_Tree, right = TRUE)

# Read data into R
SCCS_PBM_PA_Tight <- read.csv("~/SCCS_PBM_Data_Tight.csv")

# Set row names as 'node' column
rownames(SCCS_PBM_PA_Tight) = SCCS_PBM_PA_Tight$Node

# Create subsets of data for correlated evolution (PBM Types)
TAT_MF_subset <- SCCS_PBM_PA_Tight %>%
  filter(TAT_M %in% c(0, 1),
         TAT_F %in% c(0, 1)) %>%
  select(Node, TAT_M, TAT_F)

SCA_MF_subset <- SCCS_PBM_PA_Tight %>%
  filter(SCA_M %in% c(0, 1),
         SCA_F %in% c(0, 1)) %>%
  select(Node, SCA_M, SCA_F)

PIE_MF_subset <- SCCS_PBM_PA_Tight %>%
  filter(PIE_M %in% c(0, 1),
         PIE_F %in% c(0, 1)) %>%
  select(Node, PIE_M, PIE_F)

GEN_MF_subset <- SCCS_PBM_PA_Tight %>%
  filter(GEN_M %in% c(0, 1),
         GEN_F %in% c(0, 1)) %>%
  select(Node, GEN_M, GEN_F)

AMP_MF_subset <- SCCS_PBM_PA_Tight %>%
  filter(AMP_M %in% c(0, 1),
         AMP_F %in% c(0, 1)) %>%
  select(Node, AMP_M, AMP_F)

DEN_MF_subset <- SCCS_PBM_PA_Tight %>%
  filter(DEN_M %in% c(0, 1),
         DEN_F %in% c(0, 1)) %>%
  select(Node, DEN_M, DEN_F)

SHA_MF_subset <- SCCS_PBM_PA_Tight %>%
  filter(SHA_M %in% c(0, 1),
         SHA_F %in% c(0, 1)) %>%
  select(Node, SHA_M, SHA_F)


## Tattooing
#Prune SCCS_Tree to  sample
TAT_MF_subset_Sample <- keep.tip(SCCS_Tree, TAT_MF_subset$Node)
#Avoid zero edge lengths
TAT_MF_subset_Sample$edge.length[TAT_MF_subset_Sample$edge.length==0] = 0.000001
#Order the data frame by the order of tree tiplabels
TAT_MF_subset <- TAT_MF_subset[TAT_MF_subset_Sample$tip.label,]

TAT_Male<-setNames(TAT_MF_subset$TAT_M,rownames(TAT_MF_subset))
TAT_Female<-setNames(TAT_MF_subset$TAT_F,rownames(TAT_MF_subset))

#Run pagel's discrete test of correlated evolution
fit.TAT_MF <- fitPagel(TAT_MF_subset_Sample,TAT_Male,TAT_Female)

saveRDS(fit.TAT_MF, "~/fit.TAT_MF_Tight.rds")

## Puncture tattooing
# Create subsets of data for correlated evolution
TAT_Punture_MF_subset <- SCCS_PBM_PA_Tight %>%
  filter(TAT_M_Puncture %in% c(0, 1),
         TAT_F_Puncture %in% c(0, 1)) %>%
  select(Node, TAT_M_Puncture, TAT_F_Puncture)

#Prune SCCS_Tree to  sample
TAT_Puncture_MF_subset_Sample <- keep.tip(SCCS_Tree, TAT_Punture_MF_subset$Node)
#Avoid zero edge lengths
TAT_Puncture_MF_subset_Sample$edge.length[TAT_Puncture_MF_subset_Sample$edge.length==0] = 0.000001
#Order the data frame by the order of tree tiplabels
TAT_Punture_MF_subset <- TAT_Punture_MF_subset[TAT_Puncture_MF_subset_Sample$tip.label,]

TAT_Puncture_Male<-setNames(TAT_Punture_MF_subset$TAT_M_Puncture,rownames(TAT_Punture_MF_subset))
TAT_Puncture_Female<-setNames(TAT_Punture_MF_subset$TAT_F_Puncture,rownames(TAT_Punture_MF_subset))

#Run pagel's discrete test of correlated evolution
fit.TAT_Puncture_MF <- fitPagel(TAT_Puncture_MF_subset_Sample,TAT_Puncture_Male,TAT_Puncture_Female)

saveRDS(fit.TAT_Puncture_MF, "~/fit.TAT_Puncture_MF_Tight.rds")

## Scarification
#Prune SCCS_Tree to  sample
SCA_MF_subset_Sample <- keep.tip(SCCS_Tree, SCA_MF_subset$Node)
#Avoid zero edge lengths
SCA_MF_subset_Sample$edge.length[SCA_MF_subset_Sample$edge.length==0] = 0.000001
#Order the data frame by the order of tree tiplabels
SCA_MF_subset <- SCA_MF_subset[SCA_MF_subset_Sample$tip.label,]

SCA_Male<-setNames(SCA_MF_subset$SCA_M,rownames(SCA_MF_subset))
SCA_Female<-setNames(SCA_MF_subset$SCA_F,rownames(SCA_MF_subset))

#Run pagel's discrete test of correlated evolution
fit.SCA_MF <- fitPagel(SCA_MF_subset_Sample,SCA_Male,SCA_Female)

saveRDS(fit.SCA_MF, "~/fit.SCA_MF_Tight.rds")


## Slicing scarification
# Create subsets of data for correlated evolution
SCA_Slice_MF_subset <- SCCS_PBM_PA_Tight %>%
  filter(SCA_M_Cut_Slice %in% c(0, 1),
         SCA_F_Cut_Slice %in% c(0, 1)) %>%
  select(Node, SCA_M_Cut_Slice, SCA_F_Cut_Slice)

#Prune SCCS_Tree to  sample
SCA_Slice_MF_subset_Sample <- keep.tip(SCCS_Tree, SCA_Slice_MF_subset$Node)
#Avoid zero edge lengths
SCA_Slice_MF_subset_Sample$edge.length[SCA_Slice_MF_subset_Sample$edge.length==0] = 0.000001
#Order the data frame by the order of tree tiplabels
SCA_Slice_MF_subset <- SCA_Slice_MF_subset[SCA_Slice_MF_subset_Sample$tip.label,]

SCA_Slice_Male<-setNames(SCA_Slice_MF_subset$SCA_M_Cut_Slice,rownames(SCA_Slice_MF_subset))
SCA_Slice_Female<-setNames(SCA_Slice_MF_subset$SCA_F_Cut_Slice,rownames(SCA_Slice_MF_subset))

#Run pagel's discrete test of correlated evolution
fit.SCA_Slice_MF <- fitPagel(SCA_Slice_MF_subset_Sample,SCA_Slice_Male,SCA_Slice_Female)

saveRDS(fit.SCA_Slice_MF, "~/fit.SCA_Slice_MF_Tight.rds")

## Piercing
#Prune SCCS_Tree to  sample
PIE_MF_subset_Sample <- keep.tip(SCCS_Tree, PIE_MF_subset$Node)
#Avoid zero edge lengths
PIE_MF_subset_Sample$edge.length[PIE_MF_subset_Sample$edge.length==0] = 0.000001
#Order the data frame by the order of tree tiplabels
PIE_MF_subset <- PIE_MF_subset[PIE_MF_subset_Sample$tip.label,]

PIE_Male<-setNames(PIE_MF_subset$PIE_M,rownames(PIE_MF_subset))
PIE_Female<-setNames(PIE_MF_subset$PIE_F,rownames(PIE_MF_subset))

#Run pagel's discrete test of correlated evolution
fit.PIE_MF <- fitPagel(PIE_MF_subset_Sample,PIE_Male,PIE_Female)

saveRDS(fit.PIE_MF, "~/fit.PIE_MF_Tight.rds")

## Enlarged earlobe piercing
# Create subsets of data for correlated evolution
PIE_LobeEnlarge_MF_subset <- SCCS_PBM_PA_Tight %>%
  filter(PIE_M_Earlobe_Enlarge %in% c(0, 1),
         PIE_F_Earlobe_Enlarge %in% c(0, 1)) %>%
  select(Node, PIE_M_Earlobe_Enlarge, PIE_F_Earlobe_Enlarge)

#Prune SCCS_Tree to  sample
PIE_LobeEnlarge_MF_subset_Sample <- keep.tip(SCCS_Tree, PIE_LobeEnlarge_MF_subset$Node)
#Avoid zero edge lengths
PIE_LobeEnlarge_MF_subset_Sample$edge.length[PIE_LobeEnlarge_MF_subset_Sample$edge.length==0] = 0.000001
#Order the data frame by the order of tree tiplabels
PIE_LobeEnlarge_MF_subset <- PIE_LobeEnlarge_MF_subset[PIE_LobeEnlarge_MF_subset_Sample$tip.label,]

PIE_LobeEnlarge_Male<-setNames(PIE_LobeEnlarge_MF_subset$PIE_M_Earlobe_Enlarge,rownames(PIE_LobeEnlarge_MF_subset))
PIE_LobeEnlarge_Female<-setNames(PIE_LobeEnlarge_MF_subset$PIE_F_Earlobe_Enlarge,rownames(PIE_LobeEnlarge_MF_subset))

#Run pagel's discrete test of correlated evolution
fit.PIE_LobeEnlarge_MF <- fitPagel(PIE_LobeEnlarge_MF_subset_Sample,PIE_LobeEnlarge_Male,PIE_LobeEnlarge_Female)

saveRDS(fit.PIE_LobeEnlarge_MF, "~/fit.PIE_LobeEnlarge_MF_Tight.rds")

## Septum piercing
# Create subsets of data for correlated evolution
PIE_Septum_MF_subset <- SCCS_PBM_PA_Tight %>%
  filter(PIE_M_Nose_Septum %in% c(0, 1),
         PIE_F_Nose_Septum %in% c(0, 1)) %>%
  select(Node, PIE_M_Nose_Septum, PIE_F_Nose_Septum)

#Prune SCCS_Tree to  sample
PIE_Septum_MF_subset_Sample <- keep.tip(SCCS_Tree, PIE_Septum_MF_subset$Node)
#Avoid zero edge lengths
PIE_Septum_MF_subset_Sample$edge.length[PIE_Septum_MF_subset_Sample$edge.length==0] = 0.000001
#Order the data frame by the order of tree tiplabels
PIE_Septum_MF_subset <- PIE_Septum_MF_subset[PIE_Septum_MF_subset_Sample$tip.label,]

PIE_Septum_Male<-setNames(PIE_Septum_MF_subset$PIE_M_Nose_Septum,rownames(PIE_Septum_MF_subset))
PIE_Septum_Female<-setNames(PIE_Septum_MF_subset$PIE_F_Nose_Septum,rownames(PIE_Septum_MF_subset))

#Run pagel's discrete test of correlated evolution
fit.PIE_Septum_MF <- fitPagel(PIE_Septum_MF_subset_Sample,PIE_Septum_Male,PIE_Septum_Female)

saveRDS(fit.PIE_Septum_MF, "~/fit.PIE_Septum_MF_Tight.rds")

## Labret piercing
# Create subsets of data for correlated evolution
PIE_Labret_MF_subset <- SCCS_PBM_PA_Tight %>%
  filter(PIE_M_Lip_Labret %in% c(0, 1),
         PIE_F_Lip_Labret %in% c(0, 1)) %>%
  select(Node, PIE_M_Lip_Labret, PIE_F_Lip_Labret)

#Prune SCCS_Tree to  sample
PIE_Labret_MF_subset_Sample <- keep.tip(SCCS_Tree, PIE_Labret_MF_subset$Node)
#Avoid zero edge lengths
PIE_Labret_MF_subset_Sample$edge.length[PIE_Labret_MF_subset_Sample$edge.length==0] = 0.000001
#Order the data frame by the order of tree tiplabels
PIE_Labret_MF_subset <- PIE_Labret_MF_subset[PIE_Labret_MF_subset_Sample$tip.label,]

PIE_Labret_Male<-setNames(PIE_Labret_MF_subset$PIE_M_Lip_Labret,rownames(PIE_Labret_MF_subset))
PIE_Labret_Female<-setNames(PIE_Labret_MF_subset$PIE_F_Lip_Labret,rownames(PIE_Labret_MF_subset))

#Run pagel's discrete test of correlated evolution
fit.PIE_Labret_MF <- fitPagel(PIE_Labret_MF_subset_Sample,PIE_Labret_Male,PIE_Labret_Female)

saveRDS(fit.PIE_Labret_MF, "~/fit.PIE_Labret_MF_Tight.rds")

## Genital Modification
#Prune SCCS_Tree to  sample
GEN_MF_subset_Sample <- keep.tip(SCCS_Tree, GEN_MF_subset$Node)
#Avoid zero edge lengths
GEN_MF_subset_Sample$edge.length[GEN_MF_subset_Sample$edge.length==0] = 0.000001
#Order the data frame by the order of tree tiplabels
GEN_MF_subset <- GEN_MF_subset[GEN_MF_subset_Sample$tip.label,]

GEN_Male<-setNames(GEN_MF_subset$GEN_M,rownames(GEN_MF_subset))
GEN_Female<-setNames(GEN_MF_subset$GEN_F,rownames(GEN_MF_subset))

#Run pagel's discrete test of correlated evolution
fit.GEN_MF <- fitPagel(GEN_MF_subset_Sample,GEN_Male,GEN_Female)

saveRDS(fit.GEN_MF, "~/fit.GEN_MF_Tight.rds")

## Circumcision and Clitoridectomy
# Create subsets of data for correlated evolution
GEN_CircumClitoridectomy_MF_subset <- SCCS_PBM_PA_Tight %>%
  filter(GEN_M_Circumcision %in% c(0, 1),
         GEN_F_Clitoridectomy %in% c(0, 1)) %>%
  select(Node, GEN_M_Circumcision, GEN_F_Clitoridectomy)

#Prune SCCS_Tree to  sample
GEN_CircumClitoridectomy_MF_subset_Sample <- keep.tip(SCCS_Tree, GEN_CircumClitoridectomy_MF_subset$Node)
#Avoid zero edge lengths
GEN_CircumClitoridectomy_MF_subset_Sample$edge.length[GEN_CircumClitoridectomy_MF_subset_Sample$edge.length==0] = 0.000001
#Order the data frame by the order of tree tiplabels
GEN_CircumClitoridectomy_MF_subset <- GEN_CircumClitoridectomy_MF_subset[GEN_CircumClitoridectomy_MF_subset_Sample$tip.label,]

GEN_CircumClitoridectomy_Male<-setNames(GEN_CircumClitoridectomy_MF_subset$GEN_M_Circumcision,rownames(GEN_CircumClitoridectomy_MF_subset))
GEN_CircumClitoridectomy_Female<-setNames(GEN_CircumClitoridectomy_MF_subset$GEN_F_Clitoridectomy,rownames(GEN_CircumClitoridectomy_MF_subset))

#Run pagel's discrete test of correlated evolution
fit.GEN_CircumClitoridectomy_MF <- fitPagel(GEN_CircumClitoridectomy_MF_subset_Sample,GEN_CircumClitoridectomy_Male,GEN_CircumClitoridectomy_Female)

saveRDS(fit.GEN_CircumClitoridectomy_MF, "~/fit.GEN_CircumClitoridectomy_MF_Tight.rds")

## Circumcision and Clitoral Excision
# Create subsets of data for correlated evolution
GEN_CircumClitExcis_MF_subset <- SCCS_PBM_PA_Tight %>%
  filter(GEN_M_Circumcision %in% c(0, 1),
         GEN_F_ClitoralExcision %in% c(0, 1)) %>%
  select(Node, GEN_M_Circumcision, GEN_F_ClitoralExcision)

#Prune SCCS_Tree to  sample
GEN_CircumClitExcis_MF_subset_Sample <- keep.tip(SCCS_Tree, GEN_CircumClitExcis_MF_subset$Node)
#Avoid zero edge lengths
GEN_CircumClitExcis_MF_subset_Sample$edge.length[GEN_CircumClitExcis_MF_subset_Sample$edge.length==0] = 0.000001
#Order the data frame by the order of tree tiplabels
GEN_CircumClitExcis_MF_subset <- GEN_CircumClitExcis_MF_subset[GEN_CircumClitExcis_MF_subset_Sample$tip.label,]

GEN_CircumClitExcis_Male<-setNames(GEN_CircumClitExcis_MF_subset$GEN_M_Circumcision,rownames(GEN_CircumClitExcis_MF_subset))
GEN_CircumClitExcis_Female<-setNames(GEN_CircumClitExcis_MF_subset$GEN_F_ClitoralExcision,rownames(GEN_CircumClitExcis_MF_subset))

#Run pagel's discrete test of correlated evolution
fit.GEN_CircumClitExcis_MF <- fitPagel(GEN_CircumClitExcis_MF_subset_Sample,GEN_CircumClitExcis_Male,GEN_CircumClitExcis_Female)

saveRDS(fit.GEN_CircumClitExcis_MF, "~/fit.GEN_CircumClitExcis_MF_Tight.rds")

## Circumcision and Clitoral Cutting
# Create new column that merges clitoridectomy and clitoral excision into one
SCCS_PBM_PA_Tight <- SCCS_PBM_PA_Tight %>%
  mutate(
    GEN_F_ClitorisCut = case_when(
      GEN_F_Clitoridectomy == 1 | GEN_F_ClitoralExcision == 1 ~ 1,
      GEN_F_Clitoridectomy == 0 | GEN_F_ClitoralExcision == 0 ~ 0,
      TRUE ~ NA_real_
    )
  )

# Create subsets of data for correlated evolution
GEN_CircumClitorisCut_MF_subset <- SCCS_PBM_PA_Tight %>%
  filter(GEN_M_Circumcision %in% c(0, 1),
         GEN_F_ClitorisCut %in% c(0, 1)) %>%
  select(Node, GEN_M_Circumcision, GEN_F_ClitorisCut)

#Prune SCCS_Tree to  sample
GEN_CircumClitorisCut_MF_subset_Sample <- keep.tip(SCCS_Tree, GEN_CircumClitorisCut_MF_subset$Node)
#Avoid zero edge lengths
GEN_CircumClitorisCut_MF_subset_Sample$edge.length[GEN_CircumClitorisCut_MF_subset_Sample$edge.length==0] = 0.000001
#Order the data frame by the order of tree tiplabels
GEN_CircumClitorisCut_MF_subset <- GEN_CircumClitorisCut_MF_subset[GEN_CircumClitorisCut_MF_subset_Sample$tip.label,]

GEN_CircumClitorisCut_Male<-setNames(GEN_CircumClitorisCut_MF_subset$GEN_M_Circumcision,rownames(GEN_CircumClitorisCut_MF_subset))
GEN_CircumClitorisCut_Female<-setNames(GEN_CircumClitorisCut_MF_subset$GEN_F_ClitorisCut,rownames(GEN_CircumClitorisCut_MF_subset))

#Run pagel's discrete test of correlated evolution
fit.GEN_CircumClitorisCut_MF <- fitPagel(GEN_CircumClitorisCut_MF_subset_Sample,GEN_CircumClitorisCut_Male,GEN_CircumClitorisCut_Female)

saveRDS(fit.GEN_CircumClitorisCut_MF, "~/fit.GEN_CircumClitorisCut_MF_Tight.rds")

## Circumcision and Vulva Cutting
# Create new column that merges clitoridectomy, clitoral excision, labial excision, introcision, and infibulation into one
vars <- c(
  "GEN_F_Clitoridectomy",
  "GEN_F_ClitoralExcision",
  "GEN_F_LabialExcision",
  "GEN_F_Infibulation",
  "GEN_F_Introcision"
)

SCCS_PBM_PA_Tight <- SCCS_PBM_PA_Tight %>%
  mutate(
    GEN_F_AllCutting = case_when(
      if_any(all_of(vars), ~ . == 1) ~ 1,
      if_any(all_of(vars), ~ . == 0) ~ 0,
      TRUE ~ NA_real_
    )
  )

# Create subsets of data for correlated evolution
GEN_CircumAllCutting_MF_subset <- SCCS_PBM_PA_Tight %>%
  filter(GEN_M_Circumcision %in% c(0, 1),
         GEN_F_AllCutting %in% c(0, 1)) %>%
  select(Node, GEN_M_Circumcision, GEN_F_AllCutting)

#Prune SCCS_Tree to  sample
GEN_CircumAllCutting_MF_subset_Sample <- keep.tip(SCCS_Tree, GEN_CircumAllCutting_MF_subset$Node)
#Avoid zero edge lengths
GEN_CircumAllCutting_MF_subset_Sample$edge.length[GEN_CircumAllCutting_MF_subset_Sample$edge.length==0] = 0.000001
#Order the data frame by the order of tree tiplabels
GEN_CircumAllCutting_MF_subset <- GEN_CircumAllCutting_MF_subset[GEN_CircumAllCutting_MF_subset_Sample$tip.label,]

GEN_CircumAllCutting_Male<-setNames(GEN_CircumAllCutting_MF_subset$GEN_M_Circumcision,rownames(GEN_CircumAllCutting_MF_subset))
GEN_CircumAllCutting_Female<-setNames(GEN_CircumAllCutting_MF_subset$GEN_F_AllCutting,rownames(GEN_CircumAllCutting_MF_subset))

#Run pagel's discrete test of correlated evolution
fit.GEN_CircumAllCutting_MF <- fitPagel(GEN_CircumAllCutting_MF_subset_Sample,GEN_CircumAllCutting_Male,GEN_CircumAllCutting_Female)

saveRDS(fit.GEN_CircumAllCutting_MF, "~/fit.GEN_CircumAllCutting_MF_Tight.rds")

## Circumcision and Labial Excision
# Create subsets of data for correlated evolution
GEN_CircumLabExcise_MF_subset <- SCCS_PBM_PA_Tight %>%
  filter(GEN_M_Circumcision %in% c(0, 1),
         GEN_F_LabialExcision %in% c(0, 1)) %>%
  select(Node, GEN_M_Circumcision, GEN_F_LabialExcision)

#Prune SCCS_Tree to  sample
GEN_CircumLabExcise_MF_subset_Sample <- keep.tip(SCCS_Tree, GEN_CircumLabExcise_MF_subset$Node)
#Avoid zero edge lengths
GEN_CircumLabExcise_MF_subset_Sample$edge.length[GEN_CircumLabExcise_MF_subset_Sample$edge.length==0] = 0.000001
#Order the data frame by the order of tree tiplabels
GEN_CircumLabExcise_MF_subset <- GEN_CircumLabExcise_MF_subset[GEN_CircumLabExcise_MF_subset_Sample$tip.label,]

GEN_CircumLabExcise_Male<-setNames(GEN_CircumLabExcise_MF_subset$GEN_M_Circumcision,rownames(GEN_CircumLabExcise_MF_subset))
GEN_CircumLabExcise_Female<-setNames(GEN_CircumLabExcise_MF_subset$GEN_F_LabialExcision,rownames(GEN_CircumLabExcise_MF_subset))

#Run pagel's discrete test of correlated evolution
fit.GEN_CircumLabExcise_MF <- fitPagel(GEN_CircumLabExcise_MF_subset_Sample,GEN_CircumLabExcise_Male,GEN_CircumLabExcise_Female)

saveRDS(fit.GEN_CircumLabExcise_MF, "~/fit.GEN_CircumLabExcise_MF_Tight.rds")

## Circumcision and Vulva Stretching
# Create subsets of data for correlated evolution
GEN_CircumStretch_MF_subset <- SCCS_PBM_PA_Tight %>%
  filter(GEN_M_Circumcision %in% c(0, 1),
         GEN_F_Stretch %in% c(0, 1)) %>%
  select(Node, GEN_M_Circumcision, GEN_F_Stretch)

#Prune SCCS_Tree to  sample
GEN_CircumStretch_MF_subset_Sample <- keep.tip(SCCS_Tree, GEN_CircumStretch_MF_subset$Node)
#Avoid zero edge lengths
GEN_CircumStretch_MF_subset_Sample$edge.length[GEN_CircumStretch_MF_subset_Sample$edge.length==0] = 0.000001
#Order the data frame by the order of tree tiplabels
GEN_CircumStretch_MF_subset <- GEN_CircumStretch_MF_subset[GEN_CircumStretch_MF_subset_Sample$tip.label,]

GEN_CircumStretch_Male<-setNames(GEN_CircumStretch_MF_subset$GEN_M_Circumcision,rownames(GEN_CircumStretch_MF_subset))
GEN_CircumStretch_Female<-setNames(GEN_CircumStretch_MF_subset$GEN_F_Stretch,rownames(GEN_CircumStretch_MF_subset))

#Run pagel's discrete test of correlated evolution
fit.GEN_CircumStretch_MF <- fitPagel(GEN_CircumStretch_MF_subset_Sample,GEN_CircumStretch_Male,GEN_CircumStretch_Female)

saveRDS(fit.GEN_CircumStretch_MF, "~/fit.GEN_CircumStretch_MF_Tight.rds")

## Amputation
#Prune SCCS_Tree to  sample
AMP_MF_subset_Sample <- keep.tip(SCCS_Tree, AMP_MF_subset$Node)
#Avoid zero edge lengths
AMP_MF_subset_Sample$edge.length[AMP_MF_subset_Sample$edge.length==0] = 0.000001
#Order the data frame by the order of tree tiplabels
AMP_MF_subset <- AMP_MF_subset[AMP_MF_subset_Sample$tip.label,]

AMP_Male<-setNames(AMP_MF_subset$AMP_M,rownames(AMP_MF_subset))
AMP_Female<-setNames(AMP_MF_subset$AMP_F,rownames(AMP_MF_subset))

#Run pagel's discrete test of correlated evolution
fit.AMP_MF <- fitPagel(AMP_MF_subset_Sample,AMP_Male,AMP_Female)

saveRDS(fit.AMP_MF, "~/fit.AMP_MF_Tight.rds")

## Dental Modification
#Prune SCCS_Tree to  sample
DEN_MF_subset_Sample <- keep.tip(SCCS_Tree, DEN_MF_subset$Node)
#Avoid zero edge lengths
DEN_MF_subset_Sample$edge.length[DEN_MF_subset_Sample$edge.length==0] = 0.000001
#Order the data frame by the order of tree tiplabels
DEN_MF_subset <- DEN_MF_subset[DEN_MF_subset_Sample$tip.label,]

DEN_Male<-setNames(DEN_MF_subset$DEN_M,rownames(DEN_MF_subset))
DEN_Female<-setNames(DEN_MF_subset$DEN_F,rownames(DEN_MF_subset))

#Run pagel's discrete test of correlated evolution
fit.DEN_MF <- fitPagel(DEN_MF_subset_Sample,DEN_Male,DEN_Female)

saveRDS(fit.DEN_MF, "~/fit.DEN_MF_Tight.rds")

## Bone Shaping
#Prune SCCS_Tree to  sample
SHA_MF_subset_Sample <- keep.tip(SCCS_Tree, SHA_MF_subset$Node)
#Avoid zero edge lengths
SHA_MF_subset_Sample$edge.length[SHA_MF_subset_Sample$edge.length==0] = 0.000001
#Order the data frame by the order of tree tiplabels
SHA_MF_subset <- SHA_MF_subset[SHA_MF_subset_Sample$tip.label,]

SHA_Male<-setNames(SHA_MF_subset$SHA_M,rownames(SHA_MF_subset))
SHA_Female<-setNames(SHA_MF_subset$SHA_F,rownames(SHA_MF_subset))

#Run pagel's discrete test of correlated evolution
fit.SHA_MF <- fitPagel(SHA_MF_subset_Sample,SHA_Male,SHA_Female)

saveRDS(fit.SHA_MF, "~/fit.SHA_MF_Tight.rds")

## Cranial modification (actually exact same as general bone shaping)
# Create subsets of data for correlated evolution
SHA_Cranial_MF_subset <- SCCS_PBM_PA_Tight %>%
  filter(SHA_M_Cranial %in% c(0, 1),
         SHA_F_Cranial %in% c(0, 1)) %>%
  select(Node, SHA_M_Cranial, SHA_F_Cranial)

#Prune SCCS_Tree to  sample
SHA_Cranial_MF_subset_Sample <- keep.tip(SCCS_Tree, SHA_Cranial_MF_subset$Node)
#Avoid zero edge lengths
SHA_Cranial_MF_subset_Sample$edge.length[SHA_Cranial_MF_subset_Sample$edge.length==0] = 0.000001
#Order the data frame by the order of tree tiplabels
SHA_Cranial_MF_subset <- SHA_Cranial_MF_subset[SHA_Cranial_MF_subset_Sample$tip.label,]

SHA_Cranial_Male<-setNames(SHA_Cranial_MF_subset$SHA_M_Cranial,rownames(SHA_Cranial_MF_subset))
SHA_Cranial_Female<-setNames(SHA_Cranial_MF_subset$SHA_F_Cranial,rownames(SHA_Cranial_MF_subset))

#Run pagel's discrete test of correlated evolution
fit.SHA_Cranial_MF <- fitPagel(SHA_Cranial_MF_subset_Sample,SHA_Cranial_Male,SHA_Cranial_Female)

saveRDS(fit.SHA_Cranial_MF, "~/fit.SHA_Cranial_MF_Tight.rds")

# LOOSE SAMPLE
# Read data into R
SCCS_PBM_PA_Loose <- read.csv("~/SCCS_PBM_Data_Loose.csv")

# Set row names as 'node' column
rownames(SCCS_PBM_PA_Loose) = SCCS_PBM_PA_Loose$Node

# Create subsets of data for correlated evolution (PBM Types)
TAT_MF_subset <- SCCS_PBM_PA_Loose %>%
  filter(TAT_M %in% c(0, 1),
         TAT_F %in% c(0, 1)) %>%
  select(Node, TAT_M, TAT_F)

SCA_MF_subset <- SCCS_PBM_PA_Loose %>%
  filter(SCA_M %in% c(0, 1),
         SCA_F %in% c(0, 1)) %>%
  select(Node, SCA_M, SCA_F)

PIE_MF_subset <- SCCS_PBM_PA_Loose %>%
  filter(PIE_M %in% c(0, 1),
         PIE_F %in% c(0, 1)) %>%
  select(Node, PIE_M, PIE_F)

GEN_MF_subset <- SCCS_PBM_PA_Loose %>%
  filter(GEN_M %in% c(0, 1),
         GEN_F %in% c(0, 1)) %>%
  select(Node, GEN_M, GEN_F)

AMP_MF_subset <- SCCS_PBM_PA_Loose %>%
  filter(AMP_M %in% c(0, 1),
         AMP_F %in% c(0, 1)) %>%
  select(Node, AMP_M, AMP_F)

DEN_MF_subset <- SCCS_PBM_PA_Loose %>%
  filter(DEN_M %in% c(0, 1),
         DEN_F %in% c(0, 1)) %>%
  select(Node, DEN_M, DEN_F)

SHA_MF_subset <- SCCS_PBM_PA_Loose %>%
  filter(SHA_M %in% c(0, 1),
         SHA_F %in% c(0, 1)) %>%
  select(Node, SHA_M, SHA_F)

## Tattooing
#Prune SCCS_Tree to  sample
TAT_MF_subset_Sample <- keep.tip(SCCS_Tree, TAT_MF_subset$Node)
#Avoid zero edge lengths
TAT_MF_subset_Sample$edge.length[TAT_MF_subset_Sample$edge.length==0] = 0.000001
#Order the data frame by the order of tree tiplabels
TAT_MF_subset <- TAT_MF_subset[TAT_MF_subset_Sample$tip.label,]

TAT_Male<-setNames(TAT_MF_subset$TAT_M,rownames(TAT_MF_subset))
TAT_Female<-setNames(TAT_MF_subset$TAT_F,rownames(TAT_MF_subset))

#Run pagel's discrete test of correlated evolution
fit.TAT_MF_Loose <- fitPagel(TAT_MF_subset_Sample,TAT_Male,TAT_Female)

saveRDS(fit.TAT_MF_Loose, "~/fit.TAT_MF_Loose.rds")

# Create subsets of data for correlated evolution
TAT_Punture_MF_subset <- SCCS_PBM_PA_Loose %>%
  filter(TAT_M_Puncture %in% c(0, 1),
         TAT_F_Puncture %in% c(0, 1)) %>%
  select(Node, TAT_M_Puncture, TAT_F_Puncture)

#Prune SCCS_Tree to  sample
TAT_Puncture_MF_subset_Sample <- keep.tip(SCCS_Tree, TAT_Punture_MF_subset$Node)
#Avoid zero edge lengths
TAT_Puncture_MF_subset_Sample$edge.length[TAT_Puncture_MF_subset_Sample$edge.length==0] = 0.000001
#Order the data frame by the order of tree tiplabels
TAT_Punture_MF_subset <- TAT_Punture_MF_subset[TAT_Puncture_MF_subset_Sample$tip.label,]

TAT_Puncture_Male<-setNames(TAT_Punture_MF_subset$TAT_M_Puncture,rownames(TAT_Punture_MF_subset))
TAT_Puncture_Female<-setNames(TAT_Punture_MF_subset$TAT_F_Puncture,rownames(TAT_Punture_MF_subset))

#Run pagel's discrete test of correlated evolution
fit.TAT_Puncture_MF <- fitPagel(TAT_Puncture_MF_subset_Sample,TAT_Puncture_Male,TAT_Puncture_Female)

saveRDS(fit.TAT_Puncture_MF, "~/fit.TAT_Puncture_MF_Loose.rds")

## Scarification
#Prune SCCS_Tree to  sample
SCA_MF_subset_Sample <- keep.tip(SCCS_Tree, SCA_MF_subset$Node)
#Avoid zero edge lengths
SCA_MF_subset_Sample$edge.length[SCA_MF_subset_Sample$edge.length==0] = 0.000001
#Order the data frame by the order of tree tiplabels
SCA_MF_subset <- SCA_MF_subset[SCA_MF_subset_Sample$tip.label,]

SCA_Male<-setNames(SCA_MF_subset$SCA_M,rownames(SCA_MF_subset))
SCA_Female<-setNames(SCA_MF_subset$SCA_F,rownames(SCA_MF_subset))

#Run pagel's discrete test of correlated evolution
fit.SCA_MF_Loose <- fitPagel(SCA_MF_subset_Sample,SCA_Male,SCA_Female)

saveRDS(fit.SCA_MF_Loose, "~/fit.SCA_MF_Loose.rds")

## Slicing scarification
# Create subsets of data for correlated evolution
SCA_Slice_MF_subset <- SCCS_PBM_PA_Loose %>%
  filter(SCA_M_Cut_Slice %in% c(0, 1),
         SCA_F_Cut_Slice %in% c(0, 1)) %>%
  select(Node, SCA_M_Cut_Slice, SCA_F_Cut_Slice)

#Prune SCCS_Tree to  sample
SCA_Slice_MF_subset_Sample <- keep.tip(SCCS_Tree, SCA_Slice_MF_subset$Node)
#Avoid zero edge lengths
SCA_Slice_MF_subset_Sample$edge.length[SCA_Slice_MF_subset_Sample$edge.length==0] = 0.000001
#Order the data frame by the order of tree tiplabels
SCA_Slice_MF_subset <- SCA_Slice_MF_subset[SCA_Slice_MF_subset_Sample$tip.label,]

SCA_Slice_Male<-setNames(SCA_Slice_MF_subset$SCA_M_Cut_Slice,rownames(SCA_Slice_MF_subset))
SCA_Slice_Female<-setNames(SCA_Slice_MF_subset$SCA_F_Cut_Slice,rownames(SCA_Slice_MF_subset))

#Run pagel's discrete test of correlated evolution
fit.SCA_Slice_MF <- fitPagel(SCA_Slice_MF_subset_Sample,SCA_Slice_Male,SCA_Slice_Female)

saveRDS(fit.SCA_Slice_MF, "~/fit.SCA_Slice_MF_Loose.rds")


## Piercing
#Prune SCCS_Tree to  sample
PIE_MF_subset_Sample <- keep.tip(SCCS_Tree, PIE_MF_subset$Node)
#Avoid zero edge lengths
PIE_MF_subset_Sample$edge.length[PIE_MF_subset_Sample$edge.length==0] = 0.000001
#Order the data frame by the order of tree tiplabels
PIE_MF_subset <- PIE_MF_subset[PIE_MF_subset_Sample$tip.label,]

PIE_Male<-setNames(PIE_MF_subset$PIE_M,rownames(PIE_MF_subset))
PIE_Female<-setNames(PIE_MF_subset$PIE_F,rownames(PIE_MF_subset))

#Run pagel's discrete test of correlated evolution
fit.PIE_MF_Loose <- fitPagel(PIE_MF_subset_Sample,PIE_Male,PIE_Female)

saveRDS(fit.PIE_MF_Loose, "~/fit.PIE_MF_Loose.rds")

## Enlarged Earlobe Piercing
# Create subsets of data for correlated evolution
PIE_LobeEnlarge_MF_subset <- SCCS_PBM_PA_Loose %>%
  filter(PIE_M_Earlobe_Enlarge %in% c(0, 1),
         PIE_F_Earlobe_Enlarge %in% c(0, 1)) %>%
  select(Node, PIE_M_Earlobe_Enlarge, PIE_F_Earlobe_Enlarge)

#Prune SCCS_Tree to  sample
PIE_LobeEnlarge_MF_subset_Sample <- keep.tip(SCCS_Tree, PIE_LobeEnlarge_MF_subset$Node)
#Avoid zero edge lengths
PIE_LobeEnlarge_MF_subset_Sample$edge.length[PIE_LobeEnlarge_MF_subset_Sample$edge.length==0] = 0.000001
#Order the data frame by the order of tree tiplabels
PIE_LobeEnlarge_MF_subset <- PIE_LobeEnlarge_MF_subset[PIE_LobeEnlarge_MF_subset_Sample$tip.label,]

PIE_LobeEnlarge_Male<-setNames(PIE_LobeEnlarge_MF_subset$PIE_M_Earlobe_Enlarge,rownames(PIE_LobeEnlarge_MF_subset))
PIE_LobeEnlarge_Female<-setNames(PIE_LobeEnlarge_MF_subset$PIE_F_Earlobe_Enlarge,rownames(PIE_LobeEnlarge_MF_subset))

#Run pagel's discrete test of correlated evolution
fit.PIE_LobeEnlarge_MF <- fitPagel(PIE_LobeEnlarge_MF_subset_Sample,PIE_LobeEnlarge_Male,PIE_LobeEnlarge_Female)

saveRDS(fit.PIE_LobeEnlarge_MF, "~/fit.PIE_LobeEnlarge_MF_Loose.rds")

## Septum Piercing
# Create subsets of data for correlated evolution
PIE_Septum_MF_subset <- SCCS_PBM_PA_Loose %>%
  filter(PIE_M_Nose_Septum %in% c(0, 1),
         PIE_F_Nose_Septum %in% c(0, 1)) %>%
  select(Node, PIE_M_Nose_Septum, PIE_F_Nose_Septum)

#Prune SCCS_Tree to  sample
PIE_Septum_MF_subset_Sample <- keep.tip(SCCS_Tree, PIE_Septum_MF_subset$Node)
#Avoid zero edge lengths
PIE_Septum_MF_subset_Sample$edge.length[PIE_Septum_MF_subset_Sample$edge.length==0] = 0.000001
#Order the data frame by the order of tree tiplabels
PIE_Septum_MF_subset <- PIE_Septum_MF_subset[PIE_Septum_MF_subset_Sample$tip.label,]

PIE_Septum_Male<-setNames(PIE_Septum_MF_subset$PIE_M_Nose_Septum,rownames(PIE_Septum_MF_subset))
PIE_Septum_Female<-setNames(PIE_Septum_MF_subset$PIE_F_Nose_Septum,rownames(PIE_Septum_MF_subset))

#Run pagel's discrete test of correlated evolution
fit.PIE_Septum_MF <- fitPagel(PIE_Septum_MF_subset_Sample,PIE_Septum_Male,PIE_Septum_Female)

saveRDS(fit.PIE_Septum_MF, "~/fit.PIE_Septum_MF_Loose.rds")

## Labret Piercing
# Create subsets of data for correlated evolution
PIE_Labret_MF_subset <- SCCS_PBM_PA_Loose %>%
  filter(PIE_M_Lip_Labret %in% c(0, 1),
         PIE_F_Lip_Labret %in% c(0, 1)) %>%
  select(Node, PIE_M_Lip_Labret, PIE_F_Lip_Labret)

#Prune SCCS_Tree to  sample
PIE_Labret_MF_subset_Sample <- keep.tip(SCCS_Tree, PIE_Labret_MF_subset$Node)
#Avoid zero edge lengths
PIE_Labret_MF_subset_Sample$edge.length[PIE_Labret_MF_subset_Sample$edge.length==0] = 0.000001
#Order the data frame by the order of tree tiplabels
PIE_Labret_MF_subset <- PIE_Labret_MF_subset[PIE_Labret_MF_subset_Sample$tip.label,]

PIE_Labret_Male<-setNames(PIE_Labret_MF_subset$PIE_M_Lip_Labret,rownames(PIE_Labret_MF_subset))
PIE_Labret_Female<-setNames(PIE_Labret_MF_subset$PIE_F_Lip_Labret,rownames(PIE_Labret_MF_subset))

#Run pagel's discrete test of correlated evolution
fit.PIE_Labret_MF <- fitPagel(PIE_Labret_MF_subset_Sample,PIE_Labret_Male,PIE_Labret_Female)

saveRDS(fit.PIE_Labret_MF, "~/fit.PIE_Labret_MF_Loose.rds")

## Genital Modification
#Prune SCCS_Tree to  sample
GEN_MF_subset_Sample <- keep.tip(SCCS_Tree, GEN_MF_subset$Node)
#Avoid zero edge lengths
GEN_MF_subset_Sample$edge.length[GEN_MF_subset_Sample$edge.length==0] = 0.000001
#Order the data frame by the order of tree tiplabels
GEN_MF_subset <- GEN_MF_subset[GEN_MF_subset_Sample$tip.label,]

GEN_Male<-setNames(GEN_MF_subset$GEN_M,rownames(GEN_MF_subset))
GEN_Female<-setNames(GEN_MF_subset$GEN_F,rownames(GEN_MF_subset))

#Run pagel's discrete test of correlated evolution
fit.GEN_MF_Loose <- fitPagel(GEN_MF_subset_Sample,GEN_Male,GEN_Female)

saveRDS(fit.GEN_MF_Loose, "~/fit.GEN_MF_Loose.rds")

## Circumcision and Clitoridectomy
# Create subsets of data for correlated evolution
GEN_CircumClitoridectomy_MF_subset <- SCCS_PBM_PA_Loose %>%
  filter(GEN_M_Circumcision %in% c(0, 1),
         GEN_F_Clitoridectomy %in% c(0, 1)) %>%
  select(Node, GEN_M_Circumcision, GEN_F_Clitoridectomy)

#Prune SCCS_Tree to  sample
GEN_CircumClitoridectomy_MF_subset_Sample <- keep.tip(SCCS_Tree, GEN_CircumClitoridectomy_MF_subset$Node)
#Avoid zero edge lengths
GEN_CircumClitoridectomy_MF_subset_Sample$edge.length[GEN_CircumClitoridectomy_MF_subset_Sample$edge.length==0] = 0.000001
#Order the data frame by the order of tree tiplabels
GEN_CircumClitoridectomy_MF_subset <- GEN_CircumClitoridectomy_MF_subset[GEN_CircumClitoridectomy_MF_subset_Sample$tip.label,]

GEN_CircumClitoridectomy_Male<-setNames(GEN_CircumClitoridectomy_MF_subset$GEN_M_Circumcision,rownames(GEN_CircumClitoridectomy_MF_subset))
GEN_CircumClitoridectomy_Female<-setNames(GEN_CircumClitoridectomy_MF_subset$GEN_F_Clitoridectomy,rownames(GEN_CircumClitoridectomy_MF_subset))

#Run pagel's discrete test of correlated evolution
fit.GEN_CircumClitoridectomy_MF <- fitPagel(GEN_CircumClitoridectomy_MF_subset_Sample,GEN_CircumClitoridectomy_Male,GEN_CircumClitoridectomy_Female)

saveRDS(fit.GEN_CircumClitoridectomy_MF, "~/fit.GEN_CircumClitoridectomy_MF_Loose.rds")

## Circumcision and Clitoral Excision
# Create subsets of data for correlated evolution
GEN_CircumClitExcis_MF_subset <- SCCS_PBM_PA_Loose %>%
  filter(GEN_M_Circumcision %in% c(0, 1),
         GEN_F_ClitoralExcision %in% c(0, 1)) %>%
  select(Node, GEN_M_Circumcision, GEN_F_ClitoralExcision)

#Prune SCCS_Tree to  sample
GEN_CircumClitExcis_MF_subset_Sample <- keep.tip(SCCS_Tree, GEN_CircumClitExcis_MF_subset$Node)
#Avoid zero edge lengths
GEN_CircumClitExcis_MF_subset_Sample$edge.length[GEN_CircumClitExcis_MF_subset_Sample$edge.length==0] = 0.000001
#Order the data frame by the order of tree tiplabels
GEN_CircumClitExcis_MF_subset <- GEN_CircumClitExcis_MF_subset[GEN_CircumClitExcis_MF_subset_Sample$tip.label,]

GEN_CircumClitExcis_Male<-setNames(GEN_CircumClitExcis_MF_subset$GEN_M_Circumcision,rownames(GEN_CircumClitExcis_MF_subset))
GEN_CircumClitExcis_Female<-setNames(GEN_CircumClitExcis_MF_subset$GEN_F_ClitoralExcision,rownames(GEN_CircumClitExcis_MF_subset))

#Run pagel's discrete test of correlated evolution
fit.GEN_CircumClitExcis_MF <- fitPagel(GEN_CircumClitExcis_MF_subset_Sample,GEN_CircumClitExcis_Male,GEN_CircumClitExcis_Female)

saveRDS(fit.GEN_CircumClitExcis_MF, "~/fit.GEN_CircumClitExcis_MF_Loose.rds")

## Circumcision and Clitoral Cutting
# Create new column that merges clitoridectomy and clitoral excision into one
SCCS_PBM_PA_Loose <- SCCS_PBM_PA_Loose %>%
  mutate(
    GEN_F_ClitorisCut = case_when(
      GEN_F_Clitoridectomy == 1 | GEN_F_ClitoralExcision == 1 ~ 1,
      GEN_F_Clitoridectomy == 0 | GEN_F_ClitoralExcision == 0 ~ 0,
      TRUE ~ NA_real_
    )
  )

# Create subsets of data for correlated evolution
GEN_CircumClitorisCut_MF_subset <- SCCS_PBM_PA_Loose %>%
  filter(GEN_M_Circumcision %in% c(0, 1),
         GEN_F_ClitorisCut %in% c(0, 1)) %>%
  select(Node, GEN_M_Circumcision, GEN_F_ClitorisCut)

#Prune SCCS_Tree to  sample
GEN_CircumClitorisCut_MF_subset_Sample <- keep.tip(SCCS_Tree, GEN_CircumClitorisCut_MF_subset$Node)
#Avoid zero edge lengths
GEN_CircumClitorisCut_MF_subset_Sample$edge.length[GEN_CircumClitorisCut_MF_subset_Sample$edge.length==0] = 0.000001
#Order the data frame by the order of tree tiplabels
GEN_CircumClitorisCut_MF_subset <- GEN_CircumClitorisCut_MF_subset[GEN_CircumClitorisCut_MF_subset_Sample$tip.label,]

GEN_CircumClitorisCut_Male<-setNames(GEN_CircumClitorisCut_MF_subset$GEN_M_Circumcision,rownames(GEN_CircumClitorisCut_MF_subset))
GEN_CircumClitorisCut_Female<-setNames(GEN_CircumClitorisCut_MF_subset$GEN_F_ClitorisCut,rownames(GEN_CircumClitorisCut_MF_subset))

#Run pagel's discrete test of correlated evolution
fit.GEN_CircumClitorisCut_MF <- fitPagel(GEN_CircumClitorisCut_MF_subset_Sample,GEN_CircumClitorisCut_Male,GEN_CircumClitorisCut_Female)

saveRDS(fit.GEN_CircumClitorisCut_MF, "~/fit.GEN_CircumClitorisCut_MF_Loose.rds")

## Circumcision and Vulva Cutting
# Create new column that merges clitoridectomy, clitoral excision, labial excision, introcision, and infibulation into one
vars <- c(
  "GEN_F_Clitoridectomy",
  "GEN_F_ClitoralExcision",
  "GEN_F_LabialExcision",
  "GEN_F_Infibulation",
  "GEN_F_Introcision"
)

SCCS_PBM_PA_Loose <- SCCS_PBM_PA_Loose %>%
  mutate(
    GEN_F_AllCutting = case_when(
      if_any(all_of(vars), ~ . == 1) ~ 1,
      if_any(all_of(vars), ~ . == 0) ~ 0,
      TRUE ~ NA_real_
    )
  )

# Create subsets of data for correlated evolution
GEN_CircumAllCutting_MF_subset <- SCCS_PBM_PA_Loose %>%
  filter(GEN_M_Circumcision %in% c(0, 1),
         GEN_F_AllCutting %in% c(0, 1)) %>%
  select(Node, GEN_M_Circumcision, GEN_F_AllCutting)

#Prune SCCS_Tree to  sample
GEN_CircumAllCutting_MF_subset_Sample <- keep.tip(SCCS_Tree, GEN_CircumAllCutting_MF_subset$Node)
#Avoid zero edge lengths
GEN_CircumAllCutting_MF_subset_Sample$edge.length[GEN_CircumAllCutting_MF_subset_Sample$edge.length==0] = 0.000001
#Order the data frame by the order of tree tiplabels
GEN_CircumAllCutting_MF_subset <- GEN_CircumAllCutting_MF_subset[GEN_CircumAllCutting_MF_subset_Sample$tip.label,]

GEN_CircumAllCutting_Male<-setNames(GEN_CircumAllCutting_MF_subset$GEN_M_Circumcision,rownames(GEN_CircumAllCutting_MF_subset))
GEN_CircumAllCutting_Female<-setNames(GEN_CircumAllCutting_MF_subset$GEN_F_AllCutting,rownames(GEN_CircumAllCutting_MF_subset))

#Run pagel's discrete test of correlated evolution
fit.GEN_CircumAllCutting_MF <- fitPagel(GEN_CircumAllCutting_MF_subset_Sample,GEN_CircumAllCutting_Male,GEN_CircumAllCutting_Female)

saveRDS(fit.GEN_CircumAllCutting_MF, "~/fit.GEN_CircumAllCutting_MF_Loose.rds")

## Circumcision and Labial Excision
# Create subsets of data for correlated evolution
GEN_CircumLabExcise_MF_subset <- SCCS_PBM_PA_Loose %>%
  filter(GEN_M_Circumcision %in% c(0, 1),
         GEN_F_LabialExcision %in% c(0, 1)) %>%
  select(Node, GEN_M_Circumcision, GEN_F_LabialExcision)

#Prune SCCS_Tree to  sample
GEN_CircumLabExcise_MF_subset_Sample <- keep.tip(SCCS_Tree, GEN_CircumLabExcise_MF_subset$Node)
#Avoid zero edge lengths
GEN_CircumLabExcise_MF_subset_Sample$edge.length[GEN_CircumLabExcise_MF_subset_Sample$edge.length==0] = 0.000001
#Order the data frame by the order of tree tiplabels
GEN_CircumLabExcise_MF_subset <- GEN_CircumLabExcise_MF_subset[GEN_CircumLabExcise_MF_subset_Sample$tip.label,]

GEN_CircumLabExcise_Male<-setNames(GEN_CircumLabExcise_MF_subset$GEN_M_Circumcision,rownames(GEN_CircumLabExcise_MF_subset))
GEN_CircumLabExcise_Female<-setNames(GEN_CircumLabExcise_MF_subset$GEN_F_LabialExcision,rownames(GEN_CircumLabExcise_MF_subset))

#Run pagel's discrete test of correlated evolution
fit.GEN_CircumLabExcise_MF <- fitPagel(GEN_CircumLabExcise_MF_subset_Sample,GEN_CircumLabExcise_Male,GEN_CircumLabExcise_Female)

saveRDS(fit.GEN_CircumLabExcise_MF, "~/fit.GEN_CircumLabExcise_MF_Loose.rds")


## Amputation
#Prune SCCS_Tree to  sample
AMP_MF_subset_Sample <- keep.tip(SCCS_Tree, AMP_MF_subset$Node)
#Avoid zero edge lengths
AMP_MF_subset_Sample$edge.length[AMP_MF_subset_Sample$edge.length==0] = 0.000001
#Order the data frame by the order of tree tiplabels
AMP_MF_subset <- AMP_MF_subset[AMP_MF_subset_Sample$tip.label,]

AMP_Male<-setNames(AMP_MF_subset$AMP_M,rownames(AMP_MF_subset))
AMP_Female<-setNames(AMP_MF_subset$AMP_F,rownames(AMP_MF_subset))

#Run pagel's discrete test of correlated evolution
fit.AMP_MF_Loose <- fitPagel(AMP_MF_subset_Sample,AMP_Male,AMP_Female)

saveRDS(fit.AMP_MF_Loose, "~/fit.AMP_MF_Loose.rds")

## Dental Modification
#Prune SCCS_Tree to  sample
DEN_MF_subset_Sample <- keep.tip(SCCS_Tree, DEN_MF_subset$Node)
#Avoid zero edge lengths
DEN_MF_subset_Sample$edge.length[DEN_MF_subset_Sample$edge.length==0] = 0.000001
#Order the data frame by the order of tree tiplabels
DEN_MF_subset <- DEN_MF_subset[DEN_MF_subset_Sample$tip.label,]

DEN_Male<-setNames(DEN_MF_subset$DEN_M,rownames(DEN_MF_subset))
DEN_Female<-setNames(DEN_MF_subset$DEN_F,rownames(DEN_MF_subset))

#Run pagel's discrete test of correlated evolution
fit.DEN_MF_Loose <- fitPagel(DEN_MF_subset_Sample,DEN_Male,DEN_Female)

saveRDS(fit.DEN_MF_Loose, "~/fit.DEN_MF_Loose.rds")

## Bone Shaping
#Prune SCCS_Tree to  sample
SHA_MF_subset_Sample <- keep.tip(SCCS_Tree, SHA_MF_subset$Node)
#Avoid zero edge lengths
SHA_MF_subset_Sample$edge.length[SHA_MF_subset_Sample$edge.length==0] = 0.000001
#Order the data frame by the order of tree tiplabels
SHA_MF_subset <- SHA_MF_subset[SHA_MF_subset_Sample$tip.label,]

SHA_Male<-setNames(SHA_MF_subset$SHA_M,rownames(SHA_MF_subset))
SHA_Female<-setNames(SHA_MF_subset$SHA_F,rownames(SHA_MF_subset))

#Run pagel's discrete test of correlated evolution
fit.SHA_MF_Loose <- fitPagel(SHA_MF_subset_Sample,SHA_Male,SHA_Female)

saveRDS(fit.SHA_MF_Loose, "~/fit.SHA_MF_Loose.rds")

## Cranial modification (exactly the same as general bone shaping)
# Create subsets of data for correlated evolution
SHA_Cranial_MF_subset <- SCCS_PBM_PA_Loose %>%
  filter(SHA_M_Cranial %in% c(0, 1),
         SHA_F_Cranial %in% c(0, 1)) %>%
  select(Node, SHA_M_Cranial, SHA_F_Cranial)

#Prune SCCS_Tree to  sample
SHA_Cranial_MF_subset_Sample <- keep.tip(SCCS_Tree, SHA_Cranial_MF_subset$Node)
#Avoid zero edge lengths
SHA_Cranial_MF_subset_Sample$edge.length[SHA_Cranial_MF_subset_Sample$edge.length==0] = 0.000001
#Order the data frame by the order of tree tiplabels
SHA_Cranial_MF_subset <- SHA_Cranial_MF_subset[SHA_Cranial_MF_subset_Sample$tip.label,]

SHA_Cranial_Male<-setNames(SHA_Cranial_MF_subset$SHA_M_Cranial,rownames(SHA_Cranial_MF_subset))
SHA_Cranial_Female<-setNames(SHA_Cranial_MF_subset$SHA_F_Cranial,rownames(SHA_Cranial_MF_subset))

#Run pagel's discrete test of correlated evolution
fit.SHA_Cranial_MF <- fitPagel(SHA_Cranial_MF_subset_Sample,SHA_Cranial_Male,SHA_Cranial_Female)

saveRDS(fit.SHA_Cranial_MF, "~/fit.SHA_Cranial_MF_Loose.rds")
