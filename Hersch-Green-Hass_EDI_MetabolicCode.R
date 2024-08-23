### Genome material cost and metabolic trade-off code

# Load packages 
library(ape)
library(phytools)
library(dplyr)

# Load in mega tree
NutNet_tree <- read.tree("NutNet Tree.tre")

# Load in species lists
all_species <- read.csv("PhotosynthesisNov2023_NOExclusions_SpeciesList.csv")
# Monarda sp. excluded as there are other Monarda species and this specific individual was genus only

#Create sp_list as a character string, not a table, so it can be matched against the tip label character string
sp_list_all <-  all_species$Taxon 

# Match species list to tips on the NutNet tree 
tip_labels = lapply( 1:length(sp_list_all), 
                     function( x ) grep( paste(sp_list_all[[x]], collapse='.+' ),
                                         NutNet_tree$tip.label, ignore.case=T ) )

# Work out which species on your list are not in the tree, so you can investigate
not_found = which( sapply( tip_labels, length )==0 )
sp_list_all[not_found] #0 not found (this is a good thing!)

# Corrections made to "Photosyn_SpeciesList_ALL_4Dec2023.csv":
# Brickellia eupatroides -> Brickellia eupatorioides (spelling error)
# Chrysanthemum leucanthemum -> Chrysanthemum indicum 
# Dicanthelium oligosanthes changed to -> Panicum_oligosanthes 
# Echinaceae angustifolia -> Echinacea angustifolia (spelling error)
# Echinaceae pallida -> Echinacea pallida (spelling error)
# Monarda sp. = EXCLUDED as it is genus only, and we have two other species
# Taraxacum officinale -> Taraxacum_campylodes (changed to closely related species)


# To make substitutions (if you want to include species rather than just leaving them out)
#Create a second row in our table with "phylo" names to match to the tree (you can then have the taxon name and phylo name in your 
# raw data to have the actual name and the name you use for phylogenetic analysis)
all_species$phylo <- all_species$Taxon

# Make sure there are no NAs in the tip labels 
tip_labels <- sapply( tip_labels, function( x ) x[ 1 ] )
not_NA <- !is.na(tip_labels )
tip_labels <- tip_labels[ not_NA ] # TRUE
NutNet_tree$tip.label[ tip_labels ] <- sp_list_all[ not_NA ]

# Construct the tree - i.e., removing all species from the larger tree not in your list
phylogeny_all <- keep.tip(NutNet_tree, tip=tip_labels )

# Save file
write.tree()
data_all <- all_species

#################################### GRASS #####################################
# Load in grass species list
grass_species <- read.csv("Photosyn_SpeciesList_GRASS_4Dec2023.csv")

#Create sp_list as a character string, not a table, so it can be matched against the tip label character string
sp_list_grass <-  grass_species$Taxon # Finally worked as a string... Need the i in front of Taxon

# Match species in your list to tips on the NutNet tree
tip_labels = lapply( 1:length(sp_list_grass), 
                     function(x) grep(paste(sp_list_grass[[x]], collapse='.+'),
                                      NutNet_tree$tip.label, ignore.case=T))

# Work out which species on your list are not in the tree, so you can investigate
not_found = which(sapply(tip_labels, length)==0 )
sp_list_grass[not_found] #0 species not found

# Corrections made to grass species list:
# Dicanthelium oligosanthes changed to -> Panicum_oligosanthes

# To make substitutions (if you want to include species rather than just leaving them out)
#Create a second row in our table with "phylo" names to match to the tree (you can then have the taxon name and phylo name in your 
# raw data to have the actual name and the name you use for phylogenetic analysis)
grass_species$phylo <- grass_species$Taxon

# Make sure there are no NAs in the tip labels 
tip_labels <- sapply(tip_labels, function(x) x[1])
not_NA <- !is.na(tip_labels)
tip_labels <- tip_labels[not_NA]
NutNet_tree$tip.label[ tip_labels] <- sp_list_grass[not_NA]

# Construct the tree - removing all species from the larger tree not in your list
phylogeny_grass <- keep.tip(NutNet_tree, tip=tip_labels )

# Save file
write.tree()
data_grass <- grass_species

#################################### FORB ######################################
# Load in your species lists
forb_species <- read.csv("Photosyn_SpeciesList_FORB_4Dec2023.csv")

#Create sp_list as a character string, not a table, so it can be matched against the tip label character string
sp_list_forb <-  forb_species$Taxon # Finally worked as a string... Need the i in front of Taxon

# Match species in your list to tips on the NutNet tree
tip_labels = lapply(1:length(sp_list_forb), 
                    function(x) grep(paste(sp_list_forb[[x]], collapse='.+'),
                                     NutNet_tree$tip.label, ignore.case=T))

# Work out which species on your list are not in the tree, so you can investigate
not_found = which(sapply(tip_labels, length)==0 )
sp_list_forb[not_found] # 0 species not found

# Corrections made to forb ppecies list:
# Brickellia eupatroides -> Brickellia eupatorioides (spelling error)
# Chrysanthemum leucanthemum -> Chrysanthemum indicum
# Echinaceae angustifolia -> Echinacea angustifolia (spelling error)
# Echinaceae pallida -> Echinacea pallida (spelling error)
# Monarda sp. = EXCLUDED
# Taraxacum officinale -> Taraxacum_campylodes (changed to a closely related species)


# To make substitutions (if you want to include species rather than just leaving them out)
#Create a second row in our table with "phylo" names to match to the tree (you can then have the taxon name and phylo name in your 
# raw data to have the actual name and the name you use for phylogenetic analysis)
forb_species$phylo <- forb_species$Taxon

# Make sure there are no NAs in the tip labels 
tip_labels <- sapply(tip_labels, function(x) x[1])
not_NA <- !is.na(tip_labels)
tip_labels <- tip_labels[not_NA]
NutNet_tree$tip.label[ tip_labels] <- sp_list_forb[not_NA]

# Construct the tree - removing all species from the larger tree not in your list
phylogeny_forb <- keep.tip(NutNet_tree, tip=tip_labels )

# Save file
write.tree()
data_forb <- forb_species
#Forb tree= pruned!

################ Testing for phylogenetic signal of GS #########################
# Load in packages
library(ape)
library(picante)
library(phytools)
library(tidyverse)

# Load the tree and data #
new_tree <- read.tree("Photosyn_AllLifeform_tree_Feb2024.tre")
new_tree

# Breaking the tree information apart into tip label, node/branch lengths, etc. 
str(new_tree)
new_tree$edge # Provides a list of beginning and end node #s
new_tree$tip.label # View each individual tip label (taxa) in the tree
new_tree$edge.length # View edge length (length of each edge in relation to the root of the tree)
new_tree$Nnode # Number of internal nodes on the tree 

#Building a phylogenetic tree image with R
plot(new_tree, edge.width = 1.0, label.offset = 0.75, 
     cex = 0.25, type = "cladogram") 

# Check to see if tree is binary. 
is.binary(new_tree)

# Check if the tree is rooted
is.rooted(new_tree)

# Change any zero-length branches to one-ten-thousandth of the tree size (Winternitz, 2016)
new_tree$edge.length[new_tree$edge.length==0] <- max(nodeHeights(new_tree))*1e-4

# Make sure tree tips match the data EXACTLY 
new_tree$tip.label <- gsub(" ", "_", new_tree$tip.label) # Replaces white space with _

# We want to test phylogenetic signal using Pagel's lambda
# This has the smallest Type I error, and small Type II errors for phylogenies >20 species (Munkemuller et al. 2012)

# Testing Pagel's lambda and Bloomberg's K

#Load in library
library(phytools) # measure phylo signal

# IMPORT DATA HERE
phydata <- read.csv("Hersch-Green and Hass_metabolicMS_EDI.csv")


# Rearrange the data so you have separate datasets for grasses and forbs
phydata_grass <- phydata %>% filter(Lifeform_2== "grass") 
phydata_forb <- phydata %>% filter(Lifeform_2== "allforb") 
# Will use "phydata" to test for signal in all lifeforms (grasses, forbs)

########################## PHYLOGENETIC SIGNAL #################################
#### Test for signal of GS in the overall model

# Isolate the species and their respective GS values. 
# The columns we are isolating from "phydata" are Taxa, GS, and logGS (this one is optional, don't need to test the log-transformed GS as untransformed GS will do)
# ROWS MAY CHANGE WITH VARYING DATA SHEETS, BE SURE TO CHECK YOUR DATA BEFORE CONTINUING

#All functional groups (grass, forbs)
unique_sp_ALL <- unique(phydata[,c(4, 12)], row.names=1) 

#Grass
unique_GS_grass <- unique(phydata_grass[,c(4, 12)], row.names=1) 

#Forb
unique_GS_forb <- unique(phydata_forb[,c(4, 12)], row.names=1) 

#Now create a list of named numbers (the species and the GS/log GS value) as a vector not a table
# ALL FUNCTIONAL GROUPS
GS <- setNames(unique_sp_ALL$GS, unique_sp_ALL$Taxa)

#GRASS
GS_grass <- setNames(unique_GS_grass$GS, unique_GS_grass$Taxa)

#FORB
GS_forb <- setNames(unique_GS_forb$GS, unique_GS_forb$Taxa)

# Now use phylosig with the tree and your vector of named GS values for each FG you want to test
############################## Pagel's lambda  ###############################
#ALL FUNCTIONAL GROUPS
phylosig(new_tree, logGS, method = "lambda", test = TRUE) 


#GRASS
phylosig(new_tree, logGS_grass, method = "lambda", test = TRUE) 

#FORB
phylosig(new_tree, logGS_forb, method = "lambda", test = TRUE) 

############################### Bloomberg's K  #################################
# Blomberg's K, based on 10,000 randomizations #
# ALL FUNCTIONAL GROUPS
phylosig(new_tree, logGS, method = "K", test = T, nsim= 10000)

# GRASS
phylosig(new_tree, logGS_grass, method = "K", test = T, nsim= 10000)

# FORB
phylosig(new_tree, logGS_forb, method = "K", test = T, nsim= 10000)

############ Phylogenetic Generalized Linear Mixed Model (PGLMM) ###############
# See Ch. 4 (https://leanpub.com/correlateddata) for more information/examples

# Load in libraries

library(phyr)
library(ape)
library(geiger)
library(nlme)
library(phytools)
library(ggplot2)
library(lme4)
library(dplyr)
library(afex)
library(lmerTest)
library(mosaic)
library(tidyverse)

# Mixed-effects model with the phylogeny as a covariate argument in the pglmm function
# This function gives you p-values for random effects and is easier to interpret.

#Load in data
bigdata <- read.csv("Hersch-Green and Hass_metabolicMS_EDI.csv") 

# Read in data for forb and grass trees
forb_tree <- read.tree("Photosyn_forb_tree_4Dec2023.tre")
grass_tree <- read.tree("Photosyn_grass_tree_4Dec2023.tre")
all_tree <- read.tree("Photosyn_AllLifeform_tree_Feb2024.tre")

forb_data <- filter(bigdata, Lifeform_2== "allforb") 
grass_data <- filter(bigdata, Lifeform_2== "grass") 

#Create a correlation matrix that explains the expected correlation from root to tip under Brownian motion.

# The expected correlation is proportional to the amount of evolutionary history, from root to tips that these species share through common descent.
 
# Data is essentially re-weighted per species

#For all species within the phylogeny
vv_ALL <- vcv(all_tree, model= "Brownian", corr= TRUE)
species_to_keep <- c("Achillea_millefolium", "Ambrosia_artemisiifolia", "Ambrosia_psilostachya", "Ambrosia_trifida", "Amorpha_canescens",
                     "Artemisia_ludoviciana", "Asclepias_syriaca", "Baptisia_alba", "Baptisia_australis", "Brassica_sp.",           
                     "Brickellia_eupatorioides", "Centaurea_stoebe",  "Chamaecrista_fasciculata" ,       
                     "Chenopodium_album",  "Chrysanthemum_indicum",   "Clinopodium_vulgare",         
                     "Coreopsis_tripteris",  "Daucus_carota", "Dalea_candida",              
                     "Echinacea_angustifolia", "Echinacea_pallida",  "Erigeron_annuus",             
                     "Erigeron_sp.",  "Erigeron_strigosus",    "Eryngium_yuccifolium",        
                     "Euthamia_graminifolia", "Fragaria_vesca",  "Galium_aparine",             
                     "Geranium_carolinianum",   "Helianthus_annuus" ,   "Lespedeza_capitata",       
                     "Liatris_aspera" , "Monarda_fistulosa" ,         
                     "Oenothera_biennis",   "Oxalis_stricta",   "Physalis_pumila",      
                     "Plantago_lanceolata",  "Potentilla_recta",            
                     "Ratibida_pinnata", "Rudbeckia_hirta",  "Ruellia_humilis"  ,           
                     "Rumex_acetosella", "Salvia_azurea"     ,          
                     "Silene_latifolia", "Solanum_carolinense",         
                     "Solidago_altissima",  "Solidago_canadensis" ,        
                     "Solidago_gigantea", "Solidago_missouriensis", "Solidago_nemoralis"   ,       
                     "Solidago_rigida",              "Solidago_speciosa",  "Symphyotrichum_ericoides",    
                     "Symphyotrichum_oblongifolium", "Symphyotrichum_pilosum",     
                     "Taraxacum_campylodes" ,        "Tragopogon_dubius",       "Trifolium_pratense",
                     "Vernonia_baldwinii", "Veronica_sp.", "Vicia_grandiflora", "Vicia_sativa", "Viola_sp.", "Andropogon_gerardii",  
                     "Anthoxanthum_odoratum",   "Bouteloua_curtipendula", 
                     "Bromus_inermis",        "Dactylis_glomerata",      "Elymus_canadensis",
                     "Elymus_repens",         "Festuca_arundinacea",     "Panicum_capillare",
                     "Panicum_oligosanthes",  "Panicum_virgatum",        "Phleum_pratense",
                     "Poa_pratensis",         "Schizachyrium_scoparium", "Setaria_pumila",
                     "Sorghastrum_nutans",    "Sorghum_halpense",        "Sporobolus_heterolepis")



vv <- vcv(forb_tree, model= "Brownian", corr= TRUE)
species_to_keep <- c("Achillea_millefolium",        "Ambrosia_artemisiifolia",      "Ambrosia_psilostachya", "Ambrosia_trifida", "Amorpha_canescens",
                     "Artemisia_ludoviciana", "Asclepias_syriaca", "Baptisia_alba", "Baptisia_australis", "Brassica_sp.",           
                     "Brickellia_eupatorioides", "Centaurea_stoebe",  "Chamaecrista_fasciculata" ,       
                     "Chenopodium_album",  "Chrysanthemum_indicum",   "Clinopodium_vulgare",         
                     "Coreopsis_tripteris",  "Daucus_carota", "Dalea_candida",              
                     "Echinacea_angustifolia", "Echinacea_pallida",  "Erigeron_annuus",             
                     "Erigeron_sp.",  "Erigeron_strigosus",    "Eryngium_yuccifolium",        
                     "Euthamia_graminifolia",     "Fragaria_vesca",              
                     "Galium_aparine",             
                     "Geranium_carolinianum",   "Helianthus_annuus" ,   "Lespedeza_capitata",       
                     "Liatris_aspera" , "Monarda_fistulosa" ,         
                     "Oenothera_biennis",   "Oxalis_stricta",   "Physalis_pumila",      
                     "Plantago_lanceolata",  "Potentilla_recta",            
                     "Ratibida_pinnata", "Rudbeckia_hirta",  "Ruellia_humilis"  ,           
                     "Rumex_acetosella", "Salvia_azurea"     ,          
                     "Silene_latifolia", "Solanum_carolinense",         
                     "Solidago_altissima", "Solidago_canadensis" ,        
                     "Solidago_gigantea",  "Solidago_missouriensis", "Solidago_nemoralis"   ,       
                     "Solidago_rigida", "Solidago_speciosa",  "Symphyotrichum_ericoides",    
                     "Symphyotrichum_oblongifolium", "Symphyotrichum_pilosum",     
                     "Taraxacum_campylodes" ,        "Tragopogon_dubius",       "Trifolium_pratense",
                     "Vernonia_baldwinii", "Veronica_sp.", 
                     "Vicia_grandiflora",  "Vicia_sativa", "Viola_sp.")

# Get the row and column indices for the selected FORB species
forb_species_traits <- read.csv("PhotosynthesisNov2023_NoExclusions_TraitForb.csv")
species_indices <- which(forb_species_traits$Taxa %in% species_to_keep)
species_indices <- as.matrix(species_indices)
is.matrix(species_indices) #TRUE 

# Subset the covariance matrix based on the selected species
reduced_cov_matrix <- cov(species_indices) # Now you can set up PGLMM for forbs

#Double checking everything; looking at tree for forbs, checking species names
plot(forb_tree)
name.check(forb_tree, forb_species)

#Order tip labels alphabetically
unique_forb <- sort(unique(forb_tree$tip.label))
unique_forb #number of unique forb species

#Ordering data
unique_forb2 <- sort(unique(forb_species_traits$Taxa))

#Renaming forb data
forb_data_final <- forb_species_traits

############################## GRASSES #########################################
vv_g <- vcv(grass_tree, model= "Brownian", corr= TRUE)
species_to_keep_g <- c("Andropogon_gerardii",   "Anthoxanthum_odoratum",   "Bouteloua_curtipendula", 
                       "Bromus_inermis",        "Dactylis_glomerata",      "Elymus_canadensis",
                       "Elymus_repens",         "Festuca_arundinacea",     "Panicum_capillare",
                       "Panicum_oligosanthes",  "Panicum_virgatum",        "Phleum_pratense",
                       "Poa_pratensis",         "Schizachyrium_scoparium", "Setaria_pumila",
                       "Sorghastrum_nutans",    "Sorghum_halpense",        "Sporobolus_heterolepis")

# Get the row and column indices for the selected GRASS species
grass_species_traits <- read.csv("PhotosynthesisNov2023_NOExclusions_TraitGrass.csv")
species_indices_g <- which(grass_species_traits$Taxa %in% species_to_keep_g)
species_indices_g <- as.matrix(species_indices_g)
is.matrix(species_indices_g) #TRUE 

# Subset the covariance matrix based on the selected species
reduced_cov_matrix_g <- cov(species_indices_g) # Now you can set up PGLMM for grasses

# Order tip labels alphabetically
unique_grass <- sort(unique(grass_tree$tip.label))
unique_grass #number of unique grass species

#Renaming the grass data
grass_data_final <- grass_species_traits 
################################################################################
############################ PGLMM MODELS ######################################

#Load in library
library(phytools)

# Testing the effect of GS, lifeform, treatment, 3-way interactions, and random effects of site and block within site, on cellular nutrient content or C, N, P. 

model1 <- pglmm(SQRT_C_per_Cell_ng ~ LOG_GS + Lifeform_2 + Treatment + LOG_GS*Lifeform_2*Treatment +
                      (1|Taxa__) + (1|Block@Site) + (1|Site), 
                    data= bigdata, 
                    cov_ranef = list(Taxa= vv_ALL),
                    family="gaussian") 
model1 #view table after cov matrix has been standardized by R

================================================================
model2 <- pglmm(Sqrt_N_per_Cell_ng ~ LOG_GS + Lifeform_2 + Treatment + LOG_GS*Lifeform_2*Treatment +
                  (1|Taxa__) + (1|Block@Site) + (1|Site), 
                data= bigdata, 
                cov_ranef = list(Taxa= vv_ALL),
                family="gaussian") 
model2 #view table after cov matrix has been standardized by R

================================================================
model3 <- pglmm(SQRT_.P_per_cell_ng ~ LOG_GS + Lifeform_2 + Treatment + LOG_GS*Lifeform_2*Treatment +
                  (1|Taxa__) + (1|Block@Site) + (1|Site), 
                data= bigdata, 
                cov_ranef = list(Taxa= vv_ALL),
                family="gaussian")
model3 #view table after cov matrix has been standardized by R

================================================================
================================================================
# Testing how stomata size and stomata density are influenced by GS, lifeform, their interaction, and random effects of site and block nested within site

model4 <- pglmm(SQRT_StomataDensity ~ LOG_GS + Lifeform_2 + LOG_GS*Lifeform_2 +
                  (1|Taxa__) + (1|Block@Site) + (1|Site), 
                data= bigdata, 
                cov_ranef = list(Taxa= vv_ALL),
                family="gaussian") 
model4 #view table after cov matrix has been standardized by R

================================================================
model5 <- pglmm(Stomata_Size_avg ~ LOG_GS + Lifeform_2 + LOG_GS*Lifeform_2 +
                  (1|Taxa__) + (1|Block@Site) + (1|Site), 
                data= bigdata, 
                cov_ranef = list(Taxa= vv_ALL),
                family="gaussian") 
model5 #view table after cov matrix has been standardized by R

################################################################################

# Testing the effets of GS, lifeform, treatment, their 3-way interaction, and random effects of site and block nested within site on cellular nutrients per mg of leaf tissue

# Nitrogen mg leaf tissue
modelA <- pglmm(Nmg.mgleaftissue ~ LOG_GS + Lifeform + Treatment + LOG_GS*Lifeform*Treatment +
                  (1|Taxa__) + (1|Block@Site) + (1|Site), 
                data= bigdata, 
                cov_ranef = list(Taxa= vv_ALL),
                family="gaussian") 
modelA #view table after cov matrix has been standardized by R

================================================================
# Carbon mg leaf tissue
modelB <- pglmm(Cmg.mgleaftissue ~ LOG_GS + Lifeform + Treatment + LOG_GS*Lifeform*Treatment +
                  (1|Taxa__) + (1|Block@Site) + (1|Site), 
                data= bigdata, 
                cov_ranef = list(Taxa= vv_ALL),
                family="gaussian") 
modelB #view table after cov matrix has been standardized by R

================================================================
# Phosphorus mg leaf tissue
modelC <- pglmm(mgP.mgleaftissue ~ LOG_GS + Lifeform_2 + Treatment + LOG_GS*Lifeform_2*Treatment +
                  (1|Taxa__) + (1|Block@Site) + (1|Site), 
                data= bigdata, 
                cov_ranef = list(Taxa= vv_ALL),
                family="gaussian") 
modelC #view table after cov matrix has been standardized by R

# END FILE #
