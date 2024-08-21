# genomesize_sem.r
## script to analyze the impact of genome size on physiology under different environmental
## contexts
## from Hersch-Green data

## load libraries
library(piecewiseSEM)
library(lme4)
library(car)

## load data
physgs_data <- read.csv('../data/PhysiologyGS_Final.csv')
head(physgs_data)

## build SEM
### simplify column names
physgs_data$E <- physgs_data$E..mmol.m.2.s.1. 

### binaries for soil nutrient treatments
physgs_data$ntrt_binary[physgs_data$N_added == 'N_YES'] <- 1 
physgs_data$ntrt_binary[physgs_data$N_added == 'N_NO'] <- 0 # binaries for soil nutrient treatments
physgs_data$ptrt_binary[physgs_data$P_added == 'P_YES'] <- 1 # binaries for soil nutrient treatments
physgs_data$ptrt_binary[physgs_data$P_added == 'P_NO'] <- 0 # binaries for soil nutrient treatments

### create composite variable for stomata
physgs_data$stomdens_stomsize <- physgs_data$Stomata_Size_avg * physgs_data$Stomata_Density

### transform and scale variables
physgs_data$sc_amax <- as.numeric(scale(sqrt(physgs_data$Amax)))
physgs_data$sc_e <- as.numeric(scale(sqrt(physgs_data$E)))
physgs_data$sc_ncell <- as.numeric(scale(sqrt(physgs_data$N_PerCell.mg)))
physgs_data$sc_pcell <- as.numeric(scale(sqrt(physgs_data$P_PerCell.mg)))
physgs_data$sc_gs <- as.numeric(scale(log(physgs_data$GS)))
physgs_data$sc_stomsize <- as.numeric(scale(physgs_data$Stomata_Size_avg))
physgs_data$sc_stomdens <- as.numeric(scale(physgs_data$Stomata_Density))
physgs_data$sc_mat <- as.numeric(scale(physgs_data$MAT_v2))
physgs_data$sc_map <- as.numeric(scale(physgs_data$MAP_v2))
physgs_data$sc_wue <- as.numeric(scale(physgs_data$WUE))
physgs_data$sc_stomdens_stomsize <- as.numeric(scale(physgs_data$stomdens_stomsize))

### create data subsets for plant lifeforms
physgs_data_forb <- subset(physgs_data, Lifeform == 'FORB')
physgs_data_grass <- subset(physgs_data, Lifeform == 'GRASS')

### remove all rows with relevant NAs for pcell only (not ncell)
physgs_data_forb_narm_ncell <- physgs_data_forb[complete.cases(physgs_data_forb[, c(45:47,49:54)]), ]
physgs_data_grass_narm_ncell <- physgs_data_grass[complete.cases(physgs_data_grass[, c(45:47,49:54)]), ]
physgs_data_all_narm_ncell <- physgs_data[complete.cases(physgs_data[, c(45:47,49:54)]), ]

### fit models without pcell
#### forb
physgs_forb_sem_ncell <- psem(
  
  amax = lmer(sc_amax ~ sc_stomdens + 
                sc_stomsize +
                sc_ncell +
                sc_mat + sc_map +
                sc_gs +
                ntrt_binary + ptrt_binary +
                (1|Taxa),
              data = physgs_data_forb_narm_ncell),
  
  transpiration = lmer(sc_e ~ sc_stomdens + 
                         sc_stomsize +
                         sc_ncell +
                         sc_mat + sc_map +
                         sc_gs +
                         ntrt_binary + ptrt_binary +
                         (1|Taxa), 
                       data = physgs_data_forb_narm_ncell),
  
  stomata_size = lmer(sc_stomsize ~ sc_gs +
                       (1|Taxa), 
                     data = physgs_data_forb_narm_ncell),
  
  stomata_density = lmer(sc_stomdens ~ sc_gs +
                          (1|Taxa), 
                        data = physgs_data_forb_narm_ncell),
  
  ncell = lmer(sc_ncell ~ sc_gs +
                 ntrt_binary + ptrt_binary +
                 (1|Taxa), 
               data = physgs_data_forb_narm_ncell),
  
  genome_size = lmer(sc_gs ~ sc_mat + sc_map +
              (1|Taxa), 
            data = physgs_data_forb_narm_ncell)
  
)

physgs_forb_sem_summary_ncell <- summary(physgs_forb_sem_ncell)

#### grass
physgs_grass_sem_ncell <- psem(
  
  amax = lmer(sc_amax ~ sc_stomdens + 
                sc_stomsize +
                sc_ncell +
                sc_mat + sc_map +
                sc_gs +
                ntrt_binary + ptrt_binary +
                (1|Taxa),
              data = physgs_data_grass_narm_ncell),
  
  transpiration = lmer(sc_e ~ sc_stomdens + 
                         sc_stomsize +
                         sc_ncell + 
                         sc_mat + sc_map +
                         sc_gs +
                         ntrt_binary + ptrt_binary +
                         (1|Taxa), 
                       data = physgs_data_grass_narm_ncell),
  
  stomata_size = lmer(sc_stomsize ~ sc_gs +
                       (1|Taxa), 
                     data = physgs_data_grass_narm_ncell),
  
  stomata_density = lmer(sc_stomdens ~ sc_gs +
                          (1|Taxa), 
                        data = physgs_data_grass_narm_ncell),
  
  ncell = lmer(sc_ncell ~ sc_gs +
                 ntrt_binary + ptrt_binary +
                 (1|Taxa), 
               data = physgs_data_grass_narm_ncell),

  genome_size = lmer(sc_gs ~ sc_mat + sc_map +
               (1|Taxa), 
            data = physgs_data_grass_narm_ncell)
  
)

physgs_grass_sem_summary_ncell <- summary(physgs_grass_sem_ncell)

## write output
# write.csv(physgs_forb_sem_summary_ncell[[6]], '../output/physgs_forb_sem_summary_ncell.csv', row.names = F)
# write.csv(physgs_grass_sem_summary_ncell[[6]], '../output/physgs_grass_sem_summary_ncell.csv', row.names = F)

# write.csv(physgs_forb_sem_summary_ncell[[7]], '../output/physgs_forb_sem_summary_ncell_r2.csv', row.names = F)
# write.csv(physgs_grass_sem_summary_ncell[[7]], '../output/physgs_grass_sem_summary_ncell_r2.csv', row.names = F)
