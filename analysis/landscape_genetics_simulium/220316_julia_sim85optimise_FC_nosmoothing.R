setwd("C:/Users/User/OneDrive - LA TROBE UNIVERSITY/Onchocerciasis/PhD/Genomics")
# Function
optimise_resistance <- function(cont.rast, write.dir = "data/210712_resistanceGA_100/"){
  GA.inputs <- GA.prep(ASCII.dir = cont.rast,
                       Results.dir = write.dir,
                       select.trans = list("A"),
                       method = "LL",
                       # max.scale = par[id,"max_scale"]/cell_size,
                       max.cat=100,
                       max.cont=100,
                       parallel = 10,
                       gaisl = FALSE,
                       # island.pop = 25,
                       # numIslands = 10,
                       run = 50,
                       maxiter = 1000,
                       # quiet = FALSE,
                       # maxiter=1000
                       pop.mult = 30
                       # scale=TRUE,
                       # scale.surfaces = parms$scale
  )
  
  jl.inputs <- jl.prep(n.Pops = length(sample.coords),
                       response = lower(genetic_distance),
                       CS_Point.File = sample.coords,
                       cholmod = T,
                       JULIA_HOME = JULIA_HOME
  )
  # Run optimization
  SS_RESULTS <- SS_optim(jl.inputs = jl.inputs,
                         GA.inputs = GA.inputs)
  return(SS_RESULTS)
}

# Running the function ----------------------------------------------------

# Load packages
library("lme4")
library("ggplot2")
library("sp")
library("raster")
library("lme4")
library(rgdal)
library(ResistanceGA)
library("parallel")
library("doParallel")
library(tictoc)


# Load data ---------------------------------------------------------------
ghana_LG <- readRDS("Code/Project codes/220128_Ghana_simulium_analysis/data/220304_GT_sim85LG.rds")

ghana_sites_selected <- ghana_LG[[2]]
coordinates(ghana_sites_selected) <- ghana_sites_selected[,c("coords.x1", "coords.x2")]
covs_simulium <- stack("Data/Ghana GIS/220228_cov_simulium.grd")
covs_selected <- c("elevation", "isothermality", "SM1315_GT_utm", "FC_GT_utm", "annual_precpitation")

covariates_selected <- covs_simulium[[covs_selected]]
covariates_list <- aggregate(covariates_selected, fact = 2, fun = mean, na.rm = TRUE)

covariates_selected %>% values %>% summary()

reynolds_fst <- ghana_LG[[1]]
genetic_distance <- reynolds_fst/(1-reynolds_fst)

# mat_geo <- ghana_LG[[3]]

## Sample locations
sample.coords <- as(ghana_sites_selected,"SpatialPoints")

# Make a folder to store everything -------------

JULIA_HOME = 'C:/Users/User/AppData/Local/Programs/Julia-1.6.1/bin'


# Flow accumulation -------------------------------------------------------
cont.rast <- covariates_list[[4]]


## Iteration 1
write.dir <- paste0("Code/Project codes/220128_Ghana_simulium_analysis/results/resistanceGA_85/220316_nosmooth__", names(covariates_list)[4],"_1/")  
if(!dir.exists(write.dir)) dir.create(write.dir)

resistance_FC1 <- optimise_resistance(cont.rast = cont.rast, write.dir = write.dir)

FC_resist_1 <- raster(paste0(write.dir,"Results/FC_GT_utm.asc"))


## Iteration 2
write.dir <- paste0("Code/Project codes/220128_Ghana_simulium_analysis/results/resistanceGA_85/220316_nosmooth__", names(covariates_list)[4],"_2/")  

if(!dir.exists(write.dir)) dir.create(write.dir)

resistance_FC2 <- optimise_resistance(cont.rast = cont.rast, write.dir = write.dir)

FC_resist_2 <- raster(paste0(write.dir,"Results/FC_GT_utm.asc"))

## Iteration 3
write.dir <- paste0("Code/Project codes/220128_Ghana_simulium_analysis/results/resistanceGA_85/220316_nosmooth__", names(covariates_list)[4],"_3/")  

if(!dir.exists(write.dir)) dir.create(write.dir)

resistance_FC3 <- optimise_resistance(cont.rast = cont.rast, write.dir = write.dir)

FC_resist_3 <- raster(paste0(write.dir,"Results/FC_GT_utm.asc"))

## Iteration 4
write.dir <- paste0("Code/Project codes/220128_Ghana_simulium_analysis/results/resistanceGA_85/220316_nosmooth__", names(covariates_list)[4],"_4/")  

if(!dir.exists(write.dir)) dir.create(write.dir)

resistance_FC4 <- optimise_resistance(cont.rast = cont.rast, write.dir = write.dir)

FC_resist_4 <- raster(paste0(write.dir,"Results/FC_GT_utm.asc"))

