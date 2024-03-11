####### Setting work directory and librarying packages ######

setwd("C:/Users/danid/OneDrive - Stanford/Otago/Lab Poulin/UORG/R/RDatas")
library(HPprediction)
library(readxl)
library(dplyr)
library(tidyr)
library(tibble)
library(malaviR)
library(data.table)
library(rgdal)
library(gdalUtils)
library(sf)
library(sp)
library(spdep)
library(ggplot2)
library(raster)
library(shapefiles)
library(letsR)
library(picante)
library(ape)
library(factoextra)
library(patchwork)
library(brms)
library(bipartite)

##### MalAvi Compiling Data #####


DF <- extract_table(table = "Hosts and Sites Table") #downloading malavi data
DF <- filter(DF,coordinates != "")
DF<- filter(DF, Host_Environment == "Wild")
DFB <- extract_table(table = "Grand Lineage Summary") #downloading malavi sequences
DFB <- DFB[c(2,25)]
DFB <- distinct(DFB) #checking for duplicate lineages
length(unique(DFB$sequence))
length(unique(DFB$Lineage_Name))
DFB1 <- with(DFB, names(table(sequence)[table(sequence) > 1])) #isolating the duplicated lineages
DFB1 <- DFB[DFB$sequence%in% DFB1, ]
DFB1$sequence <- factor(DFB1$sequence)
DFB1 <- transform(DFB1,ParasiteID=as.numeric(factor(DFB1$sequence)))
DFB1$ParasiteID <- paste("Parasite", DFB1$ParasiteID)

#removing duplicated names
DFB$Lineage_Name[DFB$Lineage_Name  == "ATLPIL05"] <-  "ATLPIL04"
DFB$Lineage_Name[DFB$Lineage_Name  == "CYACAE09"] <- "CYACAE08"
DFB$Lineage_Name[DFB$Lineage_Name  == "TURASS05"] <- "JUNPHA04"
DFB$Lineage_Name[DFB$Lineage_Name  == "NEWBR04"] <- "NEWBR01"
DFB$Lineage_Name[DFB$Lineage_Name  == "PERATE09"] <- "PERATE07"
DFB$Lineage_Name[DFB$Lineage_Name  == "TURAMA06"] <- "TURALB06"
DFB$Lineage_Name[DFB$Lineage_Name  == "TURMER12"] <- "TURMER11"
DFB$Lineage_Name[DFB$Lineage_Name  == "TURMER13"] <- "TURMER11"

length(unique(DFB$sequence))
length(unique(DFB$Lineage_Name))                

#filtering unnecessary information from dataframe
DF1 <- DF[-c(1,3,9,10,11,12,18,19,20,21)]
length(unique(DF1$Lineage_Name))

#removing duplicated names
DF1$Lineage_Name[DF1$Lineage_Name  == "ATLPIL05"] <-  "ATLPIL04"
DF1$Lineage_Name[DF1$Lineage_Name  == "CYACAE09"] <- "CYACAE08"
DF1$Lineage_Name[DF1$Lineage_Name  == "TURASS05"] <- "JUNPHA04"
DF1$Lineage_Name[DF1$Lineage_Name  == "NEWBR04"] <- "NEWBR01"
DF1$Lineage_Name[DF1$Lineage_Name  == "PERATE09"] <- "PERATE07"
DF1$Lineage_Name[DF1$Lineage_Name  == "TURAMA06"] <- "TURALB06"
DF1$Lineage_Name[DF1$Lineage_Name  == "TURMER12"] <- "TURMER11"
DF1$Lineage_Name[DF1$Lineage_Name  == "TURMER13"] <- "TURMER11"

# removing odd characters from coordinates

DF1$coordinates <- gsub("ÃÂ°", ".", DF1$coordinates)
DF1$coordinates <- gsub("ÃÂ´", "", DF1$coordinates)
DF1$coordinates <- gsub("'", "", DF1$coordinates)


#separing coordinates
DF1 <- separate(DF1, coordinates, c("Latitude", "Longitude"), sep = ",")
DF1$Latitude <- substr(DF1$Latitude, 1, nchar(DF1$Latitude) - 6)
DF1$Longitude <- substr(DF1$Longitude, 1, nchar(DF1$Longitude) - 6)

DF1[,11:12] <- sapply(DF1[,11:12], as.numeric)
sum(is.na(DF1$Latitude))
sum(is.na(DF1$Longitude))
DF1 <- na.omit(DF1)

DF1$LatLong <- paste(DF1$Latitude, DF1$Longitude, sep = "_")
length(unique(DF1$LatLong))
length(unique(DF1$Lineage_Name))
rm(DFB1, DFB)


##### compiling MalAvi dataset with other datasets ######

DF2 <- read_excel("Dufour 2020.xlsx")
DF3 <- read_excel("OpenTraits.xlsx")

DF1$species <- gsub(" ", "_", DF1$species)

data <- inner_join(DF1, dplyr::select(DF2, Match_jetz_birdlife, strategy_3, distance_4),
                   by = c("species" = "Match_jetz_birdlife"))
data1 <- inner_join(data, dplyr::select(DF3, Tree_name, Body_mass_log, Territoriality , Range_Size),
                    by = c("species" = "Tree_name"))

data2 <- with(data1, names(table(LatLong)[table(LatLong) > 9]))
data2 <- data1[data1$LatLong%in% data2, ]
data2$LatLong <- factor(data2$LatLong)

length(unique(data2$LatLong))
length(unique(data2$species))
length(unique(data2$Lineage_Name))

data2 <- transform(data2,LocID=as.numeric(factor(data2$LatLong)))
data2$LocID <- paste("Loc", data2$LocID)
data2 <- transform(data2,ParasiteID=as.numeric(factor(data2$Lineage_Name)))
data2$ParasiteID <- paste("Parasite", data2$ParasiteID)
write_csv2(data2, file = "Table1.csv") #creating file

dados <- read.csv2("Table1.csv") #data with coordinates converted
dados$Body_mass <- dados$Body_mass_log
dados$Body_mass_log <- NULL
save(dados, file = "dataframe.RData") 

###### Estimating Local Bird Richness ######

vector <- as.vector(unique(dados$species))
vector <- gsub("_", " ", vector)
write(vector, file = "species.txt")

unzip("Drop", exdir = "BOTW_2022_2")

# full path to the database required
path <- "BOTW_2022_2/Drop"

# Get a list of layer names in your geodatabase
layers <- ogrListLayers(path)

# list all feature classes in a file geodatabase
library(rgdal)
subset(ogrDrivers(), grepl("GDB", name))
ogrListLayers(path)

ogrDrivers()

# extract info
info <- ogrInfo(path, require_geomType= "wkbPolygon")
info

##### Reading ESRI polygon #####

load("UORG.RData")

# Reading filtered dataset (filtered in QGIS)

AllBirds <- readOGR(dsn = path, layer = layers) #, layer = layers) TRY REMOVE THIS 
names(AllBirds)
summary(AllBirds)

length(unique(AllBirds@data$sciname))

AllBirds1 <- sf::st_as_sf(AllBirds)
class(AllBirds1)
st_crs(AllBirds1)

ogrListLayers(path)

AllBirds_clean <- AllBirds[c(1:13201, 13203:17414),] # removing problematic layer

#creating presence/absence tableI
MB <- letsR::lets.presab(AllBirds, resol = 2)
summary(MB)

map(database = "world", bg = "lightblue", fill = TRUE, col = "darkgreen", )
points(dados$Longitude, dados$Latitude, pch = 16, col = "red")

plot(MB, xlab = "Longitude", ylab = "Latitude", main = "Bird species richness")

rm(AllBirds, AllBirds1, AllBirds_clean)

save.image("UORG.RData")


##### Bayesian forecasting of missing host-parasite links #####

#### Preparing Matrices ####

load("UORG.RData")
load("MB.RData")

MB_data <- MB$Presence_and_Absence_Matrix
class(MB_data)

# creating datadframe 
MB_data <- as.data.frame(MB_data)

# separating coordinates values

longitude <- MB_data$`Longitude(x)`
latitude <- MB_data$`Latitude(y)`

colnames(MB_data)[1] <- "Longitude"
colnames(MB_data)[2] <- "Latitude"

# subseting the matrix 5x5 degrees

occmatrix1 <- list()

Lon <- -180
Lat <- -90
i <- 1
r <- 0
for (i in 1:2592) {
  # Create a matrix for the current cell
  occmatrix1[[i]] <- MB_data[(MB_data$Longitude >= Lon & MB_data$Longitude < Lon + 5 & MB_data$Latitude >= Lat & MB_data$Latitude < Lat + 5), ]
  
  # Change value of Lon
  Lon <- Lon + 5
  
  # Change value of Lat when Lon >= 170
  if (Lon >= 175) {
    Lat <- Lat + 5
    Lon <- -180
  }
  # print loop position
  print(i)
} 

#save(occmatrix1, file = "occmatrix1.RData")

#save.image("UORG.RData")

#### Running HP Prediction models and calculating network metrics ####

# Loop with all matrices

occmatrix2 <- occmatrix1

library(bipartite)

lload("occmatrix1.RData")
load("dados.RData")
load("random_trees2.RData")


z_scores1 <- data.frame(stringsAsFactors = FALSE)
k <- 1
Lon <- -180
Lat <- -90
pred <- data.frame(stringsAsFactors = FALSE)

# Creating list with locations and host species
for (k in 1:length(occmatrix2)) {
  # print 
  print(k)
  
  occlist <- occmatrix2[[k]] %>%
    tidyr::pivot_longer(cols = -c(Latitude, Longitude),
                        names_to = "species",
                        values_to = "presence") %>%
    filter(presence >= 1) 
  
  occlist$presence <- NULL
  
  # Formatting host names to match phylogeny
  occlist$species <- gsub(" ","_", occlist$species)
  occlist$Lineage_Name <- "NA"
  
  # Filtering  parasite-host data per geographical location
  dados_com <- dados %>%
    filter(Longitude >= Lon & Longitude < Lon + 5) %>%
    filter(Latitude >= Lat & Latitude < Lat + 5)
  
  # Change value of Lon
  Lon <- Lon + 5
  
  # Change value of Lat when Lon >= 170
  if (Lon >= 175) {
    Lat <- Lat + 5
    Lon <- -180
  }
  
  dados_com <- dados_com[c(1,6)]
  
  dados_com <- rbind(dados_com, occlist[c(3,4)])
  
  # Check a condition
  if (length(dados_com$Lineage_Name[dados_com$Lineage_Name != "NA"]) < 20) {
    # Condition is true, move to the next iteration
    next
  } else {
    
    # Creating binary interaction matrix
    com <- table(dados_com$species, dados_com$Lineage_Name)
    #com[com>1] <- 1
    com <- as.matrix(unclass(com), nrow=nrow(com), ncol=ncol(com))
    
    # loading phylogeny and pruning to hosts in the interaction matrix
    index <- sample(length(random_trees2), 1)
    birds_tree <- as.phylo(random_trees2[[index]])
    
    # Getting consensus tree out of random 100 phylogenetic trees
    
    birds_tree <- drop.tip(birds_tree, birds_tree$tip.label[!birds_tree$tip.label%in%rownames(com)])
    
    # merge the tree and interaction matrix
    cleaned <- network_clean(com, birds_tree, 'full')
    com <- cleaned$Z                         # cleaned binary interaction matrix
    tree <- cleaned$tree                     # cleaned tree
    
    ## Running the model
    
    out <- network_est(Z = com, tree = tree, slices = 1000, model.type = 'full')
    
    # calculating difference between before and after predictions
    pred <- rbind.data.frame(paste(out/com))
    
    ##### Calculating Network Metrics ######
    postocc <-  out$Z
    nas <- colnames(postocc) == "NA"
    postocc <- postocc[, colnames(postocc) != "NA"]
    
    # Network Metrics
    
    Comm_properties <-networklevel(postocc, index=c("NODF", "interaction evenness", "H2"), 
                                   H2_integer=FALSE)
    Comm_properties$Community <- paste("Com", k, sep = " ")  
    Comm_properties <- as.data.frame(Comm_properties)
    
    # Parasite Metrics
    Para_properties <-grouplevel(postocc, index=c("extinction slope", "niche overlap"), level = "higher")
    Para_properties$Community <- paste("Com", k, sep = " ")  
    Para_properties <- as.data.frame(Para_properties)
    
    # joining datasets
    
    data <- merge(Para_properties, Comm_properties, by = "Community")
    
    # Calculating Z-Scores for parasite metrics
    
    null_model <- bipartite::nullmodel(postocc, method = "vaznul")
    
    ##### Use Loop to calculate it #####
    
    j <- 1
    
    NullScores <- data.frame(niche.overlap.HL = numeric(), 
                             extinction.slope.HL = numeric(), 
                             nodf = numeric(),
                             H2 = numeric(),
                             Community = factor(),
                             stringsAsFactors = FALSE)
    
    for (j in 1:length(null_model)) {
      # Calculate group-level metrics for random values
      random_metrics1 <- networklevel(null_model[[j]], index= c("NODF" , "H2"), 
                                      H2_integer=FALSE)
      random_metrics1 <- t(random_metrics1)
      random_metrics1 <- data.frame(random_metrics1, stringsAsFactors = FALSE)
      
      
      # Calculate group-level metrics for random values
      random_metrics2 <- grouplevel(null_model[[j]], index=c("extinction slope", "niche overlap"), level = "higher")
      random_metrics2 <- t(random_metrics2)
      random_metrics2 <- data.frame(random_metrics2, stringsAsFactors = FALSE)
      
      # merge datasets
      
      random_metrics <- merge(random_metrics1, random_metrics2)
      
      # Adding values to dataframe
      NullScores <- rbind.data.frame(NullScores, random_metrics)
      
      # summing J
      j <- j + 1
      
    }
    
    #Calculate mean and standard deviation for each community
    mean_random <- colMeans(NullScores)
    mean_random <- t(mean_random)
    mean_random <- as.data.frame(mean_random)
    
    sd_random <- apply(NullScores, 2, sd)
    sd_random <- t(sd_random)
    sd_random <- as.data.frame(sd_random)
    
    # Calculate Z-scores for each community
    niche <-  (Para_properties$niche.overlap.HL - mean_random$niche.overlap.HL) / sd_random$niche.overlap.HL
    extinct <- (Para_properties$extinction.slope.HL - mean_random$extinction.slope.HL)  / sd_random$extinction.slope.HL
    nodf <- (Comm_properties$NODF - mean_random$NODF) / sd_random$NODF
    H2 <- (Comm_properties$H2 - mean_random$H2) / sd_random$H2
    Community <- Para_properties$Community
    
    newrow <- c(Community, niche, extinct, nodf, H2)
    
    z_scores1 <- rbind.data.frame(z_scores1, newrow)
    k <- k + 1
  }
} 

# Renaming columns in the dataset 

colnames(z_scores1) <- c("Community", "niche_overlap", "extinction_slope", "nodf","H2_specificity")

save(z_scores1, file =  "dataset_5degree1.RData")

######### Modularity ###########

# Function to calculate molecularity according to Beckett 2016
source("C:/Users/danid/OneDrive - Stanford/Otago/Lab Poulin/UORG/Codes/LPA_wb_plus.R")  

z_scores1 <- data.frame(stringsAsFactors = FALSE)
k <- 1
Lon <- -180
Lat <- -90

# Creating list with locations and host species
for (k in 1:length(occmatrix1)) {
  # print 
  print(k)
  
  occlist <- occmatrix1[[k]] %>%
    tidyr::pivot_longer(cols = -c(Latitude, Longitude),
                        names_to = "species",
                        values_to = "presence") %>%
    filter(presence >= 1) 
  
  occlist$presence <- NULL
  
  # Formatting host names to match phylogeny
  occlist$species <- gsub(" ","_", occlist$species)
  occlist$Lineage_Name <- "NA"
  
  # Filtering  parasite-host data per geographical location
  dados_com <- dados %>%
    filter(Longitude >= Lon & Longitude < Lon + 5) %>%
    filter(Latitude >= Lat & Latitude < Lat + 5)
  
  # Change value of Lon
  Lon <- Lon + 5
  
  # Change value of Lat when Lon >= 170
  if (Lon >= 175) {
    Lat <- Lat + 5
    Lon <- -180
  }
  
  dados_com <- dados_com[c(1,6)]
  
  dados_com <- rbind(dados_com, occlist[c(3,4)])
  
  # Check a condition
  if (length(dados_com$Lineage_Name[dados_com$Lineage_Name != "NA"]) < 20) {
    # Condition is true, move to the next iteration
    next
  } else {
    #it_sum <-it_sum + 1
    #}
    #}
    # Creating binary interaction matrix
    com <- table(dados_com$species, dados_com$Lineage_Name)
    #com[com>1] <- 1
    com <- as.matrix(unclass(com), nrow=nrow(com), ncol=ncol(com))
    
    # loading phylogeny and pruning to hosts in the interaction matrix
    index <- sample(length(random_trees2), 1)
    birds_tree <- as.phylo(random_trees2[[index]])
    
    # Getting consensus tree out of random 100 phylogenetic trees
    
    birds_tree <- drop.tip(birds_tree, birds_tree$tip.label[!birds_tree$tip.label%in%rownames(com)])
    
    # merge the tree and interaction matrix
    cleaned <- network_clean(com, birds_tree, 'full')
    com <- cleaned$Z                         # cleaned binary interaction matrix
    tree <- cleaned$tree                     # cleaned tree
    
    ## Running the model
    
    out <- network_est(Z = com, tree = tree, slices = 1000, model.type = 'full')
    
    ##### Calculating Network Metrics ######
    postocc <-  out$Z
    nas <- colnames(postocc) == "NA"
    postocc <- postocc[, colnames(postocc) != "NA"]
    
    # Network Metrics
    
    MOD = LPA_wb_plus(postocc) # find labels and weighted modularity using LPAwb+
    Comm_properties <- as.data.frame(MOD$modularity)
    Comm_properties$Modularity <- Comm_properties$`MOD$modularity`
    Comm_properties$`MOD$modularity` <- NULL
    Comm_properties$Community <- paste("Com", k, sep = " ")  
    Comm_properties <- as.data.frame(Comm_properties)
    
    # Calculating Z-Scores for parasite metrics
    
    null_model <- bipartite::nullmodel(postocc, method = "vaznul")
    
    ##### Use Loop to calculate it #####
    
    j <- 1
    
    NullScores <- data.frame(Modularity = numeric(),
                             Community = factor(),
                             stringsAsFactors = FALSE)
    
    for (j in 1:length(null_model)) {
      # Calculate group-level metrics for random values
      random_metrics1 <- LPA_wb_plus(null_model[[j]])
      random_metrics1 <- as.data.frame(random_metrics1$modularity)
      random_metrics1$Modularity <- random_metrics1$`MOD$modularity`
      random_metrics1$`MOD$modularity` <- NULL
      random_metrics1 <- t(random_metrics1)
      random_metrics1 <- data.frame(random_metrics1, stringsAsFactors = FALSE)
      
      # Adding values to dataframe
      NullScores <- rbind.data.frame(NullScores, random_metrics1)
      
      # summing J
      j <- j + 1
      
    }
    
    #Calculate mean and standard deviation for each community
    mean_random <- colMeans(NullScores)
    mean_random <- t(mean_random)
    mean_random <- as.data.frame(mean_random)
    
    sd_random <- apply(NullScores, 2, sd)
    sd_random <- t(sd_random)
    sd_random <- as.data.frame(sd_random)
    
    # Calculate Z-scores for each community
    Modularity <- (Comm_properties$Modularity - mean_random$random_metrics1) / sd_random$random_metrics1
    Community <- Comm_properties$Community
    
    newrow <- c(Community, Modularity)
    
    z_scores1 <- rbind.data.frame(z_scores1, newrow)
    k <- k + 1
  }
} 

# Renaming columns in the dataset 

colnames(z_scores1) <- c("Community", "Modularity")

save(z_scores1, file =  "Modularity_5degree.RData")


######## Bayesian Modeling ########

library(brms)
library(ggplot2)
load("dataset_5degree1.RData")
data <- drop_na(z_scores1)
load("Modularity_5degree.RData")

data <- inner_join(data, z_scores1, by = "Community")

# setting data class
sapply(data[,1:11], class) # checking data class 

data[,2:11] <- sapply(data[,2:11], as.numeric) # converting it to numeric

sapply(data[,1:11], class) # checking data class

data <- data[!is.nan(data$specialization_asymmetry), ]

# checking data distribution
hist(data$niche_overlap)
hist(data$nodf)
hist(data$extinction_slope)
hist(data$H2_specificity)
hist(data$Modularity)

###### Models ######

a <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(), axis.line = element_line(colour = "black"))

var1SA <- bf(niche_overlap ~ H2_specificity, family = gaussian()) #check

prior1SA <- get_prior(var1SA, data = data)
prior1SA

mod1SA <- brms::brm(niche_overlap ~ H2_specificity, data = data,
                    family = gaussian(),  chains = 8, 
                    cores = 4, future = TRUE,
                    iter = 10000,
                    prior = c(
                      prior(student_t(3, 0, 2.5), "Intercept"),
                      prior(student_t(3, 0, 2.5), "sigma")
                    ))

plot1SA <- plot(conditional_effects(mod1SA), theme = a)

summary(mod1SA)


var2SA <- bf(nodf ~ H2_specificity, family = gaussian()) #check

prior2SA <- get_prior(var2SA, data = data)
prior2SA

mod2SA <- brms::brm(nodf ~ H2_specificity, data = data,
                    family = gaussian(),  chains = 8, 
                    cores = 4, future = TRUE,
                    iter = 10000,
                    prior = c(
                      prior(student_t(3, 0.2, 2.5), "Intercept"),
                      prior(student_t(3, 0, 2.5), "sigma")
                    ))

plot2SA <- plot(conditional_effects(mod2SA), theme = a)

summary(mod2SA)


var3SA <- bf(extinction_slope ~ H2_specificity, family = gaussian()) #check

prior3SA <- get_prior(var3SA, data = data)
prior3SA

mod3SA <- brms::brm(extinction_slope ~ H2_specificity, data = data,
                    family = gaussian(),  chains = 8, 
                    cores = 4, future = TRUE,
                    iter = 10000,
                    prior = c(
                      prior(student_t(3, -1, 2.5), "Intercept"),
                      prior(student_t(3, 0, 2.5), "sigma")
                    ))

plot3SA <- plot(conditional_effects(mod3SA), theme = a)

summary(mod3SA)

var4SA <- bf(Modularity ~ H2_specificity, family = gaussian()) #check

prior4SA <- get_prior(var4SA, data = data)
prior4SA

mod4SA <- brms::brm(Modularity ~ H2_specificity, data = data,
                    family = gaussian(),  chains = 8, 
                    cores = 4, future = TRUE,
                    iter = 10000,
                    prior = c(
                      prior(student_t(3, 0.3, 2.5), "Intercept"),
                      prior(student_t(3, 0, 2.5), "sigma")
                    ))

plot4SA <- plot(conditional_effects(mod4SA), theme = a)

summary(mod4SA)

####### Checking models #######

pp_check(mod1SA)
plot(mod1SA)
pp_check(mod2SA)
plot(mod2SA)
pp_check(mod3SA)
plot(mod3SA)
pp_check(mod4SA)
plot(mod4SA)

######## Calculating regional metrics ########

load("occmatrix1.RData")
load("dados.RData")

load("dataset_5degree1.RData")
data <- drop_na(z_scores1)
load("Modularity_5degree.RData")

data <- inner_join(data, z_scores1, by = "Community")


k <- 1
Lon <- -180
Lat <- -90

com_cells <- data.frame(com = character(), 
                        Lat = character(), 
                        Lon = character(),
                        stringsAsFactors = FALSE)

# Creating list with locations and coordinates
for (k in 1:length(occmatrix1)) {
  # print 
  print(k)
  
  occlist <- occmatrix1[[k]] %>%
    tidyr::pivot_longer(cols = -c(Latitude, Longitude),
                        names_to = "species",
                        values_to = "presence") %>%
    filter(presence >= 1) 
  
  occlist$presence <- NULL
  
  # Formatting host names to match phylogeny
  occlist$species <- gsub(" ","_", occlist$species)
  occlist$Lineage_Name <- "NA"
  
  # Filtering  parasite-host data per geographical location
  dados_com <- dados %>%
    filter(Longitude >= Lon & Longitude < Lon + 5) %>%
    filter(Latitude >= Lat & Latitude < Lat + 5)
  
  # Create columns with Lat and Lon values
  
  cells_com <- paste("Com", k, sep = " ")
  cells_Lat <- paste(Lat, Lat + 5, sep = " ")
  cells_Lon <- paste(Lon, Lon + 5, sep = " ")
  
  cells <- cbind.data.frame(cells_com, cells_Lat, cells_Lon)
  
  com_cells <- rbind.data.frame(com_cells, cells)
  
  # Change value of Lon
  Lon <- Lon + 5
  
  # Change value of Lat when Lon >= 170
  if (Lon >= 175) {
    Lat <- Lat + 5
    Lon <- -180
  }
}

data2 <- inner_join(data, com_cells, by= c("Community" = "cells_com"))

##### extracting worldclim variables #####

data2 <- separate(data2, col = cells_Lat, into = c("Lat", "Lat2"), sep = " ")
data2 <- separate(data2, col = cells_Lon, into = c("Lon", "Lon2"), sep = " ")

data2[,12:15] <- sapply(data2[,12:15], as.numeric)

# Set the extent based on the latitude and longitude boundaries

i <- 1
ext <- list()
worldclim <- data.frame(bio1 = numeric(), 
                        bio2 = numeric(), 
                        bio3 = numeric(), 
                        bio4= numeric(),
                        bio5 = numeric(),
                        bio6 = numeric(),
                        bio7 = numeric(),
                        bio8 = numeric(),
                        bio9 = numeric(),
                        bio10 = numeric(), 
                        bio11 = numeric(), 
                        bio12 = numeric(),
                        bio13 = numeric(),
                        bio14 = numeric(),
                        bio15 = numeric(),
                        bio16 = numeric(),
                        bio17 = numeric(),
                        bio18 = numeric(),
                        bio19 = numeric(),
                        stringsAsFactors = FALSE)


for (i in 1:nrow(data2)) {
  print(i)
  ext[[i]] <- extent(data2[i,]$Lon, data2[i,]$Lon2, data2[i,]$Lat, data2[i,]$Lat2)
  
  ##### getting worldclim data for each quadrant #####
  bioclim.data <- raster::getData(name = "worldclim",
                                  var = c("bio"),
                                  res = 10) # Loading BioClim Data from WorldClim
  
  coords.range <- c(data2[i,]$Lon, data2[i,]$Lon2, data2[i,]$Lat, data2[i,]$Lat2) # Creating a vector with coords range
  values <- raster::extract(bioclim.data,ext[[i]]) # Extracting the BioClim values to our Sites coords
  
  # Calculate the average for each variable
  
  column_averages <- function(data) {   # creating function
    # Create a new dataframe to store column averages
    new_data <- data.frame(data)
    
    # Calculate column averages and add them to the new dataframe
    new_data <- colMeans(data, na.rm = TRUE)
    
    # Return the new dataframe
    return(new_data)
  }
  
  averages <- column_averages(values)
  averages <- t(averages)
  worldclim[i,] <- as.data.frame(averages)
  
}


##### Getting bird taxonomic and functional diversity for each quadrant ######

b <- 1

Tax_diversity <- data.frame(stringsAsFactors = FALSE)


for (b in 1:length(occmatrix1)) {
  # print 
  print(b)
  
  occlist <- occmatrix1[[b]] %>%
    tidyr::pivot_longer(cols = -c(Latitude, Longitude),
                        names_to = "species",
                        values_to = "presence") %>%
    filter(presence >= 1) 
  
  occlist$presence <- NULL
  
  # Formatting host names to match phylogeny
  occlist$species <- gsub(" ","_", occlist$species)
  occlist$Lineage_Name <- "NA"
  
  # Filtering  parasite-host data per geographical location
  dados_com <- dados %>%
    filter(Longitude >= Lon & Longitude < Lon + 5) %>%
    filter(Latitude >= Lat & Latitude < Lat + 5)
  
  # Create columns with Lat and Lon values
  
  cells_com <- paste("Com", b, sep = " ")
  cells_Lat <- paste(Lat, Lat + 5, sep = " ")
  cells_Lon <- paste(Lon, Lon + 5, sep = " ")
  
  # Change value of Lon
  Lon <- Lon + 5
  
  # Change value of Lat when Lon >= 170
  if (Lon >= 175) {
    Lat <- Lat + 5
    Lon <- -180
  }
  
  ##### Calculating Host Taxonomic Diversity #####
  
  library(hillR)
  
  Hill_Taxa <- hill_taxa(occmatrix1[[b]][, -c(1,2)], q = 1)
  
  # creating roll per cell
  tax_row <- c(mean(Hill_Taxa), cells_com)
  tax_row <- t(tax_row)
  tax_row <- as.data.frame(tax_row)
  
  #Adding values to dataset
  Tax_diversity <- bind_rows(Tax_diversity, tax_row)
  
}

colnames(Tax_diversity) <- c("Shannon", "Community")

data2 <- inner_join(data2, Tax_diversity, by = "Community")

save.image("UORG_03-07.RData")


##### Calculating Host Functional Diversity ##### 
library(dplyr)

## Function to remove species with zero occurrences
remove_zero_species <- function(occmatrix) {
  zero_species <- character(0)  # Initialize an empty vector to store zero species
  
  for (i in 2:ncol(occmatrix)) {
    if (all(occmatrix[, i] == 0)) {
      zero_species <- c(zero_species, colnames(occmatrix[i]))
    }
  }
  
  occmatrix_filtered <- occmatrix[, !colnames(occmatrix) %in% zero_species]
  
  return(occmatrix_filtered)
}

load("UORG_03-07.RData")

Fun_diversity <- data.frame(stringsAsFactors = FALSE)

for (c in 1:length(occmatrix1)) {
  # print 
  print(c)
  
  occlist <- occmatrix1[[c]] %>%
    tidyr::pivot_longer(cols = -c(Latitude, Longitude),
                        names_to = "species",
                        values_to = "presence") %>%
    filter(presence >= 1) 
  
  occlist$presence <- NULL
  
  # Formatting host names to match phylogeny
  occlist$species <- gsub(" ","_", occlist$species)
  occlist$Lineage_Name <- "NA"
  
  # Filtering  parasite-host data per geographical location
  dados_com <- dados %>%
    filter(Longitude >= Lon & Longitude < Lon + 5) %>%
    filter(Latitude >= Lat & Latitude < Lat + 5)
  
  # Create columns with Lat and Lon values
  
  cells_com <- paste("Com", c, sep = " ")
  cells_Lat <- paste(Lat, Lat + 5, sep = " ")
  cells_Lon <- paste(Lon, Lon + 5, sep = " ")
  
  # Change value of Lon
  Lon <- Lon + 5
  
  # Change value of Lat when Lon >= 170
  if (Lon >= 175) {
    Lat <- Lat + 5
    Lon <- -180
  }
  
  
  matrix <- remove_zero_species(occmatrix1[[c]])
  length <- length(matrix)
  
  if (length == 0) {
    next
  } else {
    
    traits <- readxl::read_excel("OpenTraits.xlsx")
    
    traits1 <- traits[,c(4,8,9,12,14,15,16)]
    traits1 <- filter(traits1, `Migration-2` != "NA")
    traits1 <- filter(traits1, Territoriality != "NA")
    traits1 <- filter(traits1, `Body_mass_log` != "NA")
    traits1 <- filter(traits1, Range_Size != "NA")
    traits1 <- filter(traits1, Habitat != "NA")
    
    traits1$IUCN_name <- gsub("_", " ", traits1$IUCN_name)
    
    comm <- matrix
    comm <- comm[-c(1,2)]
    
    
    # Get the values from the "Column_Names" column in df2
    column_names <- names(comm)
    
    # Filter columns in df1 based on matching column names
    traits1 <- traits1 %>%
      filter(IUCN_name %in% column_names)
    
    # Get the values from the "Column_Names" column in df2
    column_names <- traits1$IUCN_name
    
    # Filter columns in df1 based on matching column names
    comm <- comm %>%
      dplyr::select(all_of(column_names))
    
    length2 <- length(comm)
    if (length2 < 9) {
      next
    } else {
      
      # Changing character variables to numeric
      library(fastDummies)
      sapply(traits1, class)
      traits1$Body_mass_log <- as.numeric(traits1$Body_mass_log)
      traits1$Range_Size <- as.numeric(traits1$Range_Size)
      sapply(traits1[,4:7], n_distinct)
      
      # creating dummy variables
      traits1 <- fastDummies::dummy_columns(traits1, select_columns = c("Migration-2", "Territoriality", "Diet", "Habitat"),
                                            split = TRUE)
      # Adding species to rownames
      traits1 = traits1[!duplicated(traits1$IUCN_name),]
      traits1 <- traits1 %>%
        tibble::column_to_rownames(var = "IUCN_name") 
      
      ##### Calculating Host functional Diversity #####
      
      comm = comm[rowSums(comm[])>0,] # removing communities with 0 obs
      
      Hill_Func <- hill_func(comm = comm, traits = traits1, q = 1)
      
      # creating roll per cell
      fun_row <- c(mean(Hill_Func), cells_com)
      fun_row <- t(fun_row)
      fun_row <- as.data.frame(fun_row)
      
      #Adding values to dataset
      
      Fun_diversity <- bind_rows(Fun_diversity, fun_row)
    }
  }
}

save.image("UORG-04-07.RData")


####### Combining data frames #######

load("UORG_04-07.RData")

colnames(Fun_diversity) <- c("Func_Diversity", "Community")

data2 <- inner_join(data2, Fun_diversity, by = "Community")

data2 <- data2[,-c(7,8)]

data2 <- cbind.data.frame(data2, worldclim)

# Reducing Dimensions of WorldClim Data #

data3 <- data2
sapply(data3[,c(13:34)], class)
data3[,c(13:34)] <- lapply(data3[,c(13:34)], as.numeric)
data3[,c(13:34)] <- scale(data3[,c(13:34)])

# PCA temperature
Temperature <- prcomp(data3[,c(16:26)])
Temperature
summary(Temperature)
factoextra::fviz_pca_var(Temperature, col.var = "contrib")
varTemp <- get_pca_var(Temperature)
df1 <- data.frame(contrib = varTemp$coord[c(1:11)],
                  biodata = row.names(varTemp$coord))
ggplot(df1, aes(x = biodata, y = contrib)) + geom_col()


# PCA precipitation

Precipitation <- prcomp(data3[,c(27:34)])
Precipitation
summary(Precipitation)
factoextra::fviz_pca_var(Precipitation, col.var = "contrib")
varPrec <- get_pca_var(Precipitation)
df1 <- data.frame(contrib = varPrec$coord[c(1:8)],
                  biodata = row.names(varPrec$coord))
ggplot(df1, aes(x = biodata, y = contrib)) + geom_col()


WorlClim <- prcomp(data3[,c(16:34)])

data3 <- data3[,-c(16:34)]
data3 <- cbind.data.frame(data3, Temperature$x[,1], Precipitation$x[,1], WorlClim$x[,1])
colnames(data3)[16:18] <- c("Temp", "Prec", "WorlClim")

data3$Prec <- data3$Prec * -1

data3[,c(2:9)] <- lapply(data3[,c(2:9)], as.numeric)

##### Bayesian modeling #####

library(brms)

# checking for correlation

cor(data3$Shannon, data3$Func_Diversity)
cor(data3$Temp, data3$Prec)

# Nestedness 

mod_nest <-brm(nodf ~ Shannon + Temp + Prec, data = data3,
               family = gaussian(),  chains = 8, 
               cores = 4, future = TRUE,
               iter = 10000) 

summary(mod_nest)
parameters::p_value(mod_nest)

# Niche overlap

mod_nich <-brm(niche_overlap ~ Shannon + Temp + Prec, data = data3,
               family = gaussian(),  chains = 8, 
               cores = 4, future = TRUE,
               iter = 10000) 

summary(mod_nich)
parameters::p_value(mod_nich)

# Extinction Slope 

mod_ext <-brm(extinction_slope ~ Shannon + Temp + Prec, data = data3,
              family = gaussian(),  chains = 8, 
              cores = 4, future = TRUE,
              iter = 10000) 

summary(mod_ext)
parameters::p_value(mod_ext)

# Modularity 

mod_mod <-brm(Modularity ~ Shannon + Temp + Prec, data = data3,
              family = gaussian(),  chains = 8, 
              cores = 4, future = TRUE,
              iter = 10000) 

summary(mod_mod)
parameters::p_value(mod_mod)

###### Models Parameters #####

parameters::parameters(mod2SA) # nestedness (nodf)
summary(mod2SA)
parameters::parameters(mod1SA) # niche overlap
summary(mod1SA)
parameters::parameters(mod3SA) # extinction slope
summary(mod3SA)
parameters::parameters(mod4SA) # modularity
summary(mod4SA)


parameters::parameters(mod_nest) # nestedness (nodf)
summary(mod_nest)
parameters::parameters(mod_nich) # niche overlap
summary(mod_nich)
parameters::parameters(mod_ext) #extinction slope
summary(mod_ext)
parameters::parameters(mod_mod) # modularity
summary(mod_mod)


save.image("UORG_Bayesian.RData")

##### Model stability x extionction slope #####

vul_nest <- brm(extinction_slope ~ nodf, data = data3,
                family = gaussian(),  chains = 8, 
                cores = 4, future = TRUE,
                iter = 10000) 


summary(vul_nest)
parameters::parameters(vul_nest)

###### Final Plots #######

library(patchwork)
library(brms)

load("UORG_Bayesian.RData")

# Specificity 

plot1SA <- plot(conditional_effects(mod1SA), theme = a, points = TRUE)
plot2SA <- plot(conditional_effects(mod2SA), theme = a, points = TRUE)
plot3SA <- plot(conditional_effects(mod3SA), theme = a, points = TRUE)

plot1 <- plot1SA$H2_specificity + labs(x = "Community host-parasite specificity", y = "Parasite niche overlap (competition)", tag = "B") +
  scale_y_continuous(limits = c(-3, 3))
plot1

plot2 <- plot2SA$H2_specificity + labs(x = "Community host-parasite specificity", y = "Nestedness (community stability)", tag = "A") +
  scale_y_continuous(limits = c(-3, 3))
plot2

plot3 <- plot3SA$H2_specificity + labs(x = "Community host-parasite specificity", y = "Parasite extiction slope", tag = "C")
plot3

PlotVN <- plot(conditional_effects(vul_nest), theme = a, points = TRUE)
PlotVN_2 <- PlotVN$nodf + labs(x = "Nestedness (community stability)", y = "Parasite extiction slope", tag = "C") 

plot2 + plot1 + PlotVN_2

ggsave(filename = "Fig 1.pdf", dpi = 600,  width = 10, height = 4, device = "pdf")

# Climatic and host features

plot_nich <- plot(conditional_effects(mod_nich), theme = a, points = TRUE)

nich_div <- plot_nich$Shannon + labs(x = "Host diversity", y = "Parasite niche overlap (competition)", tag = "A") +
  scale_y_continuous(limits = c(-3, 3))
nich_temp <- plot_nich$Temp + labs(x = "Temperature (first PCA axis)", y = "Parasite niche overlap (competition)", tag = "B") +
  scale_y_continuous(limits = c(-3, 3))
niche_prec <- plot_nich$Prec + labs(x = "Precipitation (first PCA axis)", y = "Parasite niche overlap (competition)", tag = "C") +
  scale_y_continuous(limits = c(-3, 3))


(nich_div + nich_temp + niche_prec)

ggsave(filename = "Fig 2.pdf", dpi = 600,  width = 10, height = 4, device = "pdf")


##### LOO cross-validation via Pareto-smoothed importance sampling (PSIS) #####

# repeat code for each model

loo_results <- brms::loo(mod_nich) 

# Summary of LOO statistics
summary(loo_results)

# LOO log-likelihood estimate
loo_results$estimate

# Pareto k diagnostic values
loo_results$diagnostics$pareto_k

# mod_nich 1 k > 0.7

brms::loo_moment_match(mod_nich,loo_results)


####### MAP ########

###### Getting community geographical data #######

load("UORG.RData")
rm(occmatrix, MB, MB_data, random_trees2, mdist)
load("occmatrix1.RData")

com_data <- data.frame(stringsAsFactors = FALSE)
k <- 1
Lon <- -180
Lat <- -90

# Creating list with locations and host species
for (k in 1:length(occmatrix1)) {
  # print 
  print(k)
  
  occlist <- occmatrix1[[k]] %>%
    tidyr::pivot_longer(cols = -c(Latitude, Longitude),
                        names_to = "species",
                        values_to = "presence") %>%
    filter(presence >= 1) 
  
  occlist$presence <- NULL
  
  # Formatting host names to match phylogeny
  occlist$species <- gsub(" ","_", occlist$species)
  occlist$Lineage_Name <- "NA"
  
  # Filtering  parasite-host data per geographical location
  dados_com <- dados %>%
    filter(Longitude >= Lon & Longitude < Lon + 5) %>%
    filter(Latitude >= Lat & Latitude < Lat + 5)
  
  # Change value of Lon
  Lon <- Lon + 5
  
  # Change value of Lat when Lon >= 170
  if (Lon >= 175) {
    Lat <- Lat + 5
    Lon <- -180
  }
  
  dados_com <- dados_com[c(1,6:12)]
  
  # Check a condition
  if (length(dados_com$Lineage_Name[dados_com$Lineage_Name != "NA"]) < 20) {
    # Condition is true, move to the next iteration
    next
  } else {
    
    # creating roll per cell
    com_Lon <- Lon
    com_Lat <- Lat
    com_cont <- dados_com[1,3]
    com_size <- length(unique(dados_com$species)) * length(unique(dados_com$Lineage_Name))
    com_links <- sum(occmatrix1[[k]])
    com_number <- paste("Com", k, sep = " ")  
    
    row <- c(com_Lon, com_Lat, com_cont, com_number, com_size, com_links)
    
    #Adding values to dataset
    
    com_data <- rbind.data.frame(com_data, row)
  }
}  

rm(occmatrix1)

colnames(com_data) <- c("Longitude", "Latitude", "Continent", "Com_ID", "Com_size", "Links")

###### Creating map #######

install.packages("maps")
library(ggplot2)
library(maps)
library(graphics)

sapply(com_data, class)
com_data$Continent <- as.factor(com_data$Continent)
com_data[,c(1,2,5,6)] <- sapply(com_data[,c(1,2,5,6)], as.numeric)
com_data$Links <- abs(com_data$Links)

write.csv2(com_data, file = "community data.csv")

newmap <- rworldmap::getMap(resolution = "high")
cex <- log(com_data$Com_size)
cex <- cex/4

pdf(file="Fig 1.pdf", width = 8, height = 6)
plot(newmap, asp = 1, col = "mistyrose1", xlim = c(-165, 165), ylim = c(-70, 70))
graphics::points(com_data$Longitude, com_data$Latitude, cex = cex, col= com_data$Continent, pch = 21, bg = bg_col[com_data$Continent],  lwd = 2)
dev.off()

