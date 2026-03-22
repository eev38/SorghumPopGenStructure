library(tidyverse)
library(eulerr)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(vegan)
library(raster)
library(sp)
library(maps)
library(dplyr)

# Load geographic distance eigen vectors

load("spac.matj.RDA")

spac.matj2_temp <- tibble::rownames_to_column(as.data.frame(spac.matj), "Accessions")  

spac.matj2 <- spac.matj2_temp[,1:159] 

# Load Language Family data

FamilyName_Coords  <-  read.csv('~/AllFamilyData_2025.csv')  %>% 
    dplyr::select(Accessions, Family_Name, Longitude_lan, Latitude_lan)  %>% 
    filter(Accessions %in% spac.matj2$Accessions)

EnvGBS_Acc <- read.csv('Acc_EnvGBS.csv')   %>% 
   rename(Accessions = Accessionsfull)  %>% 
    filter(Accessions %in% FamilyName_Coords$Accessions)

load('~/SNPTable_2025_PLINK.RData')

# Load SNP data

SNPs <- SNPTable_2025_PLINK   %>% 
    filter(Accessions %in% EnvGBS_Acc$Accessions)

# Load CHELSA Climate Variables    

#searching files that end with "tif" 
fileslist  <- list.files("~/CHELSAbioclim_V.2.1_Jan10.26", full.names=T, pattern="tif$")

## Moisture

Annual_Precipitation <- raster(paste0(fileslist[12]))
EnvGBS_Acc$Annual_Precipitation <- raster::extract(Annual_Precipitation, EnvGBS_Acc[,c("Lon","Lat")])

Precipitation_Driest_Month <- raster(paste0(fileslist[14]))
EnvGBS_Acc$Precipitation_Driest_Month <- raster::extract(Precipitation_Driest_Month, EnvGBS_Acc[,c("Lon","Lat")])
p <- (EnvGBS_Acc$Precipitation_Driest_Month + 1)
EnvGBS_Acc$Precipitation_Driest_Month_log <- log(p)

Precipitation_Seasonality <- raster(paste0(fileslist[15]))
EnvGBS_Acc$Precipitation_Seasonality <- raster::extract(Precipitation_Seasonality, EnvGBS_Acc[,c("Lon","Lat")])

Precipitation_Warmest_Quarter <- raster(paste0(fileslist[18]))
EnvGBS_Acc$Precipitation_Warmest_Quarter <- raster::extract(Precipitation_Warmest_Quarter, EnvGBS_Acc[,c("Lon","Lat")])

## Temperature

Temperature_Seasonality <- raster(paste0(fileslist[4]))
EnvGBS_Acc$Temperature_Seasonality <- raster::extract(Temperature_Seasonality, EnvGBS_Acc[,c("Lon","Lat")])

Max_Temp_Warmest_Month <- raster(paste0(fileslist[5]))
EnvGBS_Acc$Max_Temp_Warmest_Month <- raster::extract(Max_Temp_Warmest_Month, EnvGBS_Acc[,c("Lon","Lat")])

Min_Temp_Coldest_Month <- raster(paste0(fileslist[6]))
EnvGBS_Acc$Min_Temp_Coldest_Month <- raster::extract(Min_Temp_Coldest_Month, EnvGBS_Acc[,c("Lon","Lat")])

Mean_Temp_Driest_Quarter <- raster(paste0(fileslist[9]))
EnvGBS_Acc$Mean_Temp_Driest_Quarter <- raster::extract(Mean_Temp_Driest_Quarter, EnvGBS_Acc[,c("Lon","Lat")])

EnvGBS_Acc$country <- map.where("world", EnvGBS_Acc$Lon, EnvGBS_Acc$Lat)


EnvGBS_Acc3 <- EnvGBS_Acc  %>% 
    dplyr::select(-Precipitation_Driest_Quarter, -Precipitation_Driest_Month)  %>% 
   rename(
              PreW =Precipitation_Warmest_Quarter, 
              AP = Annual_Precipitation, 
              PDML = Precipitation_Driest_Month_log,
             PS = Precipitation_Seasonality,
              MTDQ = Mean_Temp_Driest_Quarter,
              TS=Temperature_Seasonality,
               MTWM=Max_Temp_Warmest_Month,
              MTCM=Min_Temp_Coldest_Month
              )


# Prepare language family data for varpart
Fam <- FamilyName_Coords  %>%  dplyr::select(Accessions, Family_Name)

#convert family df to wide format (binary presence/absence for each language)
Fam_untidy <- Fam  %>% 
  add_column(Present=1)  %>% 
    pivot_wider(names_from=Family_Name, values_from=Present)  %>% 
     replace(is.na(.), 0)  %>% 
      rename(Afro ='Afro-Asiatic', Atla= 'Atlantic-Congo', 
             Asta ='Austroasiatic',
             Blnm ='Blue Nile Mao', Book ='Bookkeeping', Cent ='Central Sudanic', Daju ='Dajuic', Dizo ='Dizoid',
             Dogo ='Dogon', Dravi ='Dravidian', East ='Eastern Jebel', Fura ='Furan', Gumu ='Gumuz', Heib='Heibanic',
             Hurr='Hurro-Urartian', Indo ='Indo-European', Kadu ='Kadugli-Krongo', Kart ='Kartvelian',
             Khoe='Khoe-Kwadi', Koma ='Koman', Kres='Kresh-Aja', Kru ='Kru', Kuli ='Kuliak', Kxa ='Kxa', Maba='Maban',
             Mand ='Mande', Narr='Narrow Talodi', Nilo ='Nilotic', Nubi ='Nubian', Nyim ='Nyimang',
             Pidg ='Pidgin', Saha ='Saharan', Sign='Sign Language', Sino ='Sino-Tibetan', Song='Songhay',
             Spee='Speech Register', Surm='Surmic', Tane ='Ta-Ne-Omotic',
             Tama='Tamaic', Tuu ='Tuu', Unat='Unattested', Uncl='Unclassifiable')  


###merging dfs
clim_lang_t1 <- merge(Fam_untidy, EnvGBS_Acc3, by= "Accessions")  

clim_lang <- merge(clim_lang_t1, spac.matj2, by= "Accessions")   %>% 
    dplyr::select(-country)

SNPs_filter <- merge(clim_lang, SNPs, by = 'Accessions') 
###

### Dataframes for varpart

# Response Data (Y) - SNPs
SNPs_filter2 <- SNPs_filter[,-c(1:211)]             

# Explanatory variables (X) - Language, Climate, Geographic distance eigen values
clim_lang2 <- clim_lang  %>% 
    dplyr::select(-Accessions)


### Run varpart

var <- varpart(Y=SNPs_filter2, 
                clim_lang2[,c('AP', 'PDML','PS','PreW')], # Moisture
                clim_lang2[,c('TS','MTWM', 'MTCM','MTDQ')], # Temp
                clim_lang2[,c(spac.matj_varIDs)], #geographic space
               clim_lang2[,c('Afro','Asta','Atla','Blnm','Book','Cent','Daju','Dizo','Dogo','Dravi','East',
                             'Fura','Gumu','Heib','Hurr','Indo','Kadu','Kart','Khoe','Koma','Kres','Kru',
                             'Kuli','Kxa','Maba','Mand','Narr','Nilo','Nubi','Nyim','Pidg','Saha','Sign',
                             'Sino','Song','Spee','Surm','Tama','Tane','Tuu', 'Unat' )],  # language removed 'Uncl'
               data= clim_lang2)
