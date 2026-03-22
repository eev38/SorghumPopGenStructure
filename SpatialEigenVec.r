library(tidyverse)
library(geosphere)
library(ggformula)
library(vegan)

load("SNPs.RDA")
load("clim_lang2.RDA")


# unique locations 
Acc_EnvGBS_temp$Lon <- round(Acc_EnvGBS_temp$Lon, 1)
Acc_EnvGBS_temp$Lat <- round(Acc_EnvGBS_temp$Lat, 1)
#unique coordinates w/ accessions
Acc_EnvGBS <- Acc_EnvGBS_temp[!duplicated(Acc_EnvGBS_temp[c("Lon","Lat")]),]
Geo_A <- Acc_EnvGBS[,c(3,2)]
Geo_B <- Acc_EnvGBS[,c(3,2)]
#create distance matrix
dist_mat <- distm(Geo_A,Geo_B,fun = distVincentyEllipsoid)
Geo_Acc <- Acc_EnvGBS$Accessionsfull
rownames(dist_mat) <- Geo_Acc #assign name to each pairwise combo
colnames(dist_mat) <- Geo_Acc
# run pcnm on dist matrix
space1j <- pcnm(dist_mat)
#make matrix with PCNM for RDA/varpart
spac.matj <- space1j$vectors[, which(space1j$values > 0)] #chose only positive eigenvalues

########

# run varpart sequentially on pcnm output 
#track adjusted r-squared for each set
    # find max adjust R-squared
        #output sequential columns with the max adjusted r-squared

#dfs from RDA files loaded above
Lan <- clim_lang2[,c(1:10)]
SNPs <- SNPs[,-1]

y = colnames(spac.matj)[0]   
TempList <- c(y)
varout = list()
for (col in colnames(spac.matj)) {   
    TempList[[length(TempList) + 1]] <- col # list of vectors with each column name added sequentially
        print(TempList)
        var <- varpart(Y=SNPs, 
                spac.matj[,c(TempList)],
                #need second explanatory variable to run - add language
                   Lan[,c("Atla", "Mand", "Nubi", "Pidg", "Afro", "Sign", "Tane", "Gumu", "Unat", "Khoe")],  
               data= spac.matj)
         tmp <- var$part[[2]]
        AdjR2 <- tmp[1,3]
        varout[[length(varout)+1]]= c(col,AdjR2) 
}

#extract max adjusted r-squared
varout_df <- do.call(rbind.data.frame, varout)
names(varout_df)[1] <- "LastCol"
names(varout_df)[2] <- "AdjR2"
#col 159 / #pcnm158 max adjusted r-squared
maxpcnm <- varout_df %>% 
   slice_max(AdjR2, n = 1)
