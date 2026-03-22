library(tidyverse)


# Read in Language group data/coordinates ffrom glottolog
geo_data <- read.csv('languages_and_dialects_geo.csv')

# filter NAs in latitude and longitude columns, then select those col and put in df
lang_latlon <- geo_data  %>%
     filter_at(vars(latitude, longitude), all_vars(!is.na(.)))  %>% 
#filter to match coords from landraces 
    filter(between(longitude, -17, 90), 
                      between(latitude, -31, 42))  %>% 
     select(longitude, latitude)

lang_matrix <- as.matrix(lang_latlon)


# Calculate distance between languages and accessions

dist_mat <- distm(Acc_matrix,lang_matrix,fun = distVincentyEllipsoid)
# sorghum accession data
Acc_vec <- TBLS1_LL$Accessions 
Lang_vec <- lang_latlon_NoNA$name
# attach row/col vector names to matrix 
rownames(dist_mat) <- Acc_vec
colnames(dist_mat) <- Lang_vec
Distmat_new <- tibble::rownames_to_column(
        as.data.frame(dist_mat), "Accessions")

# Find shortest distance language group for each accession  
#create a matrix with 3 cols and nrow(dist_mat) rows
Mat_langAcc = matrix(nrow = nrow(dist_mat) , ncol = 3 )
# loop through each row/accession as described above
 for (i in 1:nrow(dist_mat)){
    # for each row, find the lan(column) with the min. value (distance)
    closest <- which (dist_mat[i,] == min(dist_mat[i,])) 
    # col 1 should have acc rownames
     Mat_langAcc[i,1] <- rownames(dist_mat)[i]
    # col 2 should have col name of closest value
     Mat_langAcc[i,2] <- colnames(dist_mat)[closest[1]]
    # col 3 should have closest value
     Mat_langAcc[i,3] <- dist_mat[i,closest[1]]
     }

CloseLan_df <- as.data.frame(Mat_langAcc)  %>% 
    rename('Accessions'= 'V1', 'Language' = 'V2', 'Distance'='V3')

# Extract language family info
fam_name_df <- read.csv('languages.csv')
fam_name_df2 <- fam_name_df %>% rename('family_id' = 'ID', 'Family_Name' = 'Name')

# load in languages with family id 
Lan_Codes <- read.csv('languoid.csv')

Fam_Names_ID <-  left_join(Lan_Codes, fam_name_df2) %>% 
    select(family_id, Family_Name, glottocode, Macroarea)

# filtered geodata used  to construct matrix
geodata2 <- geo_data  %>%
     filter_at(vars(latitude, longitude), all_vars(!is.na(.)))  %>% 
#filter to match coords from landraces 
    filter(between(longitude, -17, 90), 
                      between(latitude, -31, 42))    

Family_ID2 <- dplyr::inner_join(geodata2, Fam_Names_ID, by='glottocode')  %>% 
    select(glottocode,family_id,Family_Name, name, level, macroarea,latitude, longitude) %>% 
    rename('Language' = name, 'Latitude_lan'= 'latitude', 'Longitude_lan'= 'longitude')                      

Family_Dist <- dplyr::inner_join(CloseLan_df, Family_ID2, by='Language')    

#join family DIST to accession coordinates
AllFamilyData_2025 <- dplyr::inner_join(Family_Dist, TBLS1_LL, by='Accessions')  %>% 
 rename('Latitude_Acc'='Lat', 'Longitude_Acc'='Lon')
