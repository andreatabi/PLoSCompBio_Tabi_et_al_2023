rm(list=ls())
library(dplyr)

# load data
df <- read.csv("~/data.csv")  

# filter data 
df <-  df %>% filter(SR>4 & sampling_effort>1)   

#########################################################################################################
### Create result database
#########################################################################################################

# number of sites
sites <- length(unique(df$site_code))
species <- length(unique(df$species_name))

# number of data points
id <- sort(unique(df$ID))
res <- data.frame()

for(i in 1:length(id)){

      data <- df[df$ID %in% id[i],]
      sp <- unique(data$species_name)
      site <- unique(data$site_code)
      se <- unique(data$sampling_effort)
      
      S <- list()
      for(j in 1:length(sp) ){
        temp <- data %>% dplyr::filter(species_name == as.character(sp[j]) ) %>% dplyr::select(species_name, biomass, total ) 
        temp <- na.omit(temp)
        S[[ sp[j] ]] <- rep((temp$biomass / temp$total), temp$total) 
      }

      res[ i, "site_code"] <- site
      res[ i, "year"] <- unique(data$year)
      res[ i, "ID"] <- id[i]
      res[ i, "latitude"] <- unique(data$lat)
      res[ i, "longitude"] <- unique(data$long)
      res[ i, "site_name"] <- unique(data$site_name )
      res[ i, "ecoregion"] <- unique(data$ecoregion )
      res[ i, "realm"] <- unique(data$realm)
      res[ i, "area"] <- unique(data$area)
      res[ i, "SumTSAfull"] <- unique(data$SumTSAfull)
      res[ i, "Buffer10Km"] <- unique(data$Buffer10Km)
      res[ i, "HumanPopDensity"] <- unique(data$HumanPopDensity)
      res[ i, "IUCN_CAT"] <- unique(data$IUCN_CAT)
      res[ i, "IUCN_Ia"] <- unique(data$IUCN_Ia)
      res[ i, "STATUS_YR"] <- unique(data$STATUS_YR)
      res[ i, "MPA"] <- unique(data$MPA)
      res[ i, "sampling_effort"] <- se
      res[ i, "sr"] <- length(sp)

      # Spearman's correlations
      AB_local <- unlist(lapply( S, length)) # abundance 
      BS_local <- unlist(lapply( S, mean)) # body size 
      BM_local <- unlist(lapply( S, sum)) # biomass 
      
      res[i, "Reg_BS_BM_pvalue"] <- summary(lm(log10(BM_local) ~ log10(BS_local)))$coef[2,4] 
      res[i, "Reg_BS_BM_coef"] <- summary(lm(log10(BM_local) ~ log10(BS_local)))$coef[2,1] 

      par(mfrow=c(1,2))
      plot(log10(BS_local), log10(BM_local), main=id[i] )
      abline(lm(log10(BM_local) ~ log10(BS_local)))

      TL <- unique(data[, c("species_name", "trophic_level")])
      plot(TL$trophic_level, log10(BS_local), main=id[i] )
      abline(lm(log10(BS_local) ~ TL$trophic_level ))

      res[ i, "rho_bs_tl"] <-  cor.test( log10(BS_local) , TL$trophic_level, method="spearman",exact=F)$est
      res[ i, "rho_bs_tl_pvalue"] <-  cor.test( log10(BS_local) , TL$trophic_level, method="spearman",exact=F)$p.value
      res[ i, "cwm_size"] <- weighted.mean(BS_local, AB_local, na.rm=T)

      print(paste0("#", i, ", site = ", id[i])) 

}

# save results -------------------------------------------------------------------------------------------------------
write.csv(res, file="~/results.csv", row.names=FALSE)



