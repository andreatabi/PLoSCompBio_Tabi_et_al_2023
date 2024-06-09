rm(list=ls())
dev.off()

# load data
df <- read.csv("~/data.csv")  

# Species richness and abundance distributions

# comm > 4 & only se > 1
df1 <-  df %>% filter(SR>4 & sampling_effort>1)   
TA1 <- df1 %>% dplyr::select(ID, Total_abundance) %>% distinct()
SR1 <- df1 %>% dplyr::select(ID, SR) %>% distinct()
  
hist( log(TA1$Total_abundance), col="#005BBB", xlab = "log total abundance", 
     ylab = "# of communities", breaks=30, main=NA, border="cornflowerblue" )
abline(v=median(log(TA1$Total_abundance), na.rm=T), col="#FFD500", lty=2, lwd=1.5)

hist( log(SR1$SR), col="#005BBB", xlab = "log richness", 
      ylab = "# of communities", breaks=30, main=NA, border="cornflowerblue" )
abline(v=median(log(SR1$SR), na.rm=T), col="#FFD500", lty=2, lwd=1.5)

# comm > 4 & SE = 1 
df2 <-  df %>% filter(SR>4  & sampling_effort==1 )   
TA2 <- df2 %>% dplyr::select(ID, Total_abundance) %>% distinct()
SR2 <- df2 %>% dplyr::select(ID, SR) %>% distinct()

hist( log(TA2$Total_abundance), col="#005BBB", xlab = "log total abundance", 
      ylab = "# of communities", breaks=30, main=NA, border="cornflowerblue" )
abline(v=median(log(TA2$Total_abundance), na.rm=T), col="#FFD500", lty=2, lwd=1.5)

hist( log(SR2$SR), col="#005BBB", xlab = "log richness", 
      ylab = "# of communities", breaks=30, main=NA, border="cornflowerblue" )
abline(v=median(log(SR2$SR), na.rm=T), col="#FFD500", lty=2, lwd=1.5)

# comm > 4 & only SE > 2 
df3 <-  df %>% filter(SR>4 & sampling_effort > 2)   
TA3 <- df3 %>% dplyr::select(ID, Total_abundance) %>% distinct()
SR3 <- df3 %>% dplyr::select(ID, SR) %>% distinct()

hist( log(TA3$Total_abundance), col="#005BBB", xlab = "log total abundance", 
      ylab = "# of communities", breaks=30, main=NA, border="cornflowerblue" )
abline(v=median(log(TA3$Total_abundance), na.rm=T), col="#FFD500", lty=2, lwd=1.5)

# comm > 4 & all
df4 <-  df %>% filter(SR>4)   
TA4 <- df4 %>% dplyr::select(ID, Total_abundance) %>% distinct()
SR4 <- df4 %>% dplyr::select(ID, SR) %>% distinct()

hist( log(TA4$Total_abundance), col="#005BBB", xlab = "log total abundance", 
      ylab = "# of communities", breaks=30, main=NA, border="cornflowerblue" )
abline(v=median(log(TA4$Total_abundance), na.rm=T), col="#FFD500", lty=2, lwd=1.5)

# comm > 4 & only se > 3
df5 <-  df %>% filter(SR>4 & sampling_effort > 3)   
TA5 <- df5 %>% dplyr::select(ID, Total_abundance) %>% distinct()
SR5 <- df5 %>% dplyr::select(ID, SR) %>% distinct()

hist( log(TA5$Total_abundance), col="#005BBB", xlab = "log total abundance", 
      ylab = "# of communities", breaks=30, main=NA, border="cornflowerblue" )
abline(v=median(log(TA5$Total_abundance), na.rm=T), col="#FFD500", lty=2, lwd=1.5)

# comm > 4 & only se > 4
df6 <-  df %>% filter(SR>4 & sampling_effort > 4)   
TA6 <- df6 %>% dplyr::select(ID, Total_abundance) %>% distinct()
SR6 <- df6 %>% dplyr::select(ID, SR) %>% distinct()

hist( log(TA6$Total_abundance), col="#005BBB", xlab = "log total abundance", 
      ylab = "# of communities", breaks=30, main=NA, border="cornflowerblue" )
abline(v=median(log(TA6$Total_abundance), na.rm=T), col="#FFD500", lty=2, lwd=1.5)

######## KS tests

# 1. Compare dist SE > 1 and SE == 1 
ks.test(TA1$Total_abundance, TA2$Total_abundance)  # not same distribution
ks.test(SR1$SR, SR2$SR)  # not same distribution 

# 2. Compare dist SE = 1 and all 
ks.test(TA2$Total_abundance, TA4$Total_abundance)  #  same distribution
ks.test(SR2$SR, SR4$SR)  #  same distribution

# 3. Compare dist SE > 1 and all 
ks.test(TA1$Total_abundance, TA4$Total_abundance)  # not same distribution
ks.test(SR1$SR, SR4$SR)  # not same distribution

# 4. Compare dist SE > 2 and all 
ks.test(TA3$Total_abundance, TA4$Total_abundance)  # not same distribution
ks.test(SR3$SR, SR4$SR)  # not same distribution

# 5. Compare dist SE > 1 and SE > 2 
ks.test(TA1$Total_abundance, TA3$Total_abundance)  # somewhat similar
ks.test(SR1$SR, SR3$SR)  # not same distribution

# 6. Compare dist SE > 2 and SE > 3 
ks.test(TA3$Total_abundance, TA5$Total_abundance)  # same distribution
ks.test(SR3$SR, SR5$SR)  # same distribution

# Compare dist SE > 3 and SE > 4 
ks.test(TA5$Total_abundance, TA6$Total_abundance)  # same distribution
ks.test(SR5$SR, SR6$SR)  # same distribution

