rm(list=ls())
dev.off()

library(reshape2)
library(dplyr)
library(ggplot2)
library(trophic)
library(ggpubr)
library(NetIndices)

# Function to simulate community structure with harvest ####################################
# Simulate #################################################################################
set.seed(19850827)
sp <- 35                                  
frac <- 0.4                               
TE.max <- rep(c(0.1, 0.3, 0.5, 0.7), 1000)    
# run the simulations ######################################################################

res <- list()
for(i in 1:length(TE.max)){
  
        print(i)
        C <- sp^(-0.65)                           
        FM <- trophic::niche(sp, C)               
        TL <- TrophInd(FM)[,1]                   
        E <- rnorm(sp, 0, 1)                     
        Q <- 10^3 
        M0 <- 1
        mean.M <- M0* Q^(TL-1+E) 
        sd.M <- 0.1 * mean.M                     
        PPMRs <- FM * matrix(rep(mean.M, sp), nrow=sp, byrow=T) / matrix(rep(mean.M, sp), nrow=sp, byrow=F)
        mPPMR <- median((PPMRs[-which(PPMRs==0)] ))

        ss <- mean.M^-0.03
        TE <- TE.max[i] * ss/max(ss)    
        
        # size-selective harvest 
        probs.M <- mean.M/sum(mean.M)              
        fished_sp <- which(mean.M %in% sample(mean.M, floor(sp*frac), prob=probs.M))  # sample larger species 
        remove <- runif(sp, 0.3, 1)
        
        harvest <- lapply(fished_sp, function(x) {
          dist <-  rnorm( 5000, mean = mean.M[[x]], sd = sd.M[[x]] )
          probs <- pnorm(dist, mean = mean.M[[x]], sd = sd.M[[x]])
          harvest.dist <- dist[-which(dist %in% sample(dist, floor(remove[[x]]*5000), prob=probs))]
          new.mean.M <- mean(harvest.dist)
          sd.new.M <- sd(harvest.dist)
          remove.bm <- (sum(dist)  - sum(harvest.dist)) / sum(dist) 
          return(data.frame(sp=x, mean.M=mean.M[[x]], new.mean.M=new.mean.M, remove.bm=remove.bm, remove.ab=remove[[x]]))
        })
        
        
        df0 <- suppressMessages({ melt(harvest) })
        df <- tidyr::spread(df0, variable, value)
        harvest.mean.M <-  mean.M   
        harvest.mean.M[fished_sp] <- df$new.mean.M
        
        
        ki <- NULL
          for(j in 1:sp){
            pp <- subset(PPMRs[,j], PPMRs[,j]>0 )
            pp[which(pp<=1 )] <- 0.9
            te <- rep(TE[j], length(pp))
            te[which( pp < 1)] <- 1-te[which( pp < 1)]
            
            if(length(pp)==0){ 
              w <- 1
              pp <- Q 
            } else { 
              w <- rep(1,length(pp)) / length(pp)
            }
            ki[j] <-  0.25 +  sum(w*(log10(te)/log10(pp)))
          }
          
          B <- (mean.M)^ki    
          if(any(B==Inf))  B[which(B==Inf)] <- max(B[-which(B==Inf)])
          k <- summary(lm( log10(B+0.001) ~ log10(mean.M)))$coef[2,1]

          harvest.B <- B
          harvest.B[fished_sp] <- (1-df$remove.bm) * B[fished_sp]
          harvest.k <- summary(lm(log10(harvest.B+0.001) ~ log10(harvest.mean.M)))$coef[2,1]

          co_bs_tl <- summary(lm(log10(mean.M) ~ TL ))$coef[2,1]
          rho_bs_tl <- cor.test(log10(mean.M) , TL, method = "spearman", exact=F)$est
          
          res[[i]] <- data.frame( n = i,
                                  sp=sp, 
                                  TE.max=TE.max[i], 
                                  k=k, 
                                  harvest.k=harvest.k,
                                  co_bs_tl = co_bs_tl,
                                  rho_bs_tl = rho_bs_tl, 
                                  medPPMR=mPPMR, 
                                  maxTL=max(TL) )
    
  
}

# put database together ####################################################################
d <- data.frame(do.call( "rbind", res))
d <- na.omit(d)
dd <- melt(d, id=c( "n", "sp", "TE.max", "co_bs_tl", "rho_bs_tl", "medPPMR",  "maxTL"))
dd$fact <- ifelse(dd$variable == "k", "Protected", "Disturbed")
te_string <- as_labeller(c( `0.1` = "TE[max] == 0.1", `0.3` = "TE[max] == 0.3",
                            `0.5` = "TE[max] == 0.5", `0.7` = "TE[max] == 0.7"), default = label_parsed)

############################################################################################
### Figure 1 ###############################################################################
############################################################################################
head(dd)

sim <- ggplot(dd)+ 
  theme_bw()+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0.25, fill = "#005BBB", alpha=0.1 ) +
  geom_violin(data=dd, aes(x=factor(fact), y=value, fill=fact), color="white"  ) +
  ylim(-1,1)+
  geom_boxplot(data=dd, aes(x=factor(fact), y=value, fill=fact), outlier.shape=16, outlier.size = 0.8, width=0.1, 
               outlier.colour = alpha("grey", 0.9) )+
  scale_fill_manual(values=alpha(c( "#FFD500","#005BBB"),0.4))+
  geom_hline(aes(yintercept=0.25), colour = "#005BBB" , linetype="dashed")+
  facet_wrap( ~ TE.max, nrow=1, labeller = te_string)+
  xlab(" ") + 
  ylab(expression(paste("Structure (", k[c], ")")) )+
  theme(legend.position = "none", axis.title = element_text(size=16),
        axis.text.x = element_text(size=14), axis.text.y = element_text(size=16),
        strip.text = element_text(size = 16), 
        strip.background = element_rect(fill="white", color="white"))
sim

mean.M <-  10^5
sd.M <- mean.M * 0.1
remove <- 0.7
dist <-  rnorm( 5000, mean = mean.M, sd = sd.M )
probs <- pnorm(dist, mean = mean.M, sd = sd.M)
harvest.dist <- dist[-which(dist %in% sample(dist, floor(remove*5000), prob=probs))]

d1 <- data.frame(bs=dist, treat="Protected")
d2 <- data.frame(bs=harvest.dist, treat="Harvested")

dist_plot <- ggplot()+ 
  geom_histogram(data=d1, aes(x=log10(bs), fill=treat), color="white", alpha=0.6)+
  geom_histogram(data=d2, aes(x=log10(bs), fill=treat), color="white", alpha=0.9)+
  scale_fill_manual(values = c("#FFD500" ,"#005BBB"))+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylab("Abundance")+
  xlab(bquote(log[10]~"(body size)" )) +
  theme(legend.title=element_blank(), legend.position=c(0.25, 0.8),
        legend.text = element_text(size=12),
        axis.title = element_text(size=16),
        axis.text.x = element_text(size=16), axis.text.y = element_text(size=16),
        strip.text = element_text(size = 16),
        strip.background = element_rect(fill="white", color="white"))

dist_plot

ppmr_plot <- ggplot(d)+ 
  geom_histogram(aes(x=log10(medPPMR)), color="white", fill="grey30", alpha=0.6)+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylab("Communities")+
  xlab(bquote(log[10]~"("~PPMR~")" )) +
  theme(legend.position = "none", axis.title = element_text(size=16),
        axis.text.x = element_text(size=16), axis.text.y = element_text(size=16),
        strip.text = element_text(size = 16),
        strip.background = element_rect(fill="white", color="white"))
ppmr_plot

BS_TL_plot <- ggplot(d)+ 
  geom_histogram(aes(x=rho_bs_tl ), color="white", fill="grey30", alpha=0.6)+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylab("Communities")+
  xlab(bquote(rho[BM:TL] )) +
  theme(legend.position = "none", axis.title = element_text(size=16),
        axis.text.x = element_text(size=16), axis.text.y = element_text(size=16),
        strip.text = element_text(size = 16),
        strip.background = element_rect(fill="white", color="white"))
BS_TL_plot

plot1 <- ggarrange( BS_TL_plot, dist_plot, ppmr_plot, ncol = 3, nrow = 1,
          font.label = list(size = 20),
          align = "hv",
          labels = "AUTO")

ggarrange(plot1, sim, ncol = 1, nrow = 2,
          font.label = list(size = 20),
          align = "hv", heights=c(0.9,1),
          labels = c("", "D"))

ggsave("~/Fig1.tiff", width=12, height=7)

############################################################################################
### Figure S2 ##############################################################################
############################################################################################

sp <- 35
C <- sp^(-0.65)                          # connectance 
FM <- trophic::niche(sp, C)              # feeding matrix 
TL <- TrophInd(FM)[,1]                   # trophic levels
E <- rnorm(sp, 0, 1)                     # noise to BS
Q <- 10^3 
M0 <- 1
mean.M <- M0* Q^(TL-1+E) 

PPMRs <- FM * matrix(rep(mean.M, sp), nrow=sp, byrow=T) / matrix(rep(mean.M, sp), nrow=sp, byrow=F)
PPMRs <- apply(PPMRs, 2, function(x) mean(x[-which(x==0)], na.rm=T)  )
PPMRs[is.na(PPMRs)] <- Q

ss <- mean.M^-0.03
TE <- 0.7 * ss/max(ss)                     # original transfer efficiency dependent on body size (Woodson et al 2018)

TE_plot <- ggplot()+ 
  geom_point(aes(x=TL, y=TE ), color="black", fill="grey30", alpha=0.6)+
  geom_smooth(aes(x=TL, y=TE ), method = "lm", color="black", fill="blue", alpha=0.3)+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylab("Trophic efficiency")+
  xlab("Trophic level")+
  theme(legend.position = "none", axis.title = element_text(size=16),
        axis.text.x = element_text(size=16), axis.text.y = element_text(size=16),
        strip.text = element_text(size = 16),
        strip.background = element_rect(fill="white", color="white"))
TE_plot

ppmr <- ggplot()+ 
  geom_point(aes(y=log10(unlist(PPMRs+1)), x=TL ), color="black", fill="grey80", alpha=0.6)+
  geom_smooth(aes(y=log10(unlist(PPMRs+1)), x=TL ), method = "lm", color="black", fill="blue", alpha=0.3)+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylab("log10(PPMR)")+
  xlab("Trophic level") +
  theme(legend.position = "none", axis.title = element_text(size=16),
        axis.text.x = element_text(size=16), axis.text.y = element_text(size=16),
        strip.text = element_text(size = 16),
        strip.background = element_rect(fill="white", color="white"))
ppmr

ggarrange( TE_plot, ppmr, ncol = 2, nrow = 1,
                    font.label = list(size = 20),
                    align = "hv",
                    labels = "AUTO")

ggsave("~/FigS2.tiff", width=8, height=4)


