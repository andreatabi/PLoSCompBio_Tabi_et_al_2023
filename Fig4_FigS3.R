rm(list=ls())
library(dplyr)
library(ggplot2) 
library(ggpattern)
library(ggpubr)
db <- read.csv('~/results.csv') 

k <- as.numeric(db$Reg_BS_BM_coef)
tsa <- as.numeric(db$SumTSAfull)
hd <- as.numeric(db$HumanPopDensity)
pro <- as.numeric(db$IUCN_Ia)
coral <- as.numeric(db$Buffer10Km)

################################################################################################## 
### Transform continuous to binary variables using median
################################################################################################## 

bk <- k
bk[k > median(k,na.rm=T)] <- 1
bk[k <= median(k,na.rm=T)] <- 0

bhd <- hd
bhd[hd > median(hd,na.rm=T)] <- 1
bhd[hd <= median(hd,na.rm=T)] <- 0

btsa <- tsa
btsa[tsa > median(tsa,na.rm=T)] <- 1
btsa[tsa <= median(tsa,na.rm=T)] <- 0

bcoral <- coral
bcoral[coral >median(coral,na.rm=T)] <- 1
bcoral[coral <=median(coral,na.rm=T)] <- 0

bpro <- pro


df <- cbind.data.frame( k=db$Reg_BS_BM_coef, 
                        mpa = ifelse(bpro==0, "Disturbed", "Protected"),
                        coral = ifelse(bbelow==0, "Rocky reef", "Coral reef"),
                        human = ifelse(bst==0, "Low human density", "High human density"),
                        tsa = ifelse(bsst==0, "Low TSA", "High TSA") )

################################################################################################## 
### Figure 4
################################################################################################## 

plotA <- ggplot(df) + 
  geom_violin( aes(y=k , x=mpa, fill=mpa),  color="white"  ) + 
  ggtitle("")+
  geom_boxplot( aes(y=k , x=mpa ), outlier.shape=1, width=0.2 ) +
  xlab("")+ ylab(expression(paste("Structure (", k[c]^e, ")")) )+
  scale_fill_manual(values=alpha(c("#FFD500","#005BBB"),0.4))+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_hline(yintercept=median(df$k), linetype="dashed", color="gray40", size=1 )+
  theme(legend.position = "none", axis.title = element_text(size=20),
        axis.text.x = element_text(size=20), axis.text.y = element_text(size=18),
        strip.text = element_text(size = 20), 
        strip.background = element_rect(fill="white", color="white"))

plotB <- ggplot(df) + 
  facet_wrap(~coral, nrow=1)+
  geom_violin( aes(y=k , x=mpa, fill=mpa),  color="white"  ) +
  geom_boxplot( aes(y=k , x=mpa ), outlier.shape = 1, width=0.2 ) +
  xlab("")+ ylab(expression(paste("Structure (", k[c]^e, ")")) )+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_hline(yintercept=median(df$k), linetype="dashed", color="gray40", size=1 )+
  scale_fill_manual(values=alpha(c("#FFD500","#005BBB"),0.4))+
  theme(legend.position = "none", axis.title = element_text(size=20),
        axis.text.x = element_text(size=20), axis.text.y = element_text(size=18),
        strip.text = element_text(size = 20), 
        strip.background = element_rect(fill="white", color="white"))

plotC <- ggplot(df, aes(y=k , x=coral, fill=coral)) + 
        geom_violin( color="white" ) +
        ggtitle("")+
        geom_boxplot( aes(y=k , x=coral ), outlier.shape=1, fill="white",  width=0.2 ) +
        xlab("")+ ylab("")+
        scale_fill_manual(values=alpha(c("#FFD500","#005BBB"),0.4))+
        theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
        geom_hline(yintercept=median(df$k), linetype="dashed", color="gray40", size=1 )+
        theme(legend.position = "none", axis.title = element_text(size=20),
              axis.text.x = element_text(size=20), axis.text.y = element_text(size=18),
              strip.text = element_text(size = 20), 
              strip.background = element_rect(fill="white", color="white"))

plotD <- ggplot(df %>% dplyr::select(k, mpa, tsa) %>% na.omit, aes(y=k , x=mpa, fill=mpa)) + 
  facet_wrap(~tsa, nrow=1)+
  geom_violin( color="white" ) +
  geom_boxplot( aes(y=k , x=mpa ), outlier.shape=1, fill="white",  width=0.2 ) +
  xlab("")+ ylab(expression(paste("Structure (", k[c]^e, ")")) )+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_hline(yintercept=median(df$k), linetype="dashed", color="gray40", size=1 )+
  scale_fill_manual(values=alpha(c("#FFD500","#005BBB"),0.4))+
  theme(legend.position = "none", axis.title = element_text(size=20),
        axis.text.x = element_text(size=20), axis.text.y = element_text(size=18),
        strip.text = element_text(size = 20), 
        strip.background = element_rect(fill="white", color="white"))

plotE <- ggplot(df %>% dplyr::select(k, mpa, human) %>% na.omit, aes(y=k , x=mpa, fill=mpa)) + 
  facet_wrap(~human, nrow=1)+
  geom_violin(color="white" ) +
  geom_boxplot( aes(y=k , x=mpa ), outlier.shape=1,  fill="white",  width=0.2 ) +
  xlab("")+ ylab(expression(paste("Structure (", k[c]^e, ")")) )+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_hline(yintercept=median(df$k), linetype="dashed", color="gray40", size=1 )+
  scale_fill_manual(values=alpha(c("#FFD500","#005BBB"),0.4))+
  theme(legend.position = "none", axis.title = element_text(size=20),
        axis.text.x = element_text(size=20), axis.text.y = element_text(size=18),
        strip.text = element_text(size = 20), 
        strip.background = element_rect(fill="white", color="white"))

ggarrange(plotA, plotB, plotC, plotD, ncol = 2, nrow = 3,
          widths = c(3, 9),
          font.label = list(size = 22),
          align = "hv",
          labels = "AUTO")

plot1 <- ggarrange(plotA, plotC,  ncol = 2, nrow = 1,
          widths = c(3, 3),
          font.label = list(size = 22),
          align = "hv",
          labels = "AUTO")

ggarrange(plot1, plotB, plotD, plotE, ncol = 1, nrow = 4,
          font.label = list(size = 22),
          align = "hv",
          labels = c("", "C", "D", "E" ))


ggsave("~/Fig4.tiff", width=9, height=14, units = "in")

################################################################################################## 
### Figure S3
################################################################################################## 

db$fact <- ifelse(db$IUCN_Ia == 1, "Protected", "Disturbed")

cwm <- ggplot(db)+ 
  geom_histogram(aes(x=log10(cwm_size), group=fact, fill=fact), color="white", alpha=0.8)+
  scale_fill_manual(values = c("#FFD500" ,"#005BBB"))+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylab("Communities")+
  xlab(bquote(log[10]~"("~CWM[BM]~") distribution" )) +
  theme(legend.title=element_blank(), legend.position=c(0.19, 0.88),
        legend.text = element_text(size=12),
        axis.title = element_text(size=16),
        axis.text.x = element_text(size=16), axis.text.y = element_text(size=16),
        strip.text = element_text(size = 16),
        strip.background = element_rect(fill="white", color="white"))

bs_tl <- ggplot(db)+ 
  geom_histogram(aes(x=rho_bs_tl, group=fact, fill=fact), color="white", alpha=0.8)+
  scale_fill_manual(values = c("#FFD500" ,"#005BBB"))+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylab("")+
  xlab(bquote(rho[BM:TL]~ "distribution")) +
  theme(legend.title=element_blank(), legend.position="none",
        legend.text = element_text(size=16),
        axis.title = element_text(size=16),
        axis.text.x = element_text(size=16), axis.text.y = element_text(size=16),
        strip.text = element_text(size = 16),
        strip.background = element_rect(fill="white", color="white"))

ggarrange( cwm, bs_tl, ncol = 2, nrow = 1,
           font.label = list(size = 20),
           align = "hv",
           labels = "AUTO")

ggsave("~/FigS3.tiff", width=8, height=4)

