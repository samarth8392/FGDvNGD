###########################################################################
###                         Samarth Mathur, PhD                         ###
###                       The Ohio State University                     ###
###                                                                     ###
###     Date Created: 04/02/22                  Last Modified: 04/05/22 ###
###########################################################################
###########################################################################
###                   ldNE.R       		                    			###
###########################################################################

#### PREREQUISITES #####
library(plyr)
library(dplyr)
library(DataCombine)
library(ggplot2)
setwd("/Users/batcomputer/Documents/Postdoc_research/WGR/Results/PostMay21/ldne")

pops <- c("BPNP","CCRO","CEBO","JENN","KBPP","KLDR","MOSQ","PRDF","ROME","SPVY","SSSP","WLRD")

## Plot LDne using only unrelated individuals

ldne <- read.table("LDNe.intergenic.geneDs.unrelInds.txt",header=T)

# Recalculate for negative sq
rows <- which(is.na(ldne$Sq))
for (i in rows)
{
  ldne$NeLd[i] <- (1/3)/(2*ldne$R2_prime[i])
}
  

ldne.mean <- ddply(ldne,c("Pop"),summarise,Size=mean(Size),mean_Nehat=mean(NeLd),sd_Nehat=sd(NeLd),
                     popmeanR2=mean(meanR2),popsdR2=sd(meanR2),
                     meanR2prime=mean(R2_prime),sdR2prime=sd(R2_prime))

write.table(ldne.mean,"MeanLDNe.intergenic.geneDs.unrelInds.txt",quote = F,row.names = F)

ggplot(ldne.mean) +
  geom_bar( aes(x=Pop, y=mean_Nehat, fill=Pop), stat="identity", alpha=0.7) +
  geom_errorbar( aes(x=Pop, ymin=mean_Nehat-sd_Nehat, ymax=mean_Nehat+sd_Nehat), width=0.4, colour="black", alpha=0.9)+
  theme_classic(base_size = 12)+ 
  theme(legend.position = "none")+
  labs(y=expression(paste(N[e(LD)])), x ="Location")+
  theme(axis.text.y =element_text(size=9))
