###########################################################################
###                           Samarth Mathur, PhD                       ###
###                       The Ohio State University                     ###
###                                                                     ###
###     Date Created: 07/15/22                  Last Modified: 07/15/22 ###
###########################################################################
###########################################################################
###                   DoS.R         		                           ###
###########################################################################

library(plyr)
library(dplyr)
library(ggplot2)
library(scales)
library(MASS)
setwd("/fs/ess/scratch/PAS1533/smathur/Scat_Popgen/tables/")

# P for polymorphic
# S for synonymous
# F for fixed (used here as same as D for divergent)
# and R for replacement.  (used here as same for N for non-synonymous)

# Case 1: KLDR as ingroup, STER as outgroup
# Case 2: STER as ingroup, KLDR as outgroup
df1.main <- read.table("KLDR.ig_STER.og.2x2table.txt",header = T) # 17257


# PR FR PS FS = Pn Dn Ps Ds

colnames(df1.main)[2:5] <- c("Pn","Dn","Ps","Ds")




# Remove genes with Pn + Dn < 4
df1.filtered <- df1.main[-which((df1.main$Pn+df1.main$Dn)<4),] #  N = 3974
dos1 <- df1.filtered$Dn/(df1.filtered$Dn+df1.filtered$Ds) - df1.filtered$Pn/(df1.filtered$Pn+df1.filtered$Ps)

# Remove Nas
df1.dos <- cbind(df1.filtered,dos1)
df1.dos <- df1.dos[-which(is.na(df1.dos$dos)),] # N = 3337
write.csv(df1.dos,"DoS.KLDRig.csv",quote = F,row.names = F)

# Choosing genes based on quantile

h1 <- hist(df1.dos$dos, breaks = 100, plot = FALSE)

# upper and lower 10%
quantile(df1.dos$dos, c(0.10,0.9))
#5%        95% 
#  -0.5  0.5


cuts1 <- cut(h1$breaks, c(-Inf,-0.524,0.5,Inf))
 
plot(h1, col=c(col=adjustcolor("coral1",alpha.f = 0.5),adjustcolor("dodgerblue4",alpha.f = 0.5),adjustcolor("darkolivegreen4",alpha.f = 0.5),
              col=adjustcolor("dodgerblue4",alpha.f = 0.5),adjustcolor("coral1",alpha.f = 0.5))[cuts1], 
     border=c(col=adjustcolor("coral1",alpha.f = 0.7),adjustcolor("dodgerblue4",alpha.f = 0.7),adjustcolor("darkolivegreen4",alpha.f = 0.7),
              col=adjustcolor("dodgerblue4",alpha.f = 0.7),adjustcolor("coral1",alpha.f = 0.7))[cuts1],
     main = "KLDR ingroup, STER outgroup",
     xlab=expression(paste("Direction of selection (DoS)")),
     ylab = "No. of genes (Total = 3,337)")
abline(v=mean(df1.dos$dos), col="dodgerblue4", lty=2, lwd=1.5) # -0.01257679

