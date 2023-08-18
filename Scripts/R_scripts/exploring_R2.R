################################################################
################################################################
# exploring R^2 in NGD/FGD analyses in response to Reviewer 1
################################################################
################################################################

################################
#	re-run analsyes presented in paper
################################
# read in data frame
d <- read.csv("FGDvNGD.pop.csv")

# run models of various measures of FGD vs. pi
m_Ndel_vs_pi <- lm(meanNdel ~ meanPi, data=d)
m_Ldrift_vs_pi <- lm(meanLdrift ~ meanPi, data=d)
m_Linbreed_vs_pi <- lm(meanHet ~ meanPi, data=d)
m_Lrealized_vs_pi <- lm(meanHom ~ meanPi, data=d)
m_piAdapt_vs_pi <- lm(pi_u10 ~ meanPi, data=d)

# run models of various measures of FGD vs. NeLD
m_Ndel_vs_NeLD <- lm(meanNdel ~ meanNeLD, data=d)
m_Ldrift_vs_NeLD <- lm(meanLdrift ~ meanNeLD, data=d)
m_Linbreed_vs_NeLD <- lm(meanHet ~ meanNeLD, data=d)
m_Lrealized_vs_NeLD <- lm(meanHom ~ meanNeLD, data=d)
m_piAdapt_vs_NeLD <- lm(pi_u10 ~ meanNeLD, data=d)

################################
#	run analyses with standardized predictors
################################

# run models of various measures of FGD vs. pi
m_Ndel_vs_pi_scl <- lm(scale(meanNdel) ~ scale(meanPi), data=d)
m_Ldrift_vs_pi_scl <- lm(scale(meanLdrift) ~ scale(meanPi), data=d)
m_Linbreed_vs_pi_scl <- lm(scale(meanHet) ~ scale(meanPi), data=d)
m_Lrealized_vs_pi_scl <- lm(scale(meanHom) ~ scale(meanPi), data=d)
m_piAdapt_vs_pi_scl <- lm(scale(pi_u10) ~ scale(meanPi), data=d)

# run models of various measures of FGD vs. NeLD
m_Ndel_vs_NeLD_scl <- lm(scale(meanNdel) ~ scale(meanNeLD), data=d)
m_Ldrift_vs_NeLD_scl <- lm(scale(meanLdrift) ~ scale(meanNeLD), data=d)
m_Linbreed_vs_NeLD_scl <- lm(scale(meanHet) ~ scale(meanNeLD), data=d)
m_Lrealized_vs_NeLD_scl <- lm(scale(meanHom) ~ scale(meanNeLD), data=d)
m_piAdapt_vs_NeLD_scl <- lm(scale(pi_u10) ~ scale(meanNeLD), data=d)

################################
#	compare standardized effect sizes
################################

#pull out the effect sizes of the standardized predictors
#	for pi and then for NeLD
pi_effect_sizes <- list("Ndel" = summary(m_Ndel_vs_pi_scl)$coefficients[2,1],
						"Ldrift" = summary(m_Ldrift_vs_pi_scl)$coefficients[2,1],
						"Linbreed" = summary(m_Linbreed_vs_pi_scl)$coefficients[2,1],
						"Lrealized" = summary(m_Lrealized_vs_pi_scl)$coefficients[2,1],
						"piAdapt" = summary(m_piAdapt_vs_pi_scl)$coefficients[2,1])


NeLD_effect_sizes <- list("Ndel" = summary(m_Ndel_vs_NeLD_scl)$coefficients[2,1],
							"Ldrift" = summary(m_Ldrift_vs_NeLD_scl)$coefficients[2,1],
							"Linbreed" = summary(m_Linbreed_vs_NeLD_scl)$coefficients[2,1],
							"Lrealized" = summary(m_Lrealized_vs_NeLD_scl)$coefficients[2,1],
							"piAdapt" = summary(m_piAdapt_vs_NeLD_scl)$coefficients[2,1])

# plot the comparison
#	note that we can take the absolute value because the 
#	sign of the standardized coefficients are the same for pi and for NeLD
#	for each predictor (i.e., both positive or both negative)
pdf(file="standardized_effect_sizes.pdf")
	plot(abs(unlist(pi_effect_sizes)),abs(unlist(NeLD_effect_sizes)),
			xlab="standardized effect sizes with pi as a predictor",
			ylab="standardized effect sizes with NeLD as a predictor",
			type="n",
			xlim=c(0.4,1.1),
			ylim=c(0.4,0.8),
			main="Comparing standardized effect sizes of pi and NeLD\non measures of functional genetic diversity")
		text(abs(unlist(pi_effect_sizes)),abs(unlist(NeLD_effect_sizes)),
				labels=c("Ndel","Ldrift","Linbreed","Lrealized","piAdapt"))
		abline(0,1,col="red",lty=2)
	legend(x="bottomright",lty=2,col="red",legend="x=y")
dev.off()
						
						
						
						
						
						

