# pcegr
Poisson Chain Event Graphs and Zero Inflated Poisson Chain Event Graphs.

This package contains the functionality to fit both Poisson Chain Event Graphs (PCEGs) and Zero-Inflated Poisson Chain Event Graphs (ZIPCEGs) using the pceg() function. To be specific, the pceg() function creates a StagedTree S3 object, which can be zero-inflated or not, while the ChainEventGraph() constructor function creates a ChainEventGraph S3 object from a StagedTree.

This package contains a data set, knee_pain_obs, which can be used to fit either a PCEG or ZIPCEG model. To fit a PCEG model, run the following code:

library(pcegr)  
ptree<-pceg(knee_pain_obs,2,TRUE,TRUE, gamma_alpha = 0.25,gamma_beta = 0.1)

This will fit a Poisson staged tree model to the data set, creating a StagedTree object. The outputs can be seen as a whole by calling ptree, but summary() contains a shortened output and $result is the stage structure.

summary(ptree)
ptree$result

We can plot this model:

plot(ptree)

and convert to a ChainEventGraph object, which can also be plotted:

pceg1<-ChainEventGraph(ptree)
plot(pceg1)

We can check the posterior estimates of the parameters using:
ptree$posterior.expectation

There are various functions available to perform analysis of the model. The quantile_band() function produces (shifted) quantile band plots which can check the fit of the model. 

quantile_band(knee_pain_obs,ptree)

The marginal_effect() function allows the user to investigate the impact of the variables on each other. For example, to see how each covariate affects the response variable:

marginal_effect(knee_pain_obs,ptree)

The vif() function identifies multicollinearity between the covariates in the model. To see how collinear the third covariate is with the first two:

vif(ptree,3)

We can also fit a ZIPCEG model to the data set. 

set.seed(13579)  
ziptree1<-zipceg(knee_pain_obs,"nlm",iter=10000,stoch_imputation = TRUE, gamma_alpha = a, gamma_beta = b)  
summary(ziptree1)
plot(ziptree1)  

We can once again check the parameter estimates.  

ziptree1$posterior.expectation

and form a CEG that can also be plotted.

plot(ChainEventGraph(ziptree1))

We check the quantile band plot again:

quantile_band(knee_pain_obs,ziptree1,zip=TRUE)

We can also compare the fit between the two models numerically using the Chi-square:  

chi_sq_calculator(knee_pain_obs,pmod,zip=FALSE,dec_place=2)$chi_sq  
chi_sq_calculator(knee_pain_obs,zipmod1,zip=TRUE,dec_place=2)$chi_sq  



