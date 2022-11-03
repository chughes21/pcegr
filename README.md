# pcegr
Poisson Chain Event Graphs and Zero Inflated Poisson Chain Event Graphs.

This package contains the functionality to fit both Poisson Chain Event Graphs (PCEGs) and Zero-Inflated Poisson Chain Event Graphs (ZIPCEGs) using the [pceg()] function. The package can also convert the outputs of this function to the S4 classes of Stratified.staged.tree and Stratified.ceg.model from the ceg package, which enables plotting. 

This package contains a data set, knee_pain_obs, which can be used to fit a PCEG and ZIPCEG model. To fit a PCEG model, run the following code:

library(pcegr)  
pmod<-pceg(knee_pain_obs,2,TRUE,TRUE, gamma_alpha = 0.25,gamma_beta = 0.1)

This will fit a Poisson staged tree model to the data set. The outputs can be seen as a whole by calling pmod, or simply the stage structure through:

pmod$result

However, usually the better way to observe a Poisson staged tree is through a plot. First, the output should be converted to the S4 class Stratified.staged.tree from the ceg package, and then plotted.

ptree<-staged.tree.creator(knee_pain_obs,pmod)  
plot(ptree)  

We can also view this as a Stratified PCEG model, an extention of the Stratified.ceg.model S4 class from the ceg package. First, the staged tree model should be converted to a PCEG, and then plotted.

pceg1<-sceg(ptree)  
plot(pceg1)  

We can check the posterior estimates of the parameters using either the original pmod, or the S4 class.
total_value_extractor(knee_pain_obs,mod,zip=FALSE)
ptree@posterior.distribution

We can also fit a ZIPCEG model to the data set. First, we fit the model, then convert it to the S4 classes, and then can plot them both.

set.seed(13579)  
zipmod1<-zipceg(knee_pain_obs,"Gibbs",iter=10000,variable_time = TRUE, stoch_imputation = TRUE, gamma_alpha = a, gamma_beta = b, beta_c = c, beta_d = d)  
ziptree1<-staged.tree.creator(knee_pain_obs,zipmod1,zip=TRUE)  
zipceg1<-sceg(ziptree1)  
plot(ziptree1)  
plot(zipceg1)  

We can once again check the parameter estimates.  
total_value_extractor(knee_pain_obs,zipmod1,zip=TRUE)  
ziptree1@posterior.distribution  

We can also compare the fit between the two models:  
chi_sq_calculator(knee_pain_obs,pmod,zip=FALSE,dec_place=2)$chi_sq  
chi_sq_calculator(knee_pain_obs,zipmod1,zip=TRUE,dec_place=2)$chi_sq  

