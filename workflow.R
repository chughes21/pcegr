usethis::use_git()
usethis::use_github(private=TRUE)
usethis::use_r("em_algorithm_functions")
devtools::document()
usethis::use_r("expected_counts_function")
usethis::use_r("gibbs_sampler_functons")
usethis::use_r("interval_creators")
usethis::use_r("imports")
usethis::use_r("numerical_optimisation_functions")
usethis::use_r("other_functions")
usethis::use_r("pceg_function")
usethis::use_r("value_extractor_function")
usethis::use_r("zip_functions")
devtools::document()
usethis::use_data_raw("knee_pain")
usethis::use_data_raw("knee_pain_obs")
devtools::document()
devtools::load_all()
usethis::use_r("prior_function")
usethis::use_r("covariate_calculator_function")
usethis::use_r("uniform_time_functions")
usethis::use_vignette("pcegr")
usethis::use_data_raw("knee_pain_param")
usethis::use_data_raw("knee_pain_age")
usethis::use_data_raw("knee_pain_age_obs")
usethis::use_r("ceg_package_updated_functions")
usethis::use_r("ceg_package_new_functions")
usethis::use_r("ceg_package_plotting_functions")
usethis::use_r("class_definitions")
usethis::use_r("deleted_function")
usethis::use_r("christmas_tree_functions")
usethis::use_r("marginal_effect_functions")
usethis::use_r("vif")
usethis::use_r("kl_divergence")
usethis::use_r("situation_stage_monitor")
usethis::use_r("binary_resizer")

devtools::load_all()
devtools::check()


pceg.obs<-pceg(knee_pain_obs,2,TRUE,TRUE)
pceg.obs$result
expected_count_calculator(knee_pain_obs,pceg.obs,zip=FALSE)

for(i in 0:3){
  print(value_extractor(knee_pain_obs,pceg.obs,level_rel_final = i-3,zip=FALSE))
}

zipceg.obs<-zipceg(knee_pain_obs,"nlm",variable_time=TRUE)
zipceg.obs$result
expected_count_calculator(knee_pain_obs,zipceg.obs)

zipceg.iter(knee_pain_obs,"nlm",iter_total=10,variable_time=TRUE)
zipceg.iter(knee_pain_obs,"nlm",iter_total=10,line=TRUE,variable_time=TRUE)
zipceg.iter(knee_pain_obs,"nlm",iter_total=10,hist=TRUE,variable_time=TRUE)
zipceg.iter(knee_pain_obs,"nlm",iter_total=10,scatter=TRUE,variable_time=TRUE)

zipceg.iter(knee_pain_obs,"nlm",iter_total=10,plot_rates=FALSE,plot_probs=TRUE,variable_time=TRUE)
zipceg.iter(knee_pain_obs,"nlm",iter_total=10,plot_rates=FALSE,plot_probs=TRUE,line=TRUE,variable_time=TRUE)
zipceg.iter(knee_pain_obs,"nlm",iter_total=10,plot_rates=FALSE,plot_probs=TRUE,hist=TRUE,variable_time=TRUE)
zipceg.iter(knee_pain_obs,"nlm",iter_total=10,plot_rates=FALSE,plot_probs=TRUE,scatter=TRUE,variable_time=TRUE)

hdi_gamma_extractor(knee_pain_obs,pceg.obs,zip=FALSE)
hdi_beta_extractor(knee_pain_obs,pceg.obs,level_rel_final=-1,zip=FALSE)
hdi_beta_extractor(knee_pain_obs,pceg.obs,level_rel_final=-2,zip=FALSE)
hdi_beta_extractor(knee_pain_obs,pceg.obs,level_rel_final=-3,zip=FALSE)

hdi_gamma_extractor(knee_pain_obs,zipceg.obs,zip=TRUE)
hdi_beta_extractor(knee_pain_obs,zipceg.obs,level_rel_final=-1,zip=TRUE)

df<-knee_pain_obs
ind<-which(df$y>0)
bin.resp<-factor(rep("No",10000),levels=c("Yes","No"))
bin.resp[ind]<-"Yes"
df<-data.frame(df,resp=bin.resp)
df<-df[,-c(4:5)]
summary(df)

mod1<-pceg(df,2)
mod1$result
tree1<-set(df)
plot(tree1)
tree2<-event.tree.creator(df,FALSE,FALSE)
plot(tree2)

stagedtree1<-staged.tree.creator(df,mod1,FALSE,FALSE)
plot(stagedtree1)

ceg1<-sceg(stagedtree1)
plot(ceg1)

mod2<-pceg(knee_pain_obs,2,TRUE,TRUE,gamma_alpha=0.25,gamma_beta=0.1)
mod2$result

tree2<-event.tree.creator(knee_pain_obs,TRUE,TRUE)
plot(tree2)

#need to fix numb output in pceg

stagedtree2<-staged.tree.creator(knee_pain_obs,mod2)
plot(stagedtree2)

ceg2<-sceg(stagedtree2)
plot(ceg2)

