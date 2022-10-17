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

devtools::load_all()
pceg.obs<-pceg(knee_pain_obs,2,TRUE,TRUE)
pceg.obs$result
expected_count_calculator(knee_pain_obs,pceg.obs,zip=FALSE)

zipceg.obs<-zipceg(knee_pain_obs,"nlm",variable_time=TRUE)
zipceg.obs$result
expected_count_calculator(knee_pain_obs,zipceg.obs)

zipceg.iter(knee_pain_obs,"nlm",iter_total=10,violin=TRUE,variable_time=TRUE)
zipceg.iter(knee_pain_obs,"nlm",iter_total=10,hist=TRUE,variable_time=TRUE)
zipceg.iter(knee_pain_obs,"nlm",iter_total=10,scatter=TRUE,variable_time=TRUE)

hdi_gamma_extractor(knee_pain_obs,pceg.obs,zip=FALSE)
hdi_beta_extractor(knee_pain_obs,pceg.obs,level_rel_final=-1,zip=FALSE)
hdi_beta_extractor(knee_pain_obs,pceg.obs,level_rel_final=-2,zip=FALSE)
hdi_beta_extractor(knee_pain_obs,pceg.obs,level_rel_final=-3,zip=FALSE)

hdi_gamma_extractor(knee_pain_obs,zipceg.obs,zip=TRUE)
hdi_beta_extractor(knee_pain_obs,zipceg.obs,level_rel_final=-1,zip=TRUE)




