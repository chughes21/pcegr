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