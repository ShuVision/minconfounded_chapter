-------------------------------------------------------------------------------
README file corresponding to submitted code for

	"Are You Normal? The Problem of Confounding in Hierarchical Linear Models"	

by 

	Adam Loy and Heike Hofmann
-------------------------------------------------------------------------------
FILE						DESCRIPTION
-------------------------------------------------------------------------------
original_radon.csv			The original radon data described by Gelman and
							Pardoe (2006) that is the motivating example
							for the paper. 
								
full_rank_radon.csv			A subset of the above data set that contains only
							those counties with full rank Z matrices.

srrs2.csv					Data file used to construct the map of radon
							in Minnesota.

radon_intro.R				R script to generate figures 1-4.

heatmap_cartoon.R			R script used to generate figure 5.

small_sim.R					R script used to perform the small simulation
							study presented in the introduction.

simulate_method.R			R script used to perform the simulation study
							for the rotated random effects.

cpp_functions.R				A few C++ functions to help speed up matrix
							algebra in R.

normality_functions.R		Functions used to test normality in the 
							simulation studies.

resid_functions.R			Functions to calculate the minimally confounded
							residuals and other residuals for a hierarchical
							linear model.

utility_function.R			Helper functions for lme4 for the simulation studies

radon_revisited.R			R script used to generate table 4 and figure 11