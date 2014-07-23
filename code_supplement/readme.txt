-------------------------------------------------------------------------------
README file corresponding to submitted code for

	"Are You Normal? The Problem of Confounding in Hierarchical Linear Models"	

by 

	Adam Loy and Heike Hofmann
-------------------------------------------------------------------------------
FILE					DESCRIPTION
-------------------------------------------------------------------------------
original_radon.csv		The original radon data described by Gelman and
						          Pardoe (2006) that is the motivating example
						          for the paper. 
								
full_rank_radon.csv		A subset of the above data set that contains only
						          those counties with full rank Z matrices.

srrs2.csv				      Data file used to construct the map of radon
						          in Minnesota.

figures1_4.R			R script to generate figures 1-4.

figure_5.R				R script to generate figure 5.

figure_6.R				R script to generate figure 6.

figure_7.R				R script to generate figure 7.

figures8_9.R			  R script to generate figures 8 and 9.

figures10_11.R			R script to generate figures 10 and 11.

small_sim.R				  R script to perform the small simulation
						        study presented in section 2.

simulate_method.R		R script to perform the simulation study
						        for the rotated random effects discussed in
						        section 4.

calcFC.R				R script to summarize the fraction of confounding
						    for the simulation study.

cpp_functions.R			A few C++ functions to help speed up matrix
						        algebra in R.

normality_functions.R	  Functions used to test normality and print results
						            for the simulation studies.

resid_functions.R		  Functions for "by hand" calculation of residuals

utility_functions.R		Helper functions for lme4 for the simulation studies