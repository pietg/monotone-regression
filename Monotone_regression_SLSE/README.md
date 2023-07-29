# monotone-regression
R scripts for confidence intervals in monotone regression, based on the SLSE

This repository acccompanies the manuscript "Confidence intervals in monotone regression", with authors Piet Groeneboom and Geurt Jongbloed (SLSE.pdf).

Coverage percentages of the SLSE (Smoothed Least Squares Estimator) and the Nadaraya-Watson estimator can be computed by running bootstrap_intervals_percentages.R.

A boxplot of 1000 values of the estimates of f_0(0.6) for both estimators are produced by boxplot_simulation.R (sample size 500).

The bandwidth selection, as described in Section 4 of the manuscript SLSE.pdf is demonstrated by the R script bandwidth.R.

In all instances one can see the output by typing "output" (for example the minimizing constant c in h=cn^{-1/5} for bandwidth.R) in the console window.
