# Debiasing_SMI
This repository contains data source and the code reproducing the example of semi-modular inference on Covid dataset.

example.RData contains two data frames. We use "southeast_df" in our example. The data was first provided in "https://github.com/alan-turing-institute/jbc-turing-rss-testdebiasing".

MCMC.R contains the code to sample from smi posterior using MCMC and to produce figures in Appendix B.

inversion method.R contains the to code to sample from smi posterior using inversion method, corresponding to Appendix E.

prior-to-posterior probability.R is the code for section 5.

WAIC.R is the code to choose the optimal influence parameter using WAIC, mentioned at the end of section 5.

synthetic data.R is the code to check elpd in Appendix C.

section6.R and section7.R are for section 6 and section 7, respectively.

Covid_SMI.pdf is the complete paper for the Covid SMI project.

There are also two additional papers related to Cut model and SMI.
Nested_tempering.pdf analyzes Normal dataset under Cut model using both nested MCMC and tempering MCMC and compares the result in terms of precision and efficiency.

Normal_SMI.pdf contains the analysis of Normal dataset under semi-modular inference scheme.
