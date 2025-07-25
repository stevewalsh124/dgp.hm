# dgp.hm

This repository contains work on Bayesian hierarchical models for cosmology data using deep Gaussian processes and principal components analysis.

Maintainers: Steve Walsh (<walshst@elms.edu>) and Annie Booth (<annie_booth@ncsu.edu>)

## Data Description 

The Mira-Titan IV dataset is comprised of output from a computer model where the inputs represent a particular cosmology (i.e., a 9-dimensional parameter setting). Three different types of model simulations were obtained (see Moran et. al 2022) with varying degrees of computational cost:

* an inexpensive model based on perturbation theory
* an expensive simulation model with lower resolution
* a very expensive model with higher resolution

Each of these three models produce an output for the dark matter power spectrum for the given inputs (cosmology). This output can be viewed as a function over the wavenumber, $k$. Each of these models have different values for $k$ where the estimated spectra will be useful:

* Perturbation theory: $0.001000 \leq k <  0.04$
* Low resolution: $0.04 \leq k < 0.25$
* High resolution: $0.04 \leq k \leq 5$

Additionally, we leverage estimates of the precision of these outputs at different $k$ values based on the findings from Moran et. al 2022 to estimate the posterior mean of the underlying power matter spectrum for a given cosmology, conditional on all of the model outputs available.

## Modeling

The modeling process is conducted in stages as follows:

1. First, download and install the `dgp.hm` R package which is stored in the `dgp.hm` folder.  Terminal commands are provided below.

2. We fit a Bayesian hierarchical model (GP or DGP) to each individual cosmology and store the predicted means and variances.  This requires reading in the cosmology data, calculating a weighted average, scaling the precisions accordingly, estimating a covariance for the low resolution runs, and building this into a complete block covariance matrix.  For DGP models, we run one cosmology for a lot of MCMC iterations (in the `fitting/initialize.R` script) and use these burned-in values to initialize the other cosmologies.  Codes for fitting the models are in the `fitting/run_fit.R` script.  This script makes use of the `dgp.hm` R package which is stored in the like-named folder.  Results are stored in the `fitting/results` folder.

3. We collect just the posterior predicted means (`fitting/collect_means.R`) and save them in two separate csv files, one for training (`fitting/results/post_means_train.csv`) and one for testing (`fitting/results/post_means_test.csv`).

4. We conduct PCA on all the posterior predictive means from the training dataset (in the `pca/fit_basis.R` script).  We keep 10 principal components, and store the basis functions to csv files (in the `pca/bases` folder).  We also store the basis function weights for each training cosmology (`pca/bases/train_weights.csv`).

5. For each of the 10 utilized principal components, we fit a GP surrogate model on the 8-dimensional cosmology parameters (in the `Mira-Titan-IV-Data/design_train.txt` file) with the corresponding basis function weights as the response.  We then use these trained surrogates to predict the corresponding basis weights for the held-out testing cosmologies (whose 8-dimensional input values are stored in `Mira-Titan-IV-Data/design_test.txt`).  We write the predicted weights to the `pca/basis/test_weights.csv` file.  Code for this is in the `pca/predict.R` script.

6. Finally, we use these predicted weights to regenerate the predicted curves for the held-out testing cosmologies.  We write the results to the `pca/results` folder.  Code for this is in the `pca/get_curves.R` script.

For all of these steps, R codes that may be used to generate relevant plots is provided in the `plotting` folder.

## Folder Organization

A summary of this repository's folder structure is provided below.

* `CosmicEmu_etc`: contains the predictions from the Cosmic Emu model for the 6 testing cosmologies (used for comparison of final predictions)
* `Mira-Tital-IV-Data`: contains the original cosmology data, including 8-dimensional specifications, pert theory, low res runs, high res runs, and precision info
* `dgp.hm`: contains the `dgp.hm` R package which is used for fitting the hierarchical models
* `fitting`: contains R scripts and results for fitting the Bayesian hierarchical models
* `paper`: contains files for the writing of our journal article
* ` pca`: contains R scripts and results for fitting the PCA models
* `plotting`: contains R scripts to generate relevant plots (for all steps above)

## Helpful Terminal Commands

To install the `dgp.hm` R package, run
```
R CMD build dgp.hm
R CMD INSTALL dgp.hm_0.1.0.tar.gz
```

To knit the paper, run
```
cd paper
pdflatex main.tex
bibtex main
pdflatex main.tex
```
You may need to repeat these commands to make sure the reference information is compiled fully.

