# dgp.hm

This repository contains work on Bayesian hierarchical models for cosmology data (including Gaussian processes and deep Gaussian processes).

Maintainers: Steve Walsh (<walshst@elms.edu>) and Annie Booth (<annie_booth@ncsu.edu>)

## Data Description 

@Steve - add some comments here on the data (pert, low res, high res, when they are used, what the precision data is, etc.)

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
* `Mira-Tital-IV-Data`: contains the original cosmology data, including 8-dimensional specifications, pert theory, low res runs, high res runs, and precision info (@Steve - not all of these scripts are actually pushed up to the repo)
* `OLD`: this is an archive folder that contains old scripts for reference only (@Steve - can we delete these now?)
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

## Things still on the TO-DO list

* Identify what we should use as the "true" value for the testing cosmologies (right now, we are using the weighted average)
* Decide on whether to incorporate other competitors besides CosmicEmu (deep process convolutions, loess smoothers?)


