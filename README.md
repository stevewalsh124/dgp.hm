# dgp.hm

# Potential Competitors
* (Kelly will share test cases with us)
* CosmicEMU
* DPC - process convolutions
* low hanging fruit - loess, independent GPs?
* for real world data - we don't have the real known truth, could we use all per t theory?  Use the pert theory from the training/testing data, then for the training data create "fake" low and high res runs from it.

# New folder structure

* **dgp.hm**: This folder contains the R package for hierarchical modeling.
* **fitting**: This folder contains the R script to fit a hierarchical model (shallow or deep) to a particular model.  A bash script offers easily parallel computation of multiple models.  Results are stored in csv files.  The `collect_means.R` script (will) collect and compile the predicted posterior means of all testing and training models.
* **OLD**: contains old scripts for possible reference
* **pca**: contains scripts for fitting PCA model (NOT DONE YET)
* **plotting**: This folder contains R scripts to generate useful figures.

To build and install the `dgp.hm` R package, run

```
R CMD build dgp.hm
R CMD INSTALL dgp.hm_0.1.0
```


