# dgp.hm

# New folder structure

* **dgp.hm**: This folder contains the R package for hierarchical modeling.
* **fitting**: This folder contains the R script to fit a hierarchical model (shallow or deep) to a particular model.  A bash script offers easily parallel computation of multiple models.  Results are stored in csv files.  The `collect_means.R` script (will) collect and compile the predicted posterior means of all testing and training models.
* **plotting**: This folder contains R scripts to generate useful figures.


To build and install the `dgp.hm` R package, run

```
R CMD build dgp.hm
R CMD INSTALL dgp.hm_0.1.0
```

Steps for reproducing results.

### 1) `1D_real_study_full_emuspace.R`

This scripts estimates the posterior mean for a cosmology using DGPs and a hierarchical model. Run this for each of the 111 training cosmologies, and the 6 test cosmologies.

### 2) `post_means.R`

Once all of the means are estimated from Step 1, combine all training and test means into two dataframes.

### 3) `predictions_from_postmeans.R`

Use principal components and GPs to predict the dark matter power spectrum for each test cosmology. 

### 4) `compare_pred_w_emu.R`

Finally, we can compare our predictions for the test cosmologies with the corresponding predictions from CosmicEmu.

