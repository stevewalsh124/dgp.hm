# dgp.hm

Steps for reproducing results:

### 1) `1D_real_study_full_emuspace.R`

Get estimates of posterior mean for each cosmology using DGPs within a hierarchical model

### 2) `post_means.R`

After running step 1 for reach of the 111 training models and 6 testing models, gather these estimated means into a two dataframes (train and test) using this script.

### 3) `predictions_from_postmeans.R`

With the trained data frame, get predictions for the 6 holdout cosmologies.

### 4) `compare_pred_w_emu.R`

After obtaining predictions with our method, look at how these predictions compare to CosmicEmu predictions.

