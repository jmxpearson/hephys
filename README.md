## Readme

### Introduction

These files perform analysis on electrophysiology data taken while subjects perform the balloon analogue risk task (BART). Most data handling is done in Python, with some plotting and specialized analysis done in R. Data are stored in an HDF5 file.

### Organization: Python

* `build_db.py` constructs the hdf5 database from Matlab and Plexon source files.

* `physutils.py` is a module containing utility functions related to data fetching and processing. Includes code for selecting all data relevant to a dataset tuple, bandlimiting, smoothing, splitting data around events, and binning spikes.

* `prep_classifier_data.py` collects multichannel field potentials for each dataset and processes them into a form suitable for fitting predictive models of behavior. This processing involves constructing a matrix of integrated power in each channel in specified frequency bands for each trial (to use as regressors). Data are written to csv files.

* `plots.py` contains two gists. One constructs a peri-event power "raster" across trials, the other a peri-event line plot off all channels, averaged across trial.

### Organization: R

* `do_lfpglm_analysis.R` calls the function in `dolfpglm.R` to run a glmnet regression on data for a specific dataset. Saves fit objects to disk.

* `process_lfpglm_analysis.R` loads the fit objects, calculates performance (captured in an ROC curve), and produces a plot (saved to disk). Uses miscellaneous functions defined in `helpers.R` and `setup_env.R`.

### Dependencies:

* Python: NumPy, SciPy, Pandas, h5py (for reading Matlab files), pytables (via Pandas), warnings, rpy2 (for calling R).

* R: glmnet, ggplot2, reshape.

