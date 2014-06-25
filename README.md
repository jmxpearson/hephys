## Readme

### Introduction

These files perform analysis on electrophysiology data taken while subjects perform the balloon analogue risk task (BART). Most data handling is done in Python, with some plotting and specialized analysis done in R. Data are stored in an HDF5 file.

### Organization: Python

* `build_db.py` defines a `DataSets` class that constructs the hdf5 database from Matlab and Plexon source files. Individual files for each task (in subdirectories) use this file in constructing databases.

* `physutils.py` is a module containing utility functions related to data fetching and processing. Includes code for selecting all data relevant to a dataset tuple, bandlimiting, smoothing, splitting data around events, and binning spikes.

* `physclasses.py` provides an LFP class wrapper for a Pandas dataframe. Includes useful functions from `physutils.py` as methods.

## bartc

* `build_bart_db.py` uses the `DataSets` class from `build_db.py` to construct the BART database.

* `correlation_in_time.py` uses local eigenvalue decomposition of the covariance matrix between channels to extract a time-varying measure of the effective dimensionality of the dynamical system underlying LFP measures.

* `corrplot.py` uses multidimensional scaling to plot the relationships among LFP channels as a function of time, using cross-channel correlation as the distance measure.

* `prep_classifier_data.py` collects multichannel field potentials for each dataset and processes them into a form suitable for fitting predictive models of behavior. This processing involves constructing a matrix of integrated power in each channel in specified frequency bands for each trial (to use as regressors). Data are written to csv files.

* `plot_channel_traces.py` constructs a peri-event line plot of all channels, averaged across trial.

* `plot_LFP_channel_raster.py` constructs a heatmap plot of peri-event LFP power for a single frequency band and channel. Each row in the plot is a trial.

### Organization: R

## bartc

* `do_lfpglm_analysis.R` calls the function in `dolfpglm.R` to run a glmnet regression on data for a specific dataset. Saves fit objects to disk.

* `process_lfpglm_analysis.R` loads the fit objects, calculates performance (captured in an ROC curve), and produces a plot (saved to disk). Details of analysis, data partitioning, and cross-validation are in `run_lfp_glm.R`. Uses miscellaneous functions defined in `helpers.R` and `setup_env.R`.

### Dependencies:

* Python: NumPy, SciPy, Pandas, h5py (for reading Matlab files), pytables (via Pandas), warnings, rpy2 (for calling R).

* R: glmnet, ggplot2, reshape.

* Task-specific code needs to import `build_db.py`, `physutils.py`, etc., so these need to be somewhere in PYTHONPATH.

