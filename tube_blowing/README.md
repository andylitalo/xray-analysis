# Guide to Analysis of X-ray Scattering Data from Argonne's Advanced Photon Source

This folder contains the folder structure and code for analyzing small-angle
(SAXS) and wide-angle (WAXS) X-ray scattering data from the Advanced Photon
Source, beamline 5-ID-D, at Argonne National Laboratory.

# Directory Structure

- `1d_integrated`: the scripts `WAXS.m` and `SAXS.m` will save two images per
scan of the 1d integrated intensity alongside the background, one vs. the
wave number q and the vs. the azimuthal angle phi.

- `1d_select`: the script `plot_1d.m` will save its plots comparing selected
frames of 1d-integrated intensity here.

- `2d_bkgd_subt`: the scripts `WAXS.m` and `SAXS.m` will save their 2D
scattering patterns with background subtracted here.

- `2d_raw`: **Put your 2D scattering patterns here before analysis.**

- `code`: **Contains all scripts and functions for the analysis.**

- `datasheets`: the scripts `WAXS.m` and `SAXS.m` will save their processed data
as `.csv` files here. The script `plot_1d.m` will use these `.csv files` to
generate its plots.

- `metadata`: contains any additional data relevant to the experiment.
