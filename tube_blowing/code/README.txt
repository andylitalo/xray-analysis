CODE 

This folder contains the scripts and functions for analyzing scattering patterns from Argonne's APS.

-IvsPhi_Calculator: function for integrating the 2D intensity with respect to q (so that the result is intensity vs. phi)

-IvsQ_Calculator: function for integrating the 2D intensity with respect to phi (so that the result is intensity vs. q)

-plot_1d: script for plotting comparisons of 1d-integrated intensity between selected frames. Uses data saved in "datasheets" folder and saves plots in "1d_select" folder.

-Q_Matrix_Calculator: function whose output is used by the two-parameter subtraction function as an input.

-Two_Para_Subtraction: function that adjusts the background to match the signal at two locations where neither the signal nor the background have a strong effect. It is useful if the true background varies from the measured background (e.g., when scanning the glass mold during the experiment we might look through a different part then we do during the background scan.

-Xray_Analysis: script for analyzing SAXS or WAXS data. Please adapt the parameters at the top to your specific dataset. Be sure to load your 2D scattering patterns in  the "2d_raw" folder and load the relevant part of the header file (from the "specfiles" folder in the beamline data) in an Excel worksheet (.xlsx) in the "datasheets" folder before starting. It saves data in the "1d_integrated" and "datasheets" folders. 