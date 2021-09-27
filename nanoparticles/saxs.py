# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 16:50:57 2019

Functions for analyzing SAXS data and background.

@author: Andy
"""

import numpy as np
import matplotlib.pyplot as plt

from scipy.signal import medfilt2d

import saxsplot



def average_intensity_1d(im, averaged_coord, center=(196, 764), lambda_=0.7293E-10, 
                         pixel_size=177.2E-6, d_detector=8.5028, r_lim=(10, 200),
                         phi_lim=(-np.pi, np.pi), show_roi=False, show_im=False):
    """
    Integrates the intensity given in the image "im" over either the azimuthal 
    angle  phi to get the intensity as a function of the wave vector q, or 
    vice-versa.

    Parameters:
        im : 2D array
            2D pixel array (image) of intensities from SAXS
        averaged_coord : string, either "phi" or "q"
            Intensity is averaged over this coordinate. The other coordinate is
            returned with the intensity.
        center : (int, int), optional
            beam center (row, col), given by Karthik for 5-ID-D, APS
        lambda_ : float, optional
            wavelength of X-ray [m], default corresponds to 17 keV
        pixel_size : float, optional
            size of pixels on detector (includes binning) [m]. Default is for 
            a binning of 4 (177.2 um)
        d_detector : float, optional
            distance to detector [m]
        r_lim : (int, int), optional
            limits of radius to examine [pixels, size given by pixel_size]
        phi_lim : (float, float), optional
            limits of angles to examine from (-pi, pi) [rad].
            *Note that the y-axis points downwards in images! (so phi=0 points
            left and phi=pi/2 points down)
        show_roi : bool, optional
            If True, shows region of interest bounded by r_lim and phi_lim.
        show_im : bool, optional
            If True, shows image being analyzed.

    Returns:
        I : 1D array of floats
            intensity [a.u.] as a function of q averaged over phi
        free_coord : 1D array of floats
            Corresponding values of free coordinate (the one that was not
            averaged over). Either q (wave vector of X-rays [1/Angstrom]) or
            phi (angle [rad]).
    """
    # Extract parameters
    xc, yc = center
#    dr = r_lim[1] - r_lim[0]
#    dphi = phi_lim[1] - phi_lim[0]
    n_rows, n_cols = im.shape
    # pre-process image
    im = im.astype(float)
    im = medfilt2d(im)

    # Average over averaged coordinate (q/r or phi)
    assert averaged_coord in ['q', 'r', 'phi'], \
        "Averaged coordinate must be q, r, or phi."
    if averaged_coord in ['q', 'r']: # same result for averaging wrt q or r
        intensity_1d, phi = phi_profile(im, center, r_lim, phi_lim)
        free_coord = phi
    elif averaged_coord == 'phi':
        intensity_1d, r = radial_profile(im, center, r_lim, phi_lim)
        free_coord = compute_q(r,  lambda_=lambda_, pixel_size=pixel_size, 
                               d_detector=d_detector)

    # Show region of interest over which averaging is performed
    if show_roi:
        saxsplot.show_roi(im, center, r_lim, phi_lim)
    if show_im:
        # show image
        plt.figure()
        plt.imshow(im)
        
    return intensity_1d, free_coord


def compute_q(r, lambda_=0.7293E-10, pixel_size=177.2E-6, d_detector=8.5028):
    """
    Computes the wave vector q given the radius [pixels] of the pixel,
    distance from detector, pixel size, and wavelength.
    
    Parameters:
        r : float
            Radius of location in image [pixels]
        ***See average_intensity_1d for details on remaining parameters.        
            
    Returns:
        q : float
            Wave vector [1/Angstrom]
    """
    # 1E-10 converts to angstroms from meters
    return (4*np.pi)*np.sin(0.5*np.arctan2(r*pixel_size,d_detector))/lambda_*1E-10


def compute_stats(file_list, averaged_coord, center=(196,764), lambda_=0.7293E-10,
                pixel_size=177.2E-6, d_detector=8.5028, r_lim=(10,200),
                phi_lim=(-np.pi,np.pi), show_roi=False):
    """
    Compute mean and standard deviation in intensity along either q or phi.
    
    Parameters:
        file_list : list of strings
            List of filepaths to 2D intensity maps to compute standard deviation
        ***See average_intensity_1d for details on remaining parameters.
    
    Returns:
        mean : 1D array of floats
            Mean intensity from the files provided.
        std : 1D array of floats
            Standard deviation of intensity from the files provided.
        free_coord : 1D array of floats
            Wavenumber q of X-rays [1/A] if averaged_coord = 'phi'; azimuthal 
            angle phi [rad] if averaged_coord = 'q' or 'r'
    """
    assert len(file_list) > 0, "Must provide non-empty file list."
    num_files = len(file_list)
    for i in range(num_files):
        im = medfilt2d(plt.imread(file_list[i]).astype(float))
        # analyze I vs. q
        intensity, free_coord = average_intensity_1d(im, averaged_coord, 
                                                     center=center,
                                                     lambda_=lambda_,
                                                     pixel_size=pixel_size,
                                                     d_detector=d_detector,
                                                     r_lim=r_lim,
                                                     phi_lim=phi_lim,
                                                     show_roi=show_roi)
        # Initialize totals of squared intensity and intensity
        if i == 0:
            sq_total = np.zeros_like(intensity)
            total = np.zeros_like(intensity)
            # stop showing roi after first time
            show_roi = False
        sq_total += intensity**2
        total += intensity
    # Average squared intensity and intensity
    sq_mean = sq_total / num_files
    mean = total / num_files
    # Compute standard deviation sigma = sqrt(<X^2> - <X>^2)
    std = np.sqrt(sq_mean - mean**2)

    return mean, std, free_coord
    

def get_scan_filenames(scan_list, file_list):
    """
    Returns list of files corresponding to given scans.
    
    Parameters:
        scan_list : list of ints
            List of scan numbers whose files are desired
        file_list : list of strings
            List of filenames to search
    
    Returns: list of filenames corresponding to given scan numbers
    """
    # Insert input into list if provide is bare int
    if not isinstance(scan_list, list):
        scan_list = [scan_list]
        
    return [f for f in file_list if any(str(i) in f for i in scan_list)]
    
    
def parse_scan_trial(filename):
    """
    Parses the filename of a WAXS/MAXS/SAXS TIFF image. Assumes the structure
    "<hdr>_hs10#_<scan###>_<frame####>.<ext>"
    <hdr> is the header, hs10# represents the detector 
    (hs102=WAXS, hs103=MAXS, hs104=SAXS), <scan###> is 3-digit number for the 
    scan and <frame####> is a 4-digit number for the frame.
    
    Parameters:
        filename : string
            Name of image file
    
    Returns:
        scan : int
            Number of the scan
        frame : int
            Number of the frame taken within the scan
    """
    frame = int(filename[-8:-4]) # last 4 digits before extension are frame
    scan = int(filename[-12:-9]) # preceding 3 digits before '_' are scan
    
    return scan, frame


def phi_profile(data, center, r_lim, phi_lim):
    """
    Radially averages 2D grid of data to determine average value at
    different values of the azimuthal angle phi.
    
    Parameters:
        data : 2D array-like of floats
            Grid of data (usually intensities)
        center : (int, int)
            Tuple of ints giving the (x,y) coordinates of the center of data
        phi_lim : (int, int)
            (Minimum phi, maximum phi) to consider in profile, within [-pi,pi]
    
    Returns:
        radial_profile : 1D array-like of floats
            Average value of data for each radius
        phi : 1D array-like of floats
            Azimuthal angles [rad] corresponding to data points
            
    TODO:
        *Adjust meshing of phi based on radius so peaks at center aren't
            missing from average just because they don't fall in fine mesh
    """
    # conversion for degrees to radians
    deg2rad = np.pi/180
    # Compute coordinates
    Y, X = np.indices((data.shape))
    R = np.sqrt((X - center[0])**2 + (Y - center[1])**2)
    Phi = np.arctan2(Y - center[1], X - center[0])
    # add pi [rad] to make positive angles for convenience
    temp_offset = np.pi
    Phi += temp_offset
    Phi /= deg2rad # convert to degrees for finer discretization into ints
    Phi = Phi.astype(int)
    # Set points outside of r range to 0
    in_r_lim = np.bitwise_and(R >= r_lim[0], R <= r_lim[1])
    data[np.bitwise_not(in_r_lim)] = 0
    # Bin data by phi and average
    bin_total = np.bincount(Phi.ravel(), data.ravel())
    num_phi = np.bincount(Phi.ravel(), in_r_lim.ravel())
    phi_profile = bin_total / num_phi # average
    # Limit indices to desired range of phi (and convert back to radians)
    phi = np.arange(0, np.max(Phi)+1)*deg2rad - temp_offset
    in_phi_lim = np.bitwise_and(phi >= phi_lim[0], phi <= phi_lim[1])
    
    return phi_profile[in_phi_lim], phi[in_phi_lim]


def radial_profile(data, center, r_lim, phi_lim):
    """
    Azimuthally averages 2D grid of data to determine average value at
    different radii.
    
    Parameters:
        data : 2D array-like of floats
            Grid of data (usually intensities)
        center : (int, int)
            Tuple of ints giving the (x,y) coordinates of the center of data
        r_lim : (int, int)
            (Minimum r, maximum r) to consider in profile
    
    Returns:
        radial_profile : 1D array-like of floats
            Average value of data for each radius
        r : 1D array-like of floats
            Radii corresponding to data points
    """
    # Compute coordinates
    data = np.copy(data)
    Y, X = np.indices((data.shape))
    R = np.sqrt((X - center[0])**2 + (Y - center[1])**2).astype(int)
    Phi = np.arctan2(Y - center[1], X - center[0])
    # Set points outside of phi range to 0
    in_phi_lim = np.bitwise_and(Phi >= phi_lim[0], Phi <= phi_lim[1])
    data[np.bitwise_not(in_phi_lim)] = 0
    # Bin data by radius and average
    weights = data.ravel()
    bin_total = np.bincount(R.ravel(), weights)
    num_r = np.bincount(R.ravel(), in_phi_lim.ravel())
    radial_profile = bin_total / num_r # average
    # Limit indices to desired range of r
    r = np.arange(0, np.max(R)+1)
    in_r_lim = np.bitwise_and(r >= r_lim[0], r <= r_lim[1])
        
    return radial_profile[in_r_lim], r[in_r_lim]