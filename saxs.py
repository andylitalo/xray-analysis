# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 16:50:57 2019

Functions for analyzing SAXS data and background.

@author: Andy
"""

import numpy as np
import matplotlib.pyplot as plt

#from scipy.interpolate import RectBivariateSpline
from scipy.signal import medfilt2d



def average_intensity_1d(im, averaged_coord, center=(196, 764), lambda_=0.7293E-10, 
                         pixel_size=177.2E-6, d_detector=8.5028, r_lim=(10, 200),
                         phi_lim=(-np.pi, np.pi), show_roi=False):
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
        center : (int, int)
            beam center (row, col), given by Karthik for 5-ID-D, APS
        lambda_ : float
            wavelength of X-ray [m], default corresponds to 17 keV
        pixel_size : float
            size of pixels on detector [m]
        d_detector : float
            distance to detector [m]
        r_lim : (int, int)
            limits of radius to examine [pixels, size given by pixel_size]
        phi_lim : (float, float)
            limits of angles to examine from (-pi, pi) [rad].
            *Note that the y-axis points downwards in images! (so phi=0 points
            left and phi=pi/2 points down)

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
        intensity_1d, phi = phi_profile(im, center, phi_lim)
        free_coord = phi
    elif averaged_coord == 'phi':
        intensity_1d, r = radial_profile(im, center, r_lim)
        free_coord = compute_q(r, pixel_size, d_detector, lambda_)

    # Show region of interest over which averaging is performed
    if show_roi:
        # Create Cartesian mesh grid
        x = np.arange(n_cols) - xc
        y = np.arange(n_rows) - yc # flip b/c y on computer is opposite
        # create meshgrid of Cartesian coordinates
        Y, X = np.meshgrid(y, x, indexing='ij')
        # Get polar coordinates of pixels
        R_roi = np.sqrt(X**2 + Y**2)
        Phi_roi = np.arctan2(Y,X)
        # black out region outside of region of interest
        im_roi = np.copy(im) # prevent overwriting original image
        black_out = np.bitwise_or(np.bitwise_or(R_roi > r_lim[1], R_roi < r_lim[0]), \
                                   np.bitwise_or(Phi_roi > phi_lim[1], Phi_roi < phi_lim[0]))
        im_roi[black_out] = 0
        im_roi[np.bitwise_not(black_out)] = 255
        # show image
        plt.figure()
        plt.imshow(im)
        plt.figure()
        plt.imshow(im_roi)
        
    return intensity_1d, free_coord


def average_intensity_scan(scan_file_list, averaged_coord, center=(196, 764), 
                           lambda_=0.7293E-10, pixel_size=177.2E-6, 
                           d_detector=8.5028, r_lim=(10, 200),
                         phi_lim=(-np.pi, np.pi), show_roi=False):
    """
    Averages the intensity for multiple scans over an averaged coordina
    (azimuthal angle  phi to get the intensity as a function of the wave vector
    q, or vice-versa).

    Parameters:
        scan_file_list: list of strings
            List of file names corresponding to a single scan
        averaged_coord : string, either "phi" or "q"
            Intensity is averaged over this coordinate. The other coordinate is
            returned with the intensity.
        center : (int, int)
            beam center (row, col), given by Karthik for 5-ID-D, APS
        lambda_ : float
            wavelength of X-ray [m], default corresponds to 17 keV
        pixel_size : float
            size of pixels on detector [m]
        d_detector : float
            distance to detector [m]
        r_lim : (int, int)
            limits of radius to examine [pixels, size given by pixel_size]
        phi_lim : (float, float)
            limits of angles to examine from (-pi, pi) [rad].
            *Note that the y-axis points downwards in images! (so phi=0 points
            left and phi=pi/2 points down)

    Returns:
        scan_mean : 1D array of floats
            intensity [a.u.] as a function of q averaged over phi
        free_coord : 1D array of floats
            Corresponding values of free coordinate (the one that was not
            averaged over). Either q (wave vector of X-rays [1/Angstrom]) or
            phi (angle [rad]).
    """    
    for i in range(len(scan_file_list)):
        file = scan_file_list[i]
        scan, frame = parse_filename(file)
        im = medfilt2d(plt.imread(file).astype(float))
        # analyze I vs. q
        intensity, free_coord = average_intensity_1d(im, averaged_coord, 
                                                     center=center,
                                                     lambda_=lambda_,
                                                     pixel_size=pixel_size,
                                                     d_detector=d_detector,
                                                     r_lim=r_lim,
                                                     phi_lim=phi_lim,
                                                     show_roi=show_roi)
        # Initialize mean intensity of the scan on first loop
        if i == 0:
            scan_total = np.zeros_like(intensity)
        # Sum intensities of all scans
        scan_total += intensity
    # Average scan
    scan_mean = scan_total / len(scan_file_list)
    
    return scan_mean, free_coord


def compute_q(r, pixel_size, d_detector, lambda_):
    """
    Computes the wave vector q given the radius [pixels] of the pixel,
    distance from detector, pixel size, and wavelength.
    
    Parameters:
        r : float
            Radius of location in image [pixels]
        pixel_size : float
            size of pixels on detector [m]
        d_detector : float
            distance to detector [m]
        lambda_ : float
            wavelength of X-ray [m], default corresponds to 17 keV
            
    Returns:
        q : float
            Wave vector [1/Angstrom]
    """
    # 1E-10 converts to angstroms from meters
    return (4*np.pi)*np.sin(0.5*np.arctan2(r*pixel_size,d_detector))/lambda_*1E-10


def parse_filename(filename):
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


def phi_profile(data, center, phi_lim):
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
    # Compute coordinates
    Y, X = np.indices((data.shape))
    Phi = np.arctan2(Y, X)
    # convert to 180 to have enough phi for each angle
    Phi *= 180/np.pi
    Phi = Phi.astype(np.int)
    # Bin data by radius and average
    bin_total = np.bincount(Phi.ravel(), data.ravel())
    num_phi = np.bincount(Phi.ravel())
    phi_profile = bin_total / num_phi # average
    # Limit indices to desired range of phi (and convert back to radians)
    phi = np.unique(Phi/(180/np.pi))
    inds = np.bitwise_and(phi >= phi_lim[0], phi <= phi_lim[1])
    
    return phi_profile[inds], phi[inds]


def radial_profile(data, center, r_lim):
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
    Y, X = np.indices((data.shape))
    R = np.sqrt((X - center[0])**2 + (Y - center[1])**2)
    R = R.astype(np.int)
    # Bin data by radius and average
    bin_total = np.bincount(R.ravel(), data.ravel())
    num_r = np.bincount(R.ravel())
    radial_profile = bin_total / num_r # average
    # Limit indices to desired range of r
    r = np.unique(R)
    inds = np.bitwise_and(r >= r_lim[0], r <= r_lim[1])
    
    return radial_profile[inds], r[inds]