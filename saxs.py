# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 16:50:57 2019

Functions for analyzing SAXS data and background.

@author: Andy
"""

import numpy as np
from scipy.interpolate import interp2d
import matplotlib.pyplot as plt


def average_intensity_1d(im, averaged_coord, bc=(196, 764), lambda_=0.7293E-10, 
                         pix_size=177.2E-6, d_detector=8.5028, r_lim=(10, 200),
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
        bc : (int, int)
            beam center (row, col), given by Karthik for 5-ID-D, APS
        lambda_ : float
            wavelength of X-ray [m], default corresponds to 17 keV
        pix_size : float
            size of pixels on detector [m]
        d_detector : float
            distance to detector [m]
        r_lim : (int, int)
            limits of radius to examine [pixels, size given by pix_size]
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
    xc, yc = bc
    dr = r_lim[1] - r_lim[0]
    dphi = phi_lim[1] - phi_lim[0]
    n_rows, n_cols = im.shape
    # Create Cartesian mesh grid
    x = np.arange(n_cols) - xc
    y = np.arange(n_rows) - yc # flip b/c y on computer is opposite
    # Create mesh of [r, phi] using heuristics to estimate number of points
    num_r = dr
    num_phi = round(n_rows*dphi)
    r = np.linspace(r_lim[0], r_lim[1], num=num_r)
    phi = np.linspace(phi_lim[0], phi_lim[1], num=num_phi)
    R, Phi = np.meshgrid(r, phi, indexing='ij')
    Xq = xc + np.multiply(R, np.cos(Phi))
    Yq = yc + np.multiply(R, np.sin(Phi))
    # interpolate 2D grid with cubic spline
    intensity_interp = interp2d(x, y, im, kind='cubic')
    intensity_grid = np.array([[intensity_interp(Xq[i,j],Yq[i,j]) 
                        for j in range(Xq.shape[1])] 
                        for i in range(Xq.shape[0])])
    # Average over averaged coordinate (q/r or phi)
    assert averaged_coord in ['q', 'r', 'phi'], \
        "Averaged coordinate must be q, r, or phi."
    if averaged_coord in ['q', 'r']: # same result for averaging wrt q or r
        intensity_1d = np.mean(intensity_grid, axis=0)
        free_coord = phi
    elif averaged_coord == 'phi':
        intensity_1d = np.mean(intensity_grid, axis=1)
        # convert radius of image to wave vector q
        q = 4E-10*np.pi*np.sin(0.5*np.arctan2(r*pix_size,d_detector))/lambda_
        free_coord = q    
    
    # Show region of interest over which averaging is performed
    if show_roi:
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


def parse_filename(filename):
    """
    Parses the filename of a WAXS/MAXS/SAXS TIFF image. Assumes the structure
    <hdr>_hs10#_<scan###>_<frame####>.<ext> where <hdr> is the header, hs10#
    represents the detector (hs102=WAXS, hs103=MAXS, hs104=SAXS), <scan###> is
    3-digit number for the scan and <frame####> is a 4-digit number for the
    frame.
    
    Parameters:
        filename : string
            Name of image file
    
    Returns:
        scan : int
            Number of the scan
        frame : int
            Number of the frame taken within the scan
    """
    scan = int(filename[-8:-4])
    frame = int(filename[-12:-9])
    
    return scan, frame