# -*- coding: utf-8 -*-
"""
Created on Tue May  7 15:12:09 2019

Plotting functions for SAXS data.

@author: Andy
"""

import numpy as np
import matplotlib.pyplot as plt


def plot(x_list, y_list, x_label='q (1/A)', y_label='intensity (a.u.)',
         label_list=[], title='SAXS Intensity vs. Wavenumber', xscale='log',
         yscale='log', ax_fs=18,t_fs=20, tk_fs=14, l_fs=14, save_path=''):
    """
    Plots the intensity [a.u.] vs. the wave number [1/A] for the given
    datasets.

    Parameters:
        x_list : list of 1D array-like of floats
            x-values (often wavenumber q [1/A]).  Each entry is a new dataset
        y_list : list of 1D array-like of floats
            y-values (often intensities [a.u.]). Each entry is a new dataset
        x_label : string, optional
            Label for x-axis. Default 'q (1/A)'
        y_label : string
            Label for y-axis. Default 'intensity (a.u.)'
        label_list : list of strings, optional
            Labels for legend (if plotting multiple datasets). Default []
        title : string, optional
            Title of plot. Default 'SAXS Intensity vs. Wavenumber'
        xscale : string, optional
            Scale of x axis: 'log' -> log scale, 'linear' -> linear scale
        yscale : string, optional
            Scale of y axis: 'log' -> log scale, 'linear' -> linear scale
        ax_fs : int, optional
            Font size of axis labels
        t_fs : int, optional
            Font size of title
        tk_fs : int, optional
            Font size of tick labels
        l_fs : int, optional
            Font size of legend
        save_path : string
            File path to save axis (w/ folders and ext). Saves if not empty.

    Returns:
        ax : axis handle
            Handle of current axis object
    """
    # Insert input into list if provided is bare array
    if not isinstance(x_list, list):
        x_list = [x_list]
    if not isinstance(y_list, list):
        y_list = [y_list]
    # Set up figure and plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel(x_label, fontsize=ax_fs)
    ax.set_ylabel(y_label, fontsize=ax_fs)
    ax.set_title(title, fontsize=t_fs)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    # Plot data
    assert len(x_list) == len(y_list), "Must have same number of x series as y series."
    for i in range(len(x_list)):
        line, = ax.plot(x_list[i], y_list[i])
        if label_list:
            line.set_label(label_list[i])
    # Show legend if labels provided
    if label_list:
        ax.legend(loc='best', fontsize=l_fs)
    # increase font size of tick labels
    ax.tick_params(axis='both', which='major', labelsize=tk_fs)
    # Save figure
    if save_path:
        plt.savefig(save_path, bbox_inches='tight')

    return ax



def show_roi(im, center, r_lim, phi_lim):
    """
    Shows the region of interest in a SAXS plot.

    Parameters:
        im : 2D array
            2D pixel array (image) of intensities from SAXS
        center : (int, int)
            beam center (row, col), given in header file (folder "processing")
         r_lim : (int, int)
            limits of radius to examine [pixels, size given by pixel_size]
        phi_lim : (float, float)
            limits of angles to examine from (-pi, pi) [rad].
            *Note that the y-axis points downwards in images! (so phi=0 points
            left and phi=pi/2 points down)
    """
    im_roi = np.copy(im) # prevents overwriting original image

    # Define ROI
    n_rows, n_cols = im.shape
    xc, yc = center
    # Create Cartesian mesh grid
    x = np.arange(n_cols) - xc
    y = np.arange(n_rows) - yc # flip b/c y on computer is opposite
    # create meshgrid of Cartesian coordinates
    Y, X = np.meshgrid(y, x, indexing='ij')
    # Get polar coordinates of pixels
    R_roi = np.sqrt(X**2 + Y**2)
    Phi_roi = np.arctan2(Y,X)
    # black out region outside of region of interest
    black_out = np.bitwise_or(np.bitwise_or(R_roi > r_lim[1], R_roi < r_lim[0]), \
                               np.bitwise_or(Phi_roi > phi_lim[1], Phi_roi < phi_lim[0]))
    im_roi[black_out] = 0
    im_roi[np.bitwise_not(black_out)] = 255

    # plot
    plt.figure()
    plt.imshow(im_roi)
    plt.title(r'$r \in$ [{r_min},{r_max}], $\phi \in$ [{phi_min:.2f}, {phi_max:.2f}]' \
              .format(r_min=r_lim[0], r_max=r_lim[1], phi_min=phi_lim[0], \
                      phi_max=phi_lim[1]))
    plt.xlabel('x [pixels]')
    plt.ylabel('y [pixels]')
