import math
import os
import sys
import fire
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
from matplotlib import cm
from scipy import signal
from skimage.feature import canny
from skimage.transform import hough_line, hough_line_peaks, rescale, resize

# TODO
# * expand rectangle around satellite trace and find associated pixels
#   * how to handle crossing events? Ideally we'd get both traces and disambiguate the crossing point
#   * how to handle flares?
# * Eliminate bad satellites based on expected trace properties
#   * get expected satellite trace length and angle (which we can get from propagating the length of the shutter speed)
# * Handle edge cases where lines extend beyond edge of image
# * Rank traces that pass the threshold and take the most central one (and mark as "ambiguous")
# * Remove stars to reduce noise
# * make sure line length algorithm is robust
#   * Take top n line start/stops and iterate through combinations, taking the one with the highest (and most consistent?) luminosity

def main(im, downscale=1, plot=False, num_satellites=1):
    # load FITS data in as an np.array
    with fits.open(im) as hdul:
        X = hdul[0].data
    # transform to log space
    X = 10 * np.log10(X)
    # handle -infs and NaNs
    # TODO just make sure there are no zeros in X
    X[~np.isfinite(X)] = np.min(X[np.isfinite(X)])
    # Shift up so all values are positive
    X -= X.min()
    image = X
    # optionally downscale
    if downscale != 1:
        image = resize(
            image,
            (image.shape[0] // downscale, image.shape[1] // downscale),
            anti_aliasing=True,
        )
    # use basic threshold to create binary image.
    # TODO is 90% background a safe assumption?
    edges = image > np.quantile(image.flatten(), 0.90)
    # test all posible angles with a 0.5 degree discretization
    # TODO limit angle search
    tested_angles = np.linspace(0, 2 * np.pi, 720, endpoint=False)
    # do hough transform and find most prominant line
    h, theta, d = hough_line(edges, theta=tested_angles)
    if plot:
        fig, axes = plt.subplots(1, 3, figsize=(15, 6))
        ax = axes.ravel()
        ax[0].imshow(image, cmap=cm.gray)
        ax[0].set_title("Input image")
        ax[0].set_axis_off()
        angle_step = 0.5 * np.diff(theta).mean()
        d_step = 0.5 * np.diff(d).mean()
        bounds = [
            np.rad2deg(theta[0] - angle_step),
            np.rad2deg(theta[-1] + angle_step),
            d[-1] + d_step,
            d[0] - d_step,
        ]
        ax[1].imshow(np.log(1 + h), extent=bounds, cmap=cm.gray, aspect=1 / 1.5)
        ax[1].set_title("Hough transform")
        ax[1].set_xlabel("Angles (degrees)")
        ax[1].set_ylabel("Distance (pixels)")
        ax[1].axis("image")
        ax[2].imshow(image, cmap=cm.gray)
        ax[2].set_ylim((image.shape[0], 0))
        ax[2].set_axis_off()
        ax[2].set_title("Detected lines")
    traces = []
    # TODO come up with reasonable values for distance and angle difference for satellites
    for _, angle, dist in zip(*hough_line_peaks(h, theta, d, num_peaks=num_satellites)):
        # get arbitrary point on line
        (x0, y0) = dist * np.array([np.cos(angle), np.sin(angle)])
        # find where the line intersects the top and bottom of image
        m = np.tan(angle + np.pi / 2)
        b = y0 - m * x0
        x_cross = -b / m
        y_max = image.shape[0]
        x_end = (y_max - b) / m
        x0, y0, x1, y1 = x_cross, 0, x_end, y_max
        # find indexes of line on image, with no repeated pixes
        # TODO are there definitely no repeated pixels?
        length = int(np.hypot(x1 - x0, y1 - y0))
        x, y = np.linspace(x0, x1, length), np.linspace(y0, y1, length)
        x, y = (x).astype(int), (y).astype(int)
        # only take indexes inside image
        mask = (x >= 0) & (x < X.shape[1]) & (y >= 0) & (y < X.shape[0])
        x, y = x[mask], y[mask]
        # Extract the values along the line
        z = X[y, x]
        # Estimate extent of trace from derivative of luminosity along line
        # TODO how sensitive are we to the window size? Leave constant for now
        window = 31
        dz = np.correlate(
            np.diff(np.correlate(z, np.ones(window), mode="same")),
            np.ones(window),
            mode="same",
        )
        start_ix = np.argmax(dz)
        end_ix = np.argmin(dz)
        x_start = x[start_ix]
        y_start = y[start_ix]
        x_end = x[end_ix-]
        y_end = y[end_ix-]
        # store indices and luminosity
        traces.append((x[start_ix:end_ix],y[start_ix:end_ix]),z[start_ix:end_ix])
        if plot:
            ax[2].plot([x_start, x_end], [y_start, y_end], marker="x", color="r")
    if plot:
        plt.tight_layout()
        plt.show()
        plt.plot(z, marker="x", color="k")

    return traces


if __name__ == "__main__":
    fire.Fire(main)
