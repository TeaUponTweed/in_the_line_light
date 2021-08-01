import math
import os
import sys
import fire
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
from dataclasses import dataclass
from matplotlib import cm
from scipy import signal
from skimage.feature import canny
from skimage.transform import hough_line, hough_line_peaks, rescale, resize
from dateutil import parser
from loadTLE import loadTLEFile
from skyfield import api
from locations import locations
from satFunctions import computeEphemeris
import datetime
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
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

@dataclass
class SatAngleParams:
    expected_angle_deg: float
    angle_tol_deg: float
    angle_disretization_deg: float

def _get_test_angles(sat_angle_params: SatAngleParams):
    if sat_angle_params.expected_angle_deg is None:
        tested_angles = np.linspace(0, 2 * np.pi, round(360/sat_angle_params.angle_disretization_deg), endpoint=False)
    else:
        assert sat_angle_params.angle_tol_deg < 180, "Angle tolerance is to high, just set expected_angle_deg to None"
        tested_angles = np.linspace(
            np.radians(sat_angle_params.expected_angle_deg-sat_angle_params.angle_tol_deg),
            np.radians(sat_angle_params.expected_angle_deg+sat_angle_params.angle_tol_deg),
            round(2*sat_angle_params.angle_tol_deg/sat_angle_params.angle_disretization_deg)
        )
    return tested_angles


def get_relevant_tle(norad_id,timestamp):
    return loadTLEFile('45748.txt')[0]

def cmd_line_main(im_file, plot=False):
    # load FITS data in as an np.array
    with fits.open(im_file) as hdul:
        assert len(hdul) == 1
        date_str = hdul[0].header['DATE-OBS']
        obs_loc = locations[hdul[0].header['OBSERVAT'].rstrip()]
        exposure_time_s = float(hdul[0].header['EXPOSURE'])
        X = hdul[0].data
        # TODO hook this up to db
        tle = get_relevant_tle(None,date_str)
        time_start = parser.parse(date_str).replace(tzinfo=api.utc)
        time_end = time_start + datetime.timedelta(seconds=exposure_time_s)
        start = computeEphemeris(tle, obs_loc, time_start)
        end = computeEphemeris(tle, obs_loc, time_end)

        w = WCS(hdul[0].header)
        from pprint import pprint
        # pprint(start)
        c_start = SkyCoord(ra=np.array([start['ra']])*u.hour, dec=np.array([start['dec']])*u.degree)
        c_end = SkyCoord(ra=np.array([end['ra']])*u.hour, dec=np.array([end['dec']])*u.degree)

        px_x_start, px_y_start = w.world_to_pixel(c_start)
        px_x_end, px_y_end = w.world_to_pixel(c_end)

    traces,X = find_sattellites_in_image(X,1,DEFAULT_SAT_ANGLE_PARAMS,plot=plot)
    plt.imshow(X,cmap=cm.gray)
    for xx,yy in traces:
        plt.plot(xx,yy,marker='x',linestyle='')
    plt.plot([px_x_start, px_x_end], [px_y_start, px_y_end],color='y')
    plt.show()

DEFAULT_SAT_ANGLE_PARAMS = SatAngleParams(None,None,0.5)
def find_sattellites_in_image(X, num_satellites, sat_angle_params,plot=False):
    # handle zeros for transformation to log space
    X[X<=0] = X[X>0].min()
    # transform to log space
    X = 10 * np.log10(X)
    assert np.all(np.isfinite(X))
    # Shift up so all values are positive
    X -= X.min()

    # use basic threshold to create binary image.
    # TODO is 95% background a safe assumption? Probably not, 90% background breaks on a least one image.
    # * There seems to be a _lot_ of discretization at the lower levels, maybe we can make use of that instead
    edges = X > np.quantile(X.flatten(), 0.95)
    # do hough transform and find most prominant line
    h, theta, d = hough_line(edges, theta=_get_test_angles(sat_angle_params))
    if plot:
        fig, axes = plt.subplots(1, 3, figsize=(15, 6))
        ax = axes.ravel()
        ax[0].imshow(X, cmap=cm.gray)
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
        ax[2].imshow(X, cmap=cm.gray)
        ax[2].set_ylim((X.shape[0], 0))
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
        y_max = X.shape[0]
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
        x_end = x[end_ix-1]
        y_end = y[end_ix-1]
        # store indices and luminosity
        if plot:
            ax[2].plot([x_start, x_end], [y_start, y_end], marker="x", color="r")

        traces.append((x[start_ix:end_ix],y[start_ix:end_ix]))

    if plot:
        plt.tight_layout()
        plt.show()
        # plt.plot(z, marker="x", color="k")
        # plt.plot(dz, marker="", color="r",linestyle=':')
        # plt.show()

    return traces,X


if __name__ == "__main__":
    fire.Fire(cmd_line_main)
