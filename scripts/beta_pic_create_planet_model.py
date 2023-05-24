import argparse
import sys
from os.path import dirname, join

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.io import fits
from astropy.time import Time
from pyreduce.cwrappers import create_spectral_model, xi_zeta_tensors
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d, maximum_filter1d
from tqdm import tqdm

from exoorbit.orbit import Orbit
from exoorbit.bodies import Star, Planet


def make_index(ymin, ymax, xmin: int = 0, xmax: int = 2048, zero: int = 0):
    """Create an index (numpy style) that will select part of an image with changing position but fixed height
    The user is responsible for making sure the height is constant, otherwise it will still work, but the subsection will not have the desired format
    Parameters
    ----------
    ymin : array[ncol](int)
        lower y border
    ymax : array[ncol](int)
        upper y border
    xmin : int
        leftmost column
    xmax : int
        rightmost colum
    zero : bool, optional
        if True count y array from 0 instead of xmin (default: False)
    Returns
    -------
    index : tuple(array[height, width], array[height, width])
        numpy index for the selection of a subsection of an image
    """

    # TODO
    # Define the indices for the pixels between two y arrays, e.g. pixels in an order
    # in x: the rows between ymin and ymax
    # in y: the column, but n times to match the x index
    ymin = np.asarray(ymin, dtype=int)
    ymax = np.asarray(ymax, dtype=int)
    xmin = int(xmin)
    xmax = int(xmax)

    if zero:
        zero = xmin

    index_x = np.array(
        [np.arange(ymin[col], ymax[col] + 1) for col in range(xmin - zero, xmax - zero)]
    )
    index_y = np.array(
        [
            np.full(ymax[col] - ymin[col] + 1, col)
            for col in range(xmin - zero, xmax - zero)
        ]
    )
    index = index_x.T, index_y.T + zero

    return index


def load_ephimerides():
    # Load Ephimerides from file
    ephimerides_fname = join(dirname(__file__), "../data/ephimerides.txt")
    ephimerides = pd.read_table(ephimerides_fname, sep="\s+", comment="#")
    ephimerides["Date"] = [
        d[3:].replace("Nov", "11").replace("Dec", "12").replace("Feb", "02")
        for d in ephimerides["Date"]
    ]
    ephimerides["Date"] = [
        "{2}-{1}-{0}".format(*d.split("-")) for d in ephimerides["Date"]
    ]
    return ephimerides


def load_planet_spectrum():
    # Load Planet spectrum from file
    spectrum_fname = join(
        dirname(__file__), "../data/beta_pic_high_res_0.8_14.0_micron.dat"
    )
    spectrum = pd.read_table(
        spectrum_fname,
        sep="\s+",
        comment="#",
        header=None,
        names=["Wavelength", "F_lambda", "F_nu"],
    )
    spectrum_units = ["nm", "W/m^2/s", "erg/cm^s/s/Hz"]
    spectrum["Wavelength"] *= 1000
    return spectrum

def load_planet_spectrum_model(model):
    # Load planet spectrum from previous fit
    hdu = fits.open(model)
    return hdu

def load_observation(date, setting, nodding):
    # Load observationsm from files
    fname_img = "/scratch/ptah/anwe5599/CRIRES/{date}_{setting}/extr/cr2res_util_calib_science_{nodding}_collapsed.fits"
    fname_trace = "/scratch/ptah/anwe5599/CRIRES/{date}_{setting}/extr/cr2res_util_calib_flat_collapsed_tw.fits"
    fname_wave = "/scratch/ptah/anwe5599/CRIRES/{date}_{setting}/extr/cr2res_obs_nodding_extracted_combined_initial.fits"
    fname_slitfunc = "/scratch/ptah/anwe5599/CRIRES/{date}_{setting}/extr/cr2res_obs_nodding_slitfunc{nodding}_initial.fits"

    hdu_img = fits.open(fname_img.format(date=date, setting=setting, nodding=nodding))
    hdu_wave = fits.open(fname_wave.format(date=date, setting=setting))
    hdu_trace = fits.open(fname_trace.format(date=date, setting=setting))
    hdu_slitfunc = fits.open(
        fname_slitfunc.format(date=date, setting=setting, nodding=nodding)
    )
    return hdu_img, hdu_wave, hdu_trace, hdu_slitfunc


def load_star_spectrum(date, setting):
    fname = join(dirname(__file__), f"../data/beta_pic_spec_{date}_{setting}.txt")
    df = pd.read_table(
        fname, sep="\s+", comment="#", header=None, names=["Wave", "Flux"]
    )
    units = ["nm", "erg/cm**2/s/Hz"]
    return df


if len(sys.argv) > 1:
    parser = argparse.ArgumentParser()
    parser.add_argument("date")
    parser.add_argument("setting")
    parser.add_argument("nodding")
    parser.add_argument("--model")
    parser.add_argument("--out")
    parser.add_argument("--out_model")
    args = parser.parse_args()
    date = args.date
    setting = args.setting
    nodding = args.nodding
    model = args.model
    outfile = args.out
    outfile_model = args.out_model
else:
    date = "2022-11-29"
    setting = "L3262"
    nodding = "A"
    outfile = None
    outfile_model = None
    model = None

star = Star("beta Pic")
planet = star.planets["b"]
orbit = Orbit(star, planet)

# From https://www.aanda.org/articles/aa/pdf/2013/07/aa20838-12.pdf
# star_magnitude = 3.454 # ± 0.003(a)
# planet_magnitude =  11.15 #± 0.2
# From https://iopscience.iop.org/article/10.1088/2041-8205/722/1/L49/pdf
contrast_magnitude = 7.7  # +- 0.3
# Convert magnitude to flux ratio
# contrast = m_p - m_s = -2.5 log10(F_p / F_s)
# F_p / F_s = 10**(-2.5 * contrast)
# Units: W/m**2
flux_ratio = 10 ** (contrast_magnitude / (-2.5))
# TODO figure out the ratio between stellar model and planet model in the L' filter
# VLT/NACO instrument
# https://www.eso.org/sci/facilities/paranal/decommissioned/naco/inst/filters.html 
# with at least 50% transmission
# or use the transmission curve of the filter?
L_bounds = 3.5, 4.1

ephimerides = load_ephimerides()

spectrum = load_planet_spectrum()
if model is not None:
    hdu_model = load_planet_spectrum_model(model)

hdu_img, hdu_extract, hdu_trace, hdu_slitfunc = load_observation(
    date, setting, nodding
)
star_spectrum = load_star_spectrum(date, setting)

# Load header information
header = hdu_slitfunc[0].header
# in pixels
nodthrow = header["ESO SEQ CUMOFFSETY"]
nodthrow_arcsec = header["ESO SEQ NODTHROW"]
# in arcseconds
slit_width = header["ESO INS SLIT1 WID"]
# arcseconds per pixel
# This is recalculated using the order traces below
pixel_scale = 0.056

# Oversampling rate
# This might depend on the order the recipe are executed in,
# so we need to search for it in the header
oversample_name = [
    n[0]
    for n in header["ESO PRO REC? PARAM? NAME"].cards
    if n[1] == "extract_oversample"
][0]
oversample_name = oversample_name[:-4] + "VALUE"
oversample = int(header[oversample_name])
# assert header["ESO PRO REC2 PARAM4 NAME"] == "extract_oversample"
# oversample = int(header["ESO PRO REC2 PARAM4 VALUE"])

# instrumental broadening (assumed Gaussian)
# TODO: how large is this? Its around 10 for the tellurics in Molecfit
instrument_broadening = 3

# As from the obs nodding recipe in the cr2res pipeline
slit_length = 10
extr_width_frac = (slit_length - nodthrow_arcsec) / slit_length
if nodding == "A":
    slit_frac_bot = 0.0
    slit_frac_mid = slit_frac_bot + extr_width_frac / 2.0
    slit_frac_top = slit_frac_bot + extr_width_frac
else:
    slit_frac_top = 1.0
    slit_frac_mid = slit_frac_top - extr_width_frac / 2.0
    slit_frac_bot = slit_frac_top - extr_width_frac

# TODO: is the offset above or below the star
# The offset should be orthogonal to the stellar spectrum
# as the slit was aligned with the orbit
# planet seperation in mas
eph = ephimerides[ephimerides["Date"] == date]["Sep"][0]

# TODO: Check that the offsets are all in the right direction
# Determine the radial velocity offset of the planet
# 1: radial velocity of beta pic
wave_offset = star.radial_velocity.to_value(u.km / u.s)
# 2: barycentric correction
paranal = EarthLocation.of_site("Paranal")
barycorr = star.coordinates.radial_velocity_correction(
    obstime=Time(date), location=paranal
)
wave_offset += barycorr.to_value(u.km / u.s)
# 3: planet velocity
wave_offset += orbit.radial_velocity_planet(Time(date)).to_value(u.km / u.s)


primary = hdu_extract[0]
hdus = [primary]

hdus_planet_model = [primary]

fig = plt.figure(figsize=(14,4))
axs = None

fig_slitfunc = plt.figure(figsize=(14,4))
axs_slitfunc = None

fig_model = plt.figure(figsize=(14,4))
axs_model = fig_model.subplots(1, 3, gridspec_kw={"wspace":0, "hspace":0})


# iterate over the detectors
for i, chip in tqdm(enumerate([1, 2, 3]), total=3, desc="Detector"):
    ext = f"CHIP{chip}.INT1"

    img = hdu_img[ext].data
    extract_table = hdu_extract[ext].data
    trace_table = hdu_trace[ext].data
    slitfunc_table = hdu_slitfunc[ext].data
    if model is not None:
        model_table = hdu_model[ext].data

    # To store the model of the planet spectrum
    # as it would look like on the detector
    planet_img = np.zeros_like(img)
    planet_table = {}

    orders = sorted(list(set(trace_table["Order"])))

    if axs is None:
        axs = fig.subplots(len(orders), 3, sharey=True, gridspec_kw={"wspace":0, "hspace":0})
        axs_slitfunc = fig_slitfunc.subplots(len(orders), 3, sharex=True, sharey=True, gridspec_kw={"wspace":0, "hspace":0})

    for j, order in tqdm(enumerate(orders), total=len(orders), desc="Order", leave=False):
        x = np.arange(1, 2049)

        wave = extract_table[f"{order:02}_01_WL"]

        # Convert the stellar spectrum to flux density units
        spec_star = extract_table[f"{order:02}_01_SPEC"]
        mask_star = np.isfinite(spec_star)
        mask_star[:20] = False
        mask_star[-20:] = False
        spec_star[~mask_star] = np.interp(
            x[~mask_star], x[mask_star], spec_star[mask_star]
        )
        smooth_star_spec = gaussian_filter1d(maximum_filter1d(spec_star, 100), 50)
        spec_star_model = np.interp(wave, star_spectrum["Wave"], star_spectrum["Flux"])

        conversion_ratio = spec_star_model / smooth_star_spec

        slitfunc = slitfunc_table[f"{order:02}_01_SLIT_FUNC"]
        height = (slitfunc.size - 1) // oversample - 1
        slitfunc_pixels = np.linspace(0, height, slitfunc.size)

        # Get the order trace
        # But we need to shift it based on the nodding
        trace_idx = trace_table["Order"] == order
        middle_coef = trace_table[trace_idx]["All"][0, ::-1]
        lower_coef = trace_table[trace_idx]["Lower"][0, ::-1]
        upper_coef = trace_table[trace_idx]["Upper"][0, ::-1]
        middle = np.polyval(middle_coef, x)
        lower = np.polyval(lower_coef, x)
        upper = np.polyval(upper_coef, x)

        # Recalculate the pixel scale based on the measured trace
        # TODO: technically the scale changes per pixel and with the curvature
        # TODO: but not enough to matter here
        pixel_scale = np.mean(slit_length / (upper - lower))
        # planet seperation in pixels
        eph_pixel = eph / 1000 / pixel_scale

        # Shift the slit fraction to the ones used by the obs_nodding recipe
        # Thus the slit function matches the image coordinates
        interp = interp1d(
            [0, 0.5, 1], np.array([lower, middle, upper]).T, kind="linear"
        )
        lower = interp(slit_frac_bot)
        middle = interp(slit_frac_mid)
        upper = interp(slit_frac_top)

        # Get the index for the big image of this trace
        width = 2048
        height = (slitfunc.size - 1) // oversample - 1
        height_upp = height_low = height // 2
        # height_upp = int(np.ceil(np.min(upper - middle)))
        # height_low = int(np.ceil(np.min(middle - lower)))
        middle_int = middle.astype(int)
        upper_int = middle_int + height_upp
        lower_int = middle_int - height_low
        img_idx = make_index(lower_int, upper_int)
        # height = height_upp + height_low + 1

        # Calculate the offset due to the slit curvature
        slit_curv_a_coef = trace_table[trace_idx]["SlitPolyA"][0, ::-1]
        slit_curv_b_coef = trace_table[trace_idx]["SlitPolyB"][0, ::-1]
        slit_curv_c_coef = trace_table[trace_idx]["SlitPolyC"][0, ::-1]
        slit_curv_a = np.polyval(slit_curv_a_coef, x) - x
        slit_curv_b = np.polyval(slit_curv_b_coef, x)
        slit_curv_c = np.polyval(slit_curv_c_coef, x)
        slit_curv_coef = [slit_curv_c, slit_curv_b, slit_curv_a]

        # Calculate the offset of the stellar spectrum
        # The offset to the planet spectrum is calculated in the planet model
        planet_offset_y = nodthrow
        planet_offset_x = np.polyval(slit_curv_coef, planet_offset_y)

        # TODO: Check that this offset is in the right direction
        wave_planet = interp1d(x, wave, fill_value="extrapolate")(x + planet_offset_x)

        # Get the sample of the planet spectrum that we need
        # With Radial velocity offset
        spec_planet = np.interp(wave_planet, spectrum["Wavelength"] + wave_offset, spectrum["F_nu"])
        # Adjust the scale to be the same as the observed star spectrum
        # Making sure that the flux contrast matches the observed one
        spec_planet /= conversion_ratio
        spec_ratio = np.sum(spec_planet) / np.sum(spec_star)
        spec_planet *= flux_ratio / spec_ratio
        # TODO How large is the instrumental broadening?
        spec_planet = gaussian_filter1d(spec_planet, instrument_broadening)
        # Add telluric absorption lines
        spec_planet_orig = np.copy(spec_planet)
        spec_tell = np.interp(wave_planet, wave, spec_star / smooth_star_spec)
        spec_planet *= spec_tell
        
        if model is not None:
            # We still calculate the planet model first for the comparison in the plot
            # TODO: load the planet model from the first iteration, instead of recalculating
            # it again every iteration
            spec_planet = model_table[f"{order:02}_01_SPEC"]

        # Use the same slitfunctiom for the planet as for the star
        # Just shifted by the expected amount based on the ephimeredes
        # shift the slit func by the expected offset due to the planet position
        # TODO: which direction is the shift?
        iy = np.arange(slitfunc.size)
        slitfunc = np.interp(iy, iy + eph_pixel * oversample, slitfunc, left=0, right=0)

        # Can I use the calculation from PyReduce extraction to simulate the
        # planet spectrum on the detector?
        # height = int((slitfunc.size - 1) / oversample - 1)
        xi, zeta, m_zeta = xi_zeta_tensors(
            width,
            height,
            middle,
            (height // 2, height // 2),
            oversample,
            slit_curv_b,
            slit_curv_c,
        )

        # Use the efficient C method
        # this uses the xi tensor to construct a model of the planet spectrum
        # accounting for the curvature in the same way as the extraction does
        planet_img_order = create_spectral_model(
            width, height, oversample, xi, spec_planet, slitfunc
        )
        planet_img[img_idx] = planet_img_order

        planet_table[f"{order:02}_01_WL"] = wave_planet
        planet_table[f"{order:02}_01_SPEC"] = spec_planet
        planet_table[f"{order:02}_01_ORIG"] = spec_planet_orig
        planet_table[f"{order:02}_01_TELL"] = spec_tell

        # Debug Plots
        ax = axs[j, i]
        ax.plot(wave, spec_star * flux_ratio, label="star")
        ax.plot(wave_planet[20:-20], spec_planet[20:-20], label="planet")
        ax.plot(wave_planet[20:-20], spec_planet_orig[20:-20], label="no tellurics")
        if j == len(orders)-1:
            ax.set_xlabel("Wavelength [nm]")
        if i == 0 and j == len(orders) // 2:
            ax.set_ylabel("Flux [W/m^2/s]")

        ax = axs_slitfunc[j, i]
        ax.plot(slitfunc_pixels, slitfunc)
        ax.axvline(slitfunc_pixels[slitfunc.size // 2], color="k", label="mid")
        ax.axvline(slitfunc_pixels[np.argmax(slitfunc)], color="r", label="peak")
        ax.axvline(
            slitfunc_pixels[np.argmax(slitfunc)] - eph_pixel,
            linestyle="--",
            color="r",
            label="peak+offset",
        )
        ax.axvline(
            slitfunc_pixels[np.argmax(slitfunc)] + eph_pixel,
            linestyle="--",
            color="r",
            label="peak+offset",
        )
        if j == len(orders)-1:
            ax.set_xlabel("y [px]")

    hdus += [fits.ImageHDU(data=planet_img, header=hdu_img[ext].header)]
    columns = [fits.Column(name=k, format="D", array=v) for k, v in planet_table.items()]
    hdus_planet_model += [fits.BinTableHDU.from_columns(columns, header=hdu_extract[ext].header)]

    ax = axs_model[i]
    ax.imshow(planet_img, aspect="auto", origin="lower", interpolation="none")

# fig.legend()
fig.savefig("beta_pic_spectrum.png", dpi=600)
fig_slitfunc.savefig("beta_pic_slitfunc.png", dpi=600)
fig_model.savefig("beta_pic_model.png", dpi=600)


if outfile is None:
    outfile = (
        f"/scratch/ptah/anwe5599/CRIRES/{date}_{setting}/extr/beta_pic_img_{nodding}.fits"
    )

hdus = fits.HDUList(hdus)
hdus.writeto(outfile, overwrite=True)

if outfile_model is None:
    outfile_model = (
        f"/scratch/ptah/anwe5599/CRIRES/{date}_{setting}/extr/beta_pic_planet_model_{nodding}.fits"
    )
hdus_planet_model = fits.HDUList(hdus_planet_model)
hdus_planet_model.writeto(outfile_model, overwrite=True)


pass
