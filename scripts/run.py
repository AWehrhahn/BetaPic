from astropy.io import fits
import numpy as np
import pandas as pd
from os.path import join, dirname
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt

from pyreduce.cwrappers import xi_zeta_tensors, create_spectral_model
from tqdm import tqdm


def make_index(
    ymin, ymax, xmin: int = 0, xmax: int = 2048, zero: int = 0
):
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
        [
            np.arange(ymin[col], ymax[col] + 1)
            for col in range(xmin - zero, xmax - zero)
        ]
    )
    index_y = np.array(
        [
            np.full(ymax[col] - ymin[col] + 1, col)
            for col in range(xmin - zero, xmax - zero)
        ]
    )
    index = index_x.T, index_y.T + zero

    return index

ephimerides_fname = join(dirname(__file__), "../data/ephimerides.txt")
ephimerides = pd.read_table(ephimerides_fname, sep="\s+", comment="#")
ephimerides["Date"] = [
    d[3:].replace("Nov", "11").replace("Dec", "12").replace("Feb", "02")
    for d in ephimerides["Date"]
]
ephimerides["Date"] = ["{2}-{1}-{0}".format(*d.split("-")) for d in ephimerides["Date"]]

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

fname_img = "/scratch/ptah/anwe5599/CRIRES/{date}_{setting}/extr/cr2res_util_calib_science_{nodding}_collapsed.fits"
fname_wave = "/scratch/ptah/anwe5599/CRIRES/{date}_{setting}/extr/cr2res_obs_nodding_extracted_combined.fits"
fname_trace = "/scratch/ptah/anwe5599/CRIRES/{date}_{setting}/extr/cr2res_util_calib_flat_collapsed_tw.fits"
fname_slitfunc = "/scratch/ptah/anwe5599/CRIRES/{date}_{setting}/extr/cr2res_obs_nodding_slitfunc{nodding}.fits"


date = "2022-11-29"
setting = "L3262"
nodding = "A"

hdu_img = fits.open(fname_img.format(date=date, setting=setting, nodding=nodding))
hdu_wave = fits.open(fname_wave.format(date=date, setting=setting))
hdu_trace = fits.open(fname_trace.format(date=date, setting=setting))
hdu_slitfunc = fits.open(
    fname_slitfunc.format(date=date, setting=setting, nodding=nodding)
)


header = hdu_slitfunc[0].header
# in pixels
nodthrow = header["ESO SEQ CUMOFFSETY"]
nodthrow_arcsec = header["ESO SEQ NODTHROW"]
# in arcseconds
slit_width = header["ESO INS SLIT1 WID"]
# arcseconds per pixel
pixel_scale = 0.056
# Oversampling rate
assert header["ESO PRO REC2 PARAM4 NAME"] == "extract_oversample"
oversample = int(header["ESO PRO REC2 PARAM4 VALUE"])

# instrumental broadening (assumed Gaussian)
# TODO: how large is this? Its around 10 for the tellurics in Molecfit
instrument_broadening = 10

# As from the obs nodding recipe
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
# planet seperation in pixels
eph_pixel = eph / 1000 / pixel_scale

# TODO iterate over the detectors
chip = 1
ext = f"CHIP{chip}.INT1"

img = hdu_img[ext].data
wave_table = hdu_wave[ext].data
trace_table = hdu_trace[ext].data
slitfunc_table = hdu_slitfunc[ext].data

# To store the model of the planet spectrum
# as it would look like on the detector
planet_img = np.zeros_like(img)

# TODO iterate over the orders
orders = sorted(list(set(trace_table["Order"])))
for order in orders:
    wave = wave_table[f"{order:02}_01_WL"]
    slitfunc = slitfunc_table[f"{order:02}_01_SLIT_FUNC"]
    slitfunc_pixels = np.arange(slitfunc.size) / oversample

    # Get the order trace
    # But we need to shift it based on the nodding
    x = np.arange(1, 2049)
    trace_idx = trace_table["Order"] == order
    middle_coef = trace_table[trace_idx]["All"][0, ::-1]
    lower_coef = trace_table[trace_idx]["Lower"][0, ::-1]
    upper_coef = trace_table[trace_idx]["Upper"][0, ::-1]
    middle = np.polyval(middle_coef, x)
    lower = np.polyval(lower_coef, x)
    upper = np.polyval(upper_coef, x)

    # Shift the slit fraction to the ones used by the obs_nodding recipe
    # Thus the slit function matches the image coordinates
    interp = interp1d([0, 0.5, 1], np.array([lower, middle, upper]).T, kind="linear")
    lower = interp(slit_frac_bot)
    middle = interp(slit_frac_mid)
    upper = interp(slit_frac_top)

    # Get the index for the big image of this trace
    width = 2048
    height = int((slitfunc.size - 1) / oversample - 1)
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
    spec = np.interp(wave_planet, spectrum["Wavelength"], spectrum["F_lambda"])
    spec = gaussian_filter1d(spec, instrument_broadening)

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
    # this uses the xi vector to construct a model of the planet spectrum
    # accounting for the curvature in the same way as the 
    # reduction does
    planet_img_order = create_spectral_model(width, height, oversample, xi, spec, slitfunc)
    planet_img[img_idx] = planet_img_order

# Debug Plots
plt.clf()
plt.imshow(planet_img, aspect="auto", origin="lower", interpolation="none")
plt.savefig("beta_pic_planet_model.png")

plt.clf()
plt.plot(wave_planet, spec)
plt.xlabel("Wavelength [nm]")
plt.ylabel("Flux [W/m^2/s]")
plt.savefig("beta_pic_spectrum.png")

plt.clf()
plt.plot(slitfunc_pixels, slitfunc)
plt.axvline(slitfunc_pixels[slitfunc.size // 2], color="k", label="mid")
plt.axvline(slitfunc_pixels[np.argmax(slitfunc)], color="r", label="peak")
plt.axvline(
    slitfunc_pixels[np.argmax(slitfunc)] - eph_pixel,
    linestyle="--",
    color="r",
    label="peak+offset",
)
plt.xlabel("y [px]")
plt.savefig("beta_pic_slitfunc.png")

pass
