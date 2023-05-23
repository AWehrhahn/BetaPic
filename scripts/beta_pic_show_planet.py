from astropy.io import fits
import argparse
import matplotlib.pyplot as plt
import numpy as np
import sys

def plot(planet, model, out=None):
    if out is None:
        out = planet[:-5] + ".png"
    planet = fits.open(planet)
    model = fits.open(model)

    for chip in [1, 2, 3]:
        ext = f"CHIP{chip}.INT1"
        idata = planet[ext].data
        mdata = model[ext].data

        orders = sorted([int(c[:2]) for c in idata.names if c[-4:] == "SPEC"])

        for order in orders:
            iwave = idata[f"{order:02}_01_WL"]
            ispec = idata[f"{order:02}_01_SPEC"]
            mwave = mdata[f"{order:02}_01_WL"]
            morig = mdata[f"{order:02}_01_ORIG"]
            mspec = mdata[f"{order:02}_01_SPEC"]

            ispec[:20] =  np.nan
            ispec[-20:] =  np.nan
            mspec[:20] = np.nan
            mspec[-20:] = np.nan

            plt.clf()
            plt.plot(iwave, ispec, label="extracted")
            # plt.plot(mwave, mspec * np.nanmedian(ispec) / np.nanmedian(mspec), label="model")
            plt.plot(mwave, morig * np.nanmedian(ispec) / np.nanmedian(morig), label="model (no tellurics)")

            plt.xlabel("Wavelength [nm]")
            plt.ylabel("F")
            plt.legend()
            plt.savefig(f"beta_pic_planet_c{chip}_o{order}.png", dpi=600)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        parser = argparse.ArgumentParser()
        parser.add_argument("planet")
        parser.add_argument("model")
        parser.add_argument("-o", "--out")
        args = parser.parse_args()
        planet = args.planet
        model = args.model
        out = args.out
    else:
        planet = "/scratch/ptah/anwe5599/CRIRES/2022-11-29_L3262/extr/cr2res_obs_nodding_extractedA_planet_0.fits"
        model = "/scratch/ptah/anwe5599/CRIRES/2022-11-29_L3262/extr/beta_pic_planet_model_A.fits"
        out = None
    plot(planet, model, out=out)