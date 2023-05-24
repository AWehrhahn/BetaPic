from astropy.io import fits
import argparse
import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy.ndimage import gaussian_filter1d

def plot(planet, model, out="beta_pic_planet.png", dpi=600):
    planet = fits.open(planet)
    model = fits.open(model)

    fig = plt.figure(figsize=(14,4))
    axs = None

    for i, chip in enumerate([1, 2, 3]):
        ext = f"CHIP{chip}.INT1"
        idata = planet[ext].data
        mdata = model[ext].data
        orders = sorted([int(c[:2]) for c in idata.names if c[-4:] == "SPEC"])

        if axs is None:
            axs = fig.subplots(len(orders), 3, sharex=True, sharey=True, gridspec_kw={"wspace":0, "hspace":0})

        for j, order in enumerate(orders):
            iwave = idata[f"{order:02}_01_WL"]
            ispec = idata[f"{order:02}_01_SPEC"]
            mwave = mdata[f"{order:02}_01_WL"]
            morig = mdata[f"{order:02}_01_ORIG"]
            mspec = mdata[f"{order:02}_01_SPEC"]
            mtell = mdata[f"{order:02}_01_TELL"]

            ispec = np.interp(mwave, iwave, ispec)

            ispec[:20] =  np.nan
            ispec[-20:] =  np.nan
            mspec[:20] = np.nan
            mspec[-20:] = np.nan


            # plt.clf()
            ax = axs[j, i]
            ax.plot(morig * np.nanmedian(mspec) / np.nanmedian(morig), label="model (no tellurics)", lw=1)
            ax.plot(morig * mtell * np.nanmedian(mspec) / np.nanmedian(morig), label="model (with tellurics)", lw=1)
            ax.plot(gaussian_filter1d(mspec, 3), label="input", lw=1)
            ax.plot(gaussian_filter1d(ispec, 3), label="extracted", lw=1)

            if j == len(orders)-1:
                ax.set_xlabel("Pixel [px]")
            if i == 0:
                ax.set_ylabel("F")

    plt.legend()
    outfile = out
    plt.savefig(outfile, dpi=dpi)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        parser = argparse.ArgumentParser()
        parser.add_argument("planet")
        parser.add_argument("model")
        parser.add_argument("-o", "--out")
        parser.add_argument("--dpi", default=600)
        args = parser.parse_args()
        planet = args.planet
        model = args.model
        out = args.out
        dpi = args.dpi
    else:
        planet = "/scratch/ptah/anwe5599/CRIRES/2022-11-29_L3262/extr/cr2res_obs_nodding_extractedA_planet_0.fits"
        model = "/scratch/ptah/anwe5599/CRIRES/2022-11-29_L3262/extr/beta_pic_model_A_0.fits"
        out = None
        dpi = 600
    plot(planet, model, out=out, dpi=600)