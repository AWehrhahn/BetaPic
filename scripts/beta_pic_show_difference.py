from astropy.io import fits
import argparse
import matplotlib.pyplot as plt
import numpy as np

def plot(initial, subtracted, out=None):
    if out is None:
        out = initial[:-5] + ".png"
    initial = fits.open(initial)
    subtracted = fits.open(subtracted)
    for chip in [1, 2, 3]:
        ext = f"CHIP{chip}.INT1"
        idata = initial[ext].data
        sdata = subtracted[ext].data 

        orders = sorted([int(c[:2]) for c in idata.names if c[-4:] == "SPEC"])

        for order in orders:
            wave = idata[f"{order:02}_01_WL"]
            ispec = idata[f"{order:02}_01_SPEC"]
            sspec = sdata[f"{order:02}_01_SPEC"]

            ispec[:20] = sspec[:20] = np.nan
            ispec[-20:] = sspec[-20:] = np.nan

            plt.clf()
            plt.plot(wave, ispec / np.nanmedian(ispec) * np.nanmedian(ispec - sspec), label="initial")
            plt.plot(wave, sspec / np.nanmedian(sspec) * np.nanmedian(ispec - sspec), "-.", label="subtracted")
            plt.plot(wave, ispec - sspec, "--", label="initial - subtracted")
            # plt.plot(wave, sspec, label="subtracted")

            plt.xlabel("Wavelength [nm]")
            plt.ylabel("$\Delta$F")
            plt.legend()
            plt.savefig(f"beta_pic_difference_c{chip}_o{order}.png", dpi=600)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("initial")
    parser.add_argument("subtracted")
    parser.add_argument("-o", "--out")
    args = parser.parse_args()
    initial = args.initial
    subtracted = args.subtracted
    out = args.out
    plot(initial, subtracted, out=out)