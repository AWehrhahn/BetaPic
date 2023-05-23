from astropy.io import fits
import argparse
import sys
import numpy as np

if len(sys.argv) > 1:
    parser = argparse.ArgumentParser()
    parser.add_argument("observation")
    parser.add_argument("model")
    parser.add_argument("--out")
    args = parser.parse_args()
    model_fname = args.model
    observation_fname = args.observation
    outfile = args.out
else:
    date = "2022-11-29"
    setting = "L3262"
    nodding = "A"
    model_fname = f"beta_pic_img_{date}_{setting}.fits"
    observation_fname = f"/scratch/ptah/anwe5599/CRIRES/{date}_{setting}/extr/cr2res_util_calib_science_{nodding}_collapsed.fits"
    outfile = None

model = fits.open(model_fname)
observation = fits.open(observation_fname)

for chip in [1, 2, 3]:
    ext = f"CHIP{chip}.INT1"
    observation[ext].data -= model[ext].data
    observation[ext].data = np.clip(observation[ext].data, 0, None)

if outfile is None:
    outfile = observation_fname[:-5] + "_subtracted_0.fits"
# fname = f"/scratch/ptah/anwe5599/CRIRES/{date}_{setting}/extr/cr2res_util_calib_science_{nodding}_collapsed_subtracted.fits"
observation.writeto(outfile, overwrite=True)