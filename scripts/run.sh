#!/bin/bash
python beta_pic_create_planet_model.py 2022-11-29 L3262 A
python diff_img.py /scratch/ptah/anwe5599/CRIRES/2022-11-29_L3262/extr/cr2res_util_calib_science_A_collapsed.fits /scratch/ptah/anwe5599/CRIRES/2022-11-29_L3262/extr/beta_pic_img_A.fits
python beta_pic_create_planet_model.py 2022-11-29 L3262 B
python diff_img.py /scratch/ptah/anwe5599/CRIRES/2022-11-29_L3262/extr/cr2res_util_calib_science_B_collapsed.fits /scratch/ptah/anwe5599/CRIRES/2022-11-29_L3262/extr/beta_pic_img_B.fits
esorex cr2res_obs_nodding --subtract_interorder_column=FALSE /scratch/ptah/anwe5599/CRIRES/2022-11-29_L3262/extr/subtracted.sof
