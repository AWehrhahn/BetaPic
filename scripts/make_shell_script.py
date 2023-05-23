from os.path import join, dirname
import subprocess
import argparse
import sys

class SofFile:
    """
    Handles the 'Set Of Files' data files that are used by esorex.
    Essentially a list of files and their descriptor that is then used by esorex.

    TODO: the class could be more sophisticated using e.g. a numpy array
    or even pandas, but that would be overkill for now
    """

    def __init__(self, data=None):
        if data is None:
            data = []
        #:list: content of the SOF file
        self.data = data

    @classmethod
    def read(cls, filename):
        """
        Reads a sof file from disk

        Parameters
        ----------
        filename : str
            name of the file to read

        Returns
        -------
        self : SofFile
            The read file
        """
        data = []
        with open(filename, "r") as f:
            for line in f.readline():
                fname, ftype = line.split(maxsplit=1)
                data += [(fname, ftype)]
        self = cls(data)
        return self

    def write(self, filename):
        """
        Writes the sof file to disk

        Parameters
        ----------
        filename : str
            The name of the file to store
        """
        content = []
        for fname, ftype in self.data:
            content += [f"{fname} {ftype}\n"]

        with open(filename, "w") as f:
            f.writelines(content)

    def append(self, filename, tag):
        self.data += [(filename, tag)]

def move_obs_nodding_files(outdir, suffix):
    
    ff_name_combined_A = join(outdir, 'cr2res_obs_nodding_combinedA.fits')
    ff_name_extracted_A = join(outdir, 'cr2res_obs_nodding_extractedA.fits')
    ff_name_slitfunc_A = join(outdir, 'cr2res_obs_nodding_slitfuncA.fits')
    ff_name_model_A = join(outdir, 'cr2res_obs_nodding_modelA.fits')
    ff_name_tracewave_A = join(outdir, 'cr2res_obs_nodding_trace_wave_A.fits')
    ff_name_combined_B = join(outdir, 'cr2res_obs_nodding_combinedB.fits')
    ff_name_extracted_B = join(outdir, 'cr2res_obs_nodding_extractedB.fits')
    ff_name_slitfunc_B = join(outdir, 'cr2res_obs_nodding_slitfuncB.fits')
    ff_name_model_B = join(outdir, 'cr2res_obs_nodding_modelB.fits')
    ff_name_tracewave_B = join(outdir, 'cr2res_obs_nodding_trace_wave_B.fits')
    ff_name_extracted_AB = join(outdir, 'cr2res_obs_nodding_extracted_combined.fits')

    ff_name_combined_A_initial = join(outdir, f'cr2res_obs_nodding_combinedA_{suffix}.fits')
    ff_name_extracted_A_initial = join(outdir, f'cr2res_obs_nodding_extractedA_{suffix}.fits')
    ff_name_slitfunc_A_initial = join(outdir, f'cr2res_obs_nodding_slitfuncA_{suffix}.fits')
    ff_name_model_A_initial = join(outdir, f'cr2res_obs_nodding_modelA_{suffix}.fits')
    ff_name_tracewave_A_initial = join(outdir, f'cr2res_obs_nodding_trace_wave_A_{suffix}.fits')
    ff_name_combined_B_initial = join(outdir, f'cr2res_obs_nodding_combinedB_{suffix}.fits')
    ff_name_extracted_B_initial = join(outdir, f'cr2res_obs_nodding_extractedB_{suffix}.fits')
    ff_name_slitfunc_B_initial = join(outdir, f'cr2res_obs_nodding_slitfuncB_{suffix}.fits')
    ff_name_model_B_initial = join(outdir, f'cr2res_obs_nodding_modelB_{suffix}.fits')
    ff_name_tracewave_B_initial = join(outdir, f'cr2res_obs_nodding_trace_wave_B_{suffix}.fits')
    ff_name_extracted_AB_initial = join(outdir, f'cr2res_obs_nodding_extracted_combined_{suffix}.fits')

    commands = [
        f"mv {ff_name_combined_A} {ff_name_combined_A_initial}\n"
        f"mv {ff_name_extracted_A} {ff_name_extracted_A_initial}\n"
        f"mv {ff_name_slitfunc_A} {ff_name_slitfunc_A_initial}\n"
        f"mv {ff_name_model_A} {ff_name_model_A_initial}\n"
        f"mv {ff_name_tracewave_A} {ff_name_tracewave_A_initial}\n"
        f"mv {ff_name_combined_B} {ff_name_combined_B_initial}\n"
        f"mv {ff_name_extracted_B} {ff_name_extracted_B_initial}\n"
        f"mv {ff_name_slitfunc_B} {ff_name_slitfunc_B_initial}\n"
        f"mv {ff_name_model_B} {ff_name_model_B_initial}\n"
        f"mv {ff_name_tracewave_B} {ff_name_tracewave_B_initial}\n"
        f"mv {ff_name_extracted_AB} {ff_name_extracted_AB_initial}\n"
    ]
    return commands[0]

if len(sys.argv) > 1:
    parser = argparse.ArgumentParser()
    parser.add_argument("date")
    parser.add_argument("setting")
    parser.add_argument("outdir")
    parser.add_argument("--iterations", default=5)
    args = parser.parse_args()
    date = args.date
    setting = args.setting
    outdir = args.outdir
    iterations = args.iterations
else:
    date = "2022-11-29"
    setting = "L3262"
    outdir = "/scratch/ptah/anwe5599/CRIRES/2022-11-29_L3262/extr/"
    iterations = 5

tw_name = join(outdir, "cr2res_util_calib_flat_collapsed_tw.fits")
ff_norm_flat = join(outdir, "cr2res_util_normflat_Open_master_flat.fits")
master_dark_bpm = join(outdir, "cr2res_cal_dark_bpm.fits")

# ff_name_combined_A = join(outdir, 'cr2res_obs_nodding_combinedA.fits')
# ff_name_extracted_A = join(outdir, 'cr2res_obs_nodding_extractedA.fits')
# ff_name_slitfunc_A = join(outdir, 'cr2res_obs_nodding_slitfuncA.fits')
# ff_name_model_A = join(outdir, 'cr2res_obs_nodding_modelA.fits')
# ff_name_tracewave_A = join(outdir, 'cr2res_obs_nodding_trace_wave_A.fits')
# ff_name_combined_B = join(outdir, 'cr2res_obs_nodding_combinedB.fits')
# ff_name_extracted_B = join(outdir, 'cr2res_obs_nodding_extractedB.fits')
# ff_name_slitfunc_B = join(outdir, 'cr2res_obs_nodding_slitfuncB.fits')
# ff_name_model_B = join(outdir, 'cr2res_obs_nodding_modelB.fits')
# ff_name_tracewave_B = join(outdir, 'cr2res_obs_nodding_trace_wave_B.fits')
# ff_name_extracted_AB = join(outdir, 'cr2res_obs_nodding_extracted_combined.fits')

# ff_name_combined_A_initial = join(outdir, 'cr2res_obs_nodding_combinedA_initial.fits')
# ff_name_extracted_A_initial = join(outdir, 'cr2res_obs_nodding_extractedA_initial.fits')
# ff_name_slitfunc_A_initial = join(outdir, 'cr2res_obs_nodding_slitfuncA_initial.fits')
# ff_name_model_A_initial = join(outdir, 'cr2res_obs_nodding_modelA_initial.fits')
# ff_name_tracewave_A_initial = join(outdir, 'cr2res_obs_nodding_trace_wave_A_initial.fits')
# ff_name_combined_B_initial = join(outdir, 'cr2res_obs_nodding_combinedB_initial.fits')
# ff_name_extracted_B_initial = join(outdir, 'cr2res_obs_nodding_extractedB_initial.fits')
# ff_name_slitfunc_B_initial = join(outdir, 'cr2res_obs_nodding_slitfuncB_initial.fits')
# ff_name_model_B_initial = join(outdir, 'cr2res_obs_nodding_modelB_initial.fits')
# ff_name_tracewave_B_initial = join(outdir, 'cr2res_obs_nodding_trace_wave_B_initial.fits')
ff_name_extracted_AB_initial = join(outdir, 'cr2res_obs_nodding_extracted_combined_initial.fits')

# filenames for each iteration
# ff_name_combined_A_subtracted = join(outdir, 'cr2res_obs_nodding_combinedA_subtracted_0.fits')
# ff_name_extracted_A_subtracted = join(outdir, 'cr2res_obs_nodding_extractedA_subtracted_0.fits')
# ff_name_slitfunc_A_subtracted = join(outdir, 'cr2res_obs_nodding_slitfuncA_subtracted_0.fits')
ff_name_model_A_subtracted = join(outdir, 'cr2res_obs_nodding_modelA_subtracted_0.fits')
# ff_name_tracewave_A_subtracted = join(outdir, 'cr2res_obs_nodding_trace_wave_A_subtracted_0.fits')
# ff_name_combined_B_subtracted = join(outdir, 'cr2res_obs_nodding_combinedB_subtracted_0.fits')
# ff_name_extracted_B_subtracted = join(outdir, 'cr2res_obs_nodding_extractedB_subtracted_0.fits')
# ff_name_slitfunc_B_subtracted = join(outdir, 'cr2res_obs_nodding_slitfuncB_subtracted_0.fits')
ff_name_model_B_subtracted = join(outdir, 'cr2res_obs_nodding_modelB_subtracted_0.fits')
# ff_name_tracewave_B_subtracted = join(outdir, 'cr2res_obs_nodding_trace_wave_B_subtracted_0.fits')
ff_name_extracted_AB_subtracted = join(outdir, 'cr2res_obs_nodding_extracted_combined_subtracted_0.fits')
ff_name_extracted_A_planet = join(outdir, f"cr2res_obs_nodding_extractedA_planet_0.fits")



planet_img_A = join(outdir, 'beta_pic_img_A_0.fits')
planet_model_A = join(outdir, 'beta_pic_model_A_0.fits')
planet_img_B = join(outdir, 'beta_pic_img_B_0.fits')
planet_model_B = join(outdir, 'beta_pic_model_B_0.fits')

ff_name_science_A = join(outdir, 'cr2res_util_calib_science_A_collapsed.fits')
ff_name_science_B = join(outdir, 'cr2res_util_calib_science_B_collapsed.fits')
ff_name_science_A_subtracted = join(outdir, 'cr2res_util_calib_science_A_collapsed_subtracted_0.fits')
ff_name_science_B_subtracted = join(outdir, 'cr2res_util_calib_science_B_collapsed_subtracted_0.fits')

# Make subtracted SOF
extract_star_sof = join(outdir, "subtracted_0.sof")
extract_star_sof_data = SofFile()
extract_star_sof_data.append(ff_name_science_A_subtracted, "OBS_NODDING_OTHER")
extract_star_sof_data.append(ff_name_science_B_subtracted, "OBS_NODDING_OTHER")
extract_star_sof_data.append(tw_name, "CAL_FLAT_TW")
extract_star_sof_data.append(ff_norm_flat, "UTIL_MASTER_FLAT")
extract_star_sof_data.append(master_dark_bpm, "CAL_DARK_BPM")
extract_star_sof_data.write(extract_star_sof)

python_dir = dirname(__file__)

commands = [
    "#!/bin/bash\n"
    # Copy the initial files to a different name
    + move_obs_nodding_files(outdir, "initial") +
    # Subtract the initial model of the planet
    f"python {join(python_dir, 'beta_pic_create_planet_model.py')} {date} {setting} A --out={planet_img_A} --out_model={planet_model_A}\n"
    f"python {join(python_dir, 'diff_img.py')} {ff_name_science_A} {planet_img_A} --out={ff_name_science_A_subtracted}\n"
    f"python {join(python_dir, 'beta_pic_create_planet_model.py')} {date} {setting} B --out={planet_img_B} --out_model={planet_model_B}\n"
    f"python {join(python_dir, 'diff_img.py')} {ff_name_science_B} {planet_img_B} --out={ff_name_science_B_subtracted}\n"
    # Run the extraction of the spectrum
    f"esorex cr2res_obs_nodding --subtract_interorder_column=FALSE {extract_star_sof}\n"
    # and save the data to a new location
    + move_obs_nodding_files(outdir, "subtracted_0") +
    # Plot the results
    f"python {join(python_dir, 'beta_pic_show_difference.py')} {ff_name_extracted_AB_initial} {ff_name_extracted_AB_subtracted}\n"
    # Remove the star model from the images
    f"python {join(python_dir, 'diff_img.py')} {ff_name_science_A} {ff_name_model_A_subtracted} --out={ff_name_science_A_subtracted}\n"
    f"python {join(python_dir, 'diff_img.py')} {ff_name_science_B} {ff_name_model_B_subtracted} --out={ff_name_science_B_subtracted}\n"
    # Run the extraction again with just the planet
    f"esorex cr2res_obs_nodding --subtract_interorder_column=FALSE {extract_star_sof}\n"
    # Move the files again
    + move_obs_nodding_files(outdir, "planet_0") +
    # Plot the extracted planet spectrum
    f"python {join(python_dir, 'beta_pic_show_planet.py')} {ff_name_extracted_A_planet} {planet_model_A}\n"
    
]

for i in range(1, iterations):    
    planet_img_A = join(outdir, f'beta_pic_img_A_{i}.fits')
    planet_model_A = join(outdir, f'beta_pic_model_A_{i}.fits')
    planet_img_B = join(outdir, f'beta_pic_img_B_{i}.fits')
    planet_model_B = join(outdir, f'beta_pic_model_B_{i}.fits')
    ff_name_science_A_subtracted = join(outdir, f'cr2res_util_calib_science_A_collapsed_subtracted_{i}.fits')
    ff_name_science_B_subtracted = join(outdir, f'cr2res_util_calib_science_B_collapsed_subtracted_{i}.fits')
    ff_name_model_A_subtracted = join(outdir, f'cr2res_obs_nodding_modelA_subtracted_{i}.fits')
    ff_name_model_B_subtracted = join(outdir, f'cr2res_obs_nodding_modelB_subtracted_{i}.fits')
    ff_name_extracted_A_planet = join(outdir, f"cr2res_obs_nodding_extractedA_planet_{i-1}.fits")
    ff_name_extracted_A_planet_new = join(outdir, f"cr2res_obs_nodding_extractedA_planet_{i}.fits")

    ff_name_extracted_B_planet = join(outdir, f"cr2res_obs_nodding_extractedB_planet_{i-1}.fits")
    ff_name_extracted_AB_subtracted = join(outdir, f'cr2res_obs_nodding_extracted_combined_subtracted_{i}.fits')


    extract_star_sof = join(outdir, f"subtracted_{i}.sof")
    extract_star_sof_data = SofFile()
    extract_star_sof_data.append(ff_name_science_A_subtracted, "OBS_NODDING_OTHER")
    extract_star_sof_data.append(ff_name_science_B_subtracted, "OBS_NODDING_OTHER")
    extract_star_sof_data.append(tw_name, "CAL_FLAT_TW")
    extract_star_sof_data.append(ff_norm_flat, "UTIL_MASTER_FLAT")
    extract_star_sof_data.append(master_dark_bpm, "CAL_DARK_BPM")
    extract_star_sof_data.write(extract_star_sof)

    commands += [
        f"python {join(python_dir, 'beta_pic_create_planet_model.py')} {date} {setting} A --model={ff_name_extracted_A_planet} --out={planet_img_A} --out_model={planet_model_A} \n"
        f"python {join(python_dir, 'diff_img.py')} {ff_name_science_A} {planet_img_A} --out={ff_name_science_A_subtracted}\n"
        f"python {join(python_dir, 'beta_pic_create_planet_model.py')} {date} {setting} B --model={ff_name_extracted_B_planet} --out={planet_img_B} --out_model={planet_model_B} \n"
        f"python {join(python_dir, 'diff_img.py')} {ff_name_science_B} {planet_img_B} --out={ff_name_science_B_subtracted}\n"
        f"esorex cr2res_obs_nodding --subtract_interorder_column=FALSE {extract_star_sof}\n"
        + move_obs_nodding_files(outdir, f"subtracted_{i}") +
        f"python {join(python_dir, 'beta_pic_show_difference.py')} {ff_name_extracted_AB_initial} {ff_name_extracted_AB_subtracted}\n"
        f"python {join(python_dir, 'diff_img.py')} {ff_name_science_A} {ff_name_model_A_subtracted} --out={ff_name_science_A_subtracted}\n"
        f"python {join(python_dir, 'diff_img.py')} {ff_name_science_B} {ff_name_model_B_subtracted} --out={ff_name_science_B_subtracted}\n"
        f"esorex cr2res_obs_nodding --subtract_interorder_column=FALSE {extract_star_sof}\n"
        + move_obs_nodding_files(outdir, f"planet_{i}") +
        f"python {join(python_dir, 'beta_pic_show_planet.py')} {ff_name_extracted_A_planet_new} {planet_model_A}\n"
    ]


script_fname = join(outdir, "beta_pic_run.sh")
with open(script_fname, "w") as f:
    f.writelines(commands)

# Make the script executable
cmd = f"chmod +x {script_fname}"
process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
output, error = process.communicate()
