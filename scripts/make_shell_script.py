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

def get_obs_nodding_filenames(outdir:str, suffix:str = "") -> list:
    """ Get the filenames of the cr2res_obs_nodding recipe (with an additional suffix) """
    if suffix != "":
        suffix = f"_{suffix}"

    files = []
    for nodding in "AB":
        files += [f'cr2res_obs_nodding_combined{nodding}{suffix}.fits']
        files += [f'cr2res_obs_nodding_extracted{nodding}{suffix}.fits']
        files += [f'cr2res_obs_nodding_slitfunc{nodding}{suffix}.fits']
        files += [f'cr2res_obs_nodding_model{nodding}{suffix}.fits']
        files += [f'cr2res_obs_nodding_trace_wave_{nodding}{suffix}.fits']
    files += [f'cr2res_obs_nodding_extracted_combined{suffix}.fits']
    files = [join(outdir, f) for f in files]

    return files

def move_obs_nodding_files(outdir:str, suffix:str, copy:bool=False) -> str:
    """ Create the commands to move (copy) the files from the cr2res_obs_nodding
        recipe to the same names with an additional suffix """
    ff_names = get_obs_nodding_filenames(outdir)
    ff_suffix = get_obs_nodding_filenames(outdir, suffix)
    cmd = "cp" if copy else "mv"
    commands = [f"{cmd} {n} {s}\n" for n, s in zip(ff_names, ff_suffix)]
    commands = "".join(commands)
    return commands

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

# Static filenames
tw_name = join(outdir, "cr2res_util_calib_flat_collapsed_tw.fits")
ff_norm_flat = join(outdir, "cr2res_util_normflat_Open_master_flat.fits")
master_dark_bpm = join(outdir, "cr2res_cal_dark_bpm.fits")
ff_name_science_A = join(outdir, 'cr2res_util_calib_science_A_collapsed.fits')
ff_name_science_B = join(outdir, 'cr2res_util_calib_science_B_collapsed.fits')
ff_name_extracted_AB_initial = join(outdir, 'cr2res_obs_nodding_extracted_combined_initial.fits')

# Iteration filenames for this iteration only
ff_name_model_A_subtracted = join(outdir, 'cr2res_obs_nodding_modelA_subtracted_0.fits')
ff_name_model_B_subtracted = join(outdir, 'cr2res_obs_nodding_modelB_subtracted_0.fits')
ff_name_extracted_AB_subtracted = join(outdir, 'cr2res_obs_nodding_extracted_combined_subtracted_0.fits')
ff_name_extracted_A_planet = join(outdir, f"cr2res_obs_nodding_extractedA_planet_0.fits")
ff_name_science_A_subtracted = join(outdir, 'cr2res_util_calib_science_A_collapsed_subtracted_0.fits')
ff_name_science_B_subtracted = join(outdir, 'cr2res_util_calib_science_B_collapsed_subtracted_0.fits')

planet_img_A = join(outdir, 'beta_pic_img_A_0.fits')
planet_model_A = join(outdir, 'beta_pic_model_A_0.fits')
planet_img_B = join(outdir, 'beta_pic_img_B_0.fits')
planet_model_B = join(outdir, 'beta_pic_model_B_0.fits')

# Plot file names
# ff_name_plot_difference = join(outdir, 'beta_pic_difference_i0_c{chip}_o{order}.png')
ff_name_plot_planet = join(outdir, 'beta_pic_planet.png')

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
python = "python"
esorex = "esorex"
create_model = join(python_dir, 'beta_pic_create_planet_model.py')
diff_img = join(python_dir, 'diff_img.py')
show_difference = join(python_dir, 'beta_pic_show_difference.py')
show_planet = join(python_dir, 'beta_pic_show_planet.py')

commands = [
    "#!/bin/bash\n"
    # Copy the initial files to a different name
    + move_obs_nodding_files(outdir, "initial") +
    # Subtract the initial model of the planet
    f"{python} {create_model} {date} {setting} A --out={planet_img_A} --out_model={planet_model_A}\n"
    f"{python} {diff_img} {ff_name_science_A} {planet_img_A} --out={ff_name_science_A_subtracted}\n"
    f"{python} {create_model} {date} {setting} B --out={planet_img_B} --out_model={planet_model_B}\n"
    f"{python} {diff_img} {ff_name_science_B} {planet_img_B} --out={ff_name_science_B_subtracted}\n"
    # Run the extraction of the spectrum
    f"{esorex} cr2res_obs_nodding --subtract_interorder_column=FALSE {extract_star_sof}\n"
    # and save the data to a new location
    + move_obs_nodding_files(outdir, "subtracted_0") +
    # Plot the results
    # f"{python} {show_difference} {ff_name_extracted_AB_initial} {ff_name_extracted_AB_subtracted} --out={ff_name_plot_difference}\n"
    # Remove the star model from the images
    f"{python} {diff_img} {ff_name_science_A} {ff_name_model_A_subtracted} --out={ff_name_science_A_subtracted}\n"
    f"{python} {diff_img} {ff_name_science_B} {ff_name_model_B_subtracted} --out={ff_name_science_B_subtracted}\n"
    # Run the extraction again with just the planet
    f"{esorex} cr2res_obs_nodding --subtract_interorder_column=FALSE {extract_star_sof}\n"
    # Move the files again
    + move_obs_nodding_files(outdir, "planet_0") +
    # Plot the extracted planet spectrum
    f"{python} {show_planet} {ff_name_extracted_A_planet} {planet_model_A} --out={ff_name_plot_planet}\n"
    
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
    ff_name_model_A_planet = join(outdir, f"cr2res_obs_nodding_modelA_planet_{i-1}.fits")
    ff_name_extracted_B_planet = join(outdir, f"cr2res_obs_nodding_extractedB_planet_{i-1}.fits")
    ff_name_model_B_planet = join(outdir, f"cr2res_obs_nodding_modelB_planet_{i-1}.fits")
    ff_name_extracted_AB_subtracted = join(outdir, f'cr2res_obs_nodding_extracted_combined_subtracted_{i}.fits')
    # ff_name_plot_difference = join(outdir, f'beta_pic_difference_i{i}_c{{chip}}_o{{order}}.png')

    extract_star_sof = join(outdir, f"subtracted_{i}.sof")
    extract_star_sof_data = SofFile()
    extract_star_sof_data.append(ff_name_science_A_subtracted, "OBS_NODDING_OTHER")
    extract_star_sof_data.append(ff_name_science_B_subtracted, "OBS_NODDING_OTHER")
    extract_star_sof_data.append(tw_name, "CAL_FLAT_TW")
    extract_star_sof_data.append(ff_norm_flat, "UTIL_MASTER_FLAT")
    extract_star_sof_data.append(master_dark_bpm, "CAL_DARK_BPM")
    extract_star_sof_data.write(extract_star_sof)

    commands += [
        # f"{python} {create_model} {date} {setting} A --model={ff_name_extracted_A_planet} --out={planet_img_A} --out_model={planet_model_A} \n"
        f"{python} {diff_img} {ff_name_science_A} {ff_name_model_A_planet} --out={ff_name_science_A_subtracted}\n"
        # f"{python} {create_model} {date} {setting} B --model={ff_name_extracted_B_planet} --out={planet_img_B} --out_model={planet_model_B} \n"
        f"{python} {diff_img} {ff_name_science_B} {ff_name_model_B_planet} --out={ff_name_science_B_subtracted}\n"
        f"{esorex} cr2res_obs_nodding --subtract_interorder_column=FALSE {extract_star_sof}\n"
        + move_obs_nodding_files(outdir, f"subtracted_{i}") +
        # f"{python} {show_difference} {ff_name_extracted_AB_initial} {ff_name_extracted_AB_subtracted} --out={ff_name_plot_difference}\n"
        f"{python} {diff_img} {ff_name_science_A} {ff_name_model_A_subtracted} --out={ff_name_science_A_subtracted}\n"
        f"{python} {diff_img} {ff_name_science_B} {ff_name_model_B_subtracted} --out={ff_name_science_B_subtracted}\n"
        f"{esorex} cr2res_obs_nodding --subtract_interorder_column=FALSE {extract_star_sof}\n"
        + move_obs_nodding_files(outdir, f"planet_{i}") +
        f"{python} {show_planet} {ff_name_extracted_A_planet_new} {planet_model_A} --out={ff_name_plot_planet}\n"
    ]


script_fname = join(outdir, "beta_pic_run.sh")
with open(script_fname, "w") as f:
    f.writelines(commands)

# Make the script executable
cmd = f"chmod +x {script_fname}"
process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
output, error = process.communicate()
