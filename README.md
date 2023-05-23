# BetaPic
Extract Beta Pic planet

To execute this package first create a shell script using:
python make_shell_script.py {date} {setting} {work_dir}

where:
    date is the date of the observation
    setting is the wavelength setting (e.g. L3262)
    work_dir is where the input and output files are located

This will create a beta_pic_run.sh script in work_dir.
Then you should run that shell script from within that directory.

This will call beta_pic_create_model.py to make a simulated model
of the planet spectrum.
This is then subtracted from the observation using diff_img.py.
The stellar spectrum is then extracted using the cr2res_obs_nodding recipe.
This stellar spectrum is then subtracted from the observation, so that only the
planet is left in the observations. These observations are then extracted again to 
get the planet spectrum. 

Using this planet spectrum we can start the process over again for a total of 5 times.

# FAQ
The model stellar spectrum is created using PySME and used to estimate the strength of 
initial planet model. Thus create_pysme_spectrum.py needs to be run for each setting.