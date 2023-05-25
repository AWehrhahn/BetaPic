# Beta Pic analysis
This package analyses observation of beta Pic b with CRIRES+. It exploits the angular seperation between the star and the planet of around 10 pixels. Here we start by subtracting an initial guess of the planet spectrum from the observations, then extract the stellar spectrum. This in turn is then subtracted from the observations to only extract the planet spectrum. By iterating this procedure a few times (3-5) the planet spectrum should converge to a consistent solution.

# Requirements
This package uses the CR2RES pipeline for the extraction of the stellar and planet spectrum (available here: https://www.eso.org/sci/software/pipelines/).
It also requires Python version 3.9 (others might work, but are untested) as well as the following packages:
astropy>=5.2.2
matplotlib>=3.3.4
numpy>=1.24.3
pandas>=1.1.5
scipy>=1.6.0
tqdm>=4.65.0
pyreduce-astro
pysme-astro
exoorbit

# How to use this package
The execution of this package is split into two parts. First it creates a shell script that will in turn call the other scripts in the correct order and with the right parameters. Additionally we need to create some data files before we can start the analysis.

## Creating the datafiles
### Planet model
Get the data file from Molliere, or run petitRADTRANS for beta Pic b.

### Stellar model
We use PySME to create a model of the stellar spectrum of beta Pic. For this analysis we will need a linelist for the L band, which can easily be created using VALD3 (http://vald.astro.uu.se/).

### Ephimeredes
We use the ephimeredes from Lacour et al. 2021 here to get the seperation in mas.

## Creating the shell script
To execute this package first create a shell script using:
`python make_shell_script.py {date} {setting} {work_dir}`

where:
    date is the date of the observation
    setting is the wavelength setting (e.g. L3262)
    work_dir is where the input and output files are located

This will create a beta_pic_run.sh script in work_dir.
Then you should run that shell script from within that directory.

## Running the script
Running the script will execute the steps of the analysis for a single dataset. The steps are described below:

### Creating a model of the planet observation
We here use a model of the planet spectrum created by Molliere using petitRADTRANS. Using the order traces and the slit illumination function of the stellar spectrum, we can then create a model observation that only contains the planet. The seperation between the star and planet is based on the ephimerides from Lacour et al. 2021. Here we also calibrate the stellar flux to W/m2/s using a model of the stellar spectrum created using PySME and fix the flux contrast between the star and planet to 7.7 magnitudes in the L band. This is done in the `beta_pic_create_planet_model.py` script.

![beta_pic_model](https://github.com/AWehrhahn/BetaPic/assets/31626864/a2fb7db5-a0b2-4cc1-8da4-1f07d5c0b629)

### Subtracting the initial model
The model is simply subtracted from the observation. Additionally we cut off negative values. This is done in the `diff_img.py` script.

### Spectrum extraction
The stellar spectra are extracted using the standard CR2RES instrument pipeline.

### Subtracting the star
The star is simply subtracted from the observation in the same way as the planet. The stellar model is one of the outputs from the CR2RES pipeline and is used directly.

### Extracting the planet
Again we use the normal CR2RES pipeline to extract the leftover planet signal. The pipeline will automatically adjust the slit illumination function to the postion of the planet.

### Subtracting the planet
The next iteration we can again use the model created by the pipeline to subtract only the planet.

### Repeat
Repeat this process until the planet spectrum converges



