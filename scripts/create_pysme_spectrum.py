from pysme.sme import SME_Structure
from pysme.linelist.vald import ValdFile
from pysme.synthesize import synthesize_spectrum
from pysme.gui.plot_plotly import FinalPlot

from astropy.io import fits
from os.path import dirname, join
import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation


date = "2022-11-29"
setting = "L3262"
# nodding = "A"

fname_wave = f"/scratch/ptah/anwe5599/CRIRES/{date}_{setting}/extr/cr2res_obs_nodding_extracted_combined.fits"
hdu_wave = fits.open(fname_wave)

wave_range = []
for chip in [1, 2, 3]:
    spectrum = hdu_wave[f"CHIP{chip}.INT1"].data
    orders = sorted([int(c[:2]) for c in spectrum.names if c[-4:] == "SPEC"])
    for order in orders:
        wave = spectrum[f"{order:02}_01_WL"]
        wmin, wmax = wave[0], wave[-1]
        wave_range += [[wmin, wmax]]

# Convert to Angstrom
wave_range = np.array(wave_range)
wave_range *= 10

sme = SME_Structure()
sme.vrad_flag = "fix"
sme.cscale_flag = "none"

sme.wran = wave_range

sme.linelist = ValdFile(join(dirname(__file__), "../data/beta_pic_L3262.lin"))

sme.teff = 8000
sme.logg = 3.83
sme.monh = -0.2

sme.vsini = 116

# Barycentric velocity correction
paranal = EarthLocation.of_site('Paranal')
sc = SkyCoord(ra=86.8211987087*u.deg, dec=-51.0665114264*u.deg)
barycorr = sc.radial_velocity_correction(obstime=Time(date), location=paranal)  
sme.vrad = barycorr.to_value(u.km/u.s)
# radial velocity of beta pic
sme.vrad += 20 

sme.normalize_by_continuum = False

fname = join(dirname(__file__), "../data/beta_pic.sme")
sme = synthesize_spectrum(sme)
sme.save(fname)

fig = FinalPlot(sme)
fig.save(filename="test.html")

wave = sme.wave.ravel()
synth = sme.synth.ravel()

# Convert to nm
wave = wave * 0.1

# Convert spectral flux density to match the 
# units for the planet spectrum
sme_unit = u.Unit("erg/s/cm2/AA/sr") * u.sr
# prt_unit = u.Unit("W/m^2/s")
nu_unit = u.Unit("erg/cm**2/s/Hz")

conversion = sme_unit.to(nu_unit, equivalencies=u.spectral_density(wave * u.nm))
synth = synth * conversion

data = np.asarray([wave, synth]).T
header = "Wave (nm),\tFlux (erg/cm**2/s/Hz)"
fname = join(dirname(__file__), f"../data/beta_pic_spec_{date}_{setting}.txt")
np.savetxt(fname, data, header=header)


pass
